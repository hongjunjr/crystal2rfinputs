"""
Crystal2HLT 结构加载、清理与组装模块

处理 PDB/CIF 文件的读取、清理、biological assembly 构建等
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure import AtomArray, get_chains

from ..config import (
    Config,
    ChainSequenceInfo,
    ResidueMapping,
    PROTEIN_RESIDUES,
)


# =============================================================================
# 结构加载
# =============================================================================

def detect_file_format(filepath: str) -> str:
    """检测文件格式 (pdb/cif)"""
    ext = Path(filepath).suffix.lower()
    if ext in ['.pdb', '.ent']:
        return 'pdb'
    elif ext in ['.cif', '.mmcif']:
        return 'cif'
    else:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('data_'):
                return 'cif'
            else:
                return 'pdb'


def load_structure(input_path: str, config: Config, logger: logging.Logger) -> AtomArray:
    """
    读取PDB/mmCIF，统一为Biotite AtomArray
    """
    logger.info(f"读取结构文件: {input_path}")
    
    file_format = detect_file_format(input_path)
    logger.info(f"检测到文件格式: {file_format}")
    
    if file_format == 'pdb':
        pdb_file = PDBFile.read(input_path)
        structure = pdb_file.get_structure(model=1)
    else:
        cif_file = CIFFile.read(input_path)
        structure = get_structure(cif_file, model=1)
    
    chains = get_chains(structure)
    unique_chains = np.unique(chains)
    
    logger.info(f"结构包含 {len(structure)} 个原子, {len(unique_chains)} 条链: {list(unique_chains)}")
    
    return structure


def filter_protein_only(structure: AtomArray, logger: logging.Logger) -> AtomArray:
    """只保留蛋白质残基"""
    protein_mask = np.isin(structure.res_name, PROTEIN_RESIDUES)
    filtered = structure[protein_mask]
    logger.info(f"过滤后保留 {len(filtered)} 个蛋白原子 (原 {len(structure)} 个)")
    return filtered


def clean_structure(structure: AtomArray, config: Config, logger: logging.Logger) -> AtomArray:
    """
    清理结构：
    - 去除水分子 (HOH)
    - 去除非蛋白HETATM
    - 处理altloc
    """
    logger.info("清理结构...")
    cleaned = filter_protein_only(structure, logger)
    return cleaned


# =============================================================================
# Biological Assembly
# =============================================================================

def extract_mmcif_assemblies(cif_file: CIFFile, logger: logging.Logger) -> dict:
    """
    从 mmCIF 提取 biological assembly 信息
    """
    from biotite.structure.io.pdbx import list_assemblies
    
    try:
        assemblies = list_assemblies(cif_file)
        logger.info(f"找到 {len(assemblies)} 个 biological assemblies")
        for assembly_id, description in assemblies.items():
            logger.debug(f"  Assembly {assembly_id}: {description}")
        return assemblies
    except Exception as e:
        logger.warning(f"无法提取 assembly 信息: {e}")
        return {}


def build_biological_assembly(
    cif_file: CIFFile,
    assembly_id: str = "1",
    logger: logging.Logger = None
) -> AtomArray:
    """
    使用 Biotite 的 get_assembly() 构建 biological assembly
    """
    from biotite.structure.io.pdbx import get_assembly
    
    logger.info(f"构建 biological assembly {assembly_id}...")
    
    try:
        assembly = get_assembly(cif_file, assembly_id=assembly_id, model=1)
        logger.info(f"  Assembly {assembly_id}: {len(assembly)} 个原子")
        return assembly
    except Exception as e:
        logger.error(f"构建 assembly {assembly_id} 失败: {e}")
        raise


# =============================================================================
# 链序列提取
# =============================================================================

AA_CODE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def extract_chain_sequences(structure: AtomArray) -> dict[str, str]:
    """从结构中提取每条链的序列（简单版本）"""
    chains = {}
    for chain_id in np.unique(structure.chain_id):
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        ca_mask = chain_atoms.atom_name == 'CA'
        ca_atoms = chain_atoms[ca_mask]
        
        seq = ""
        for atom in ca_atoms:
            res_name = atom.res_name
            if res_name in AA_CODE:
                seq += AA_CODE[res_name]
            else:
                seq += 'X'
        
        chains[chain_id] = seq
    
    return chains


def extract_chain_sequences_v2(
    structure: AtomArray,
    logger: logging.Logger
) -> dict[str, ChainSequenceInfo]:
    """
    从结构中提取每条链的序列信息
    返回 ChainSequenceInfo 对象，包含残基列表和过滤信息
    """
    chain_infos = {}
    
    for chain_id in np.unique(structure.chain_id):
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        ca_mask = chain_atoms.atom_name == 'CA'
        ca_atoms = chain_atoms[ca_mask]
        
        if len(ca_atoms) == 0:
            chain_infos[chain_id] = ChainSequenceInfo(
                chain_id=chain_id,
                sequence="",
                residue_list=[],
                is_protein=False,
                filter_reason="no CA atoms"
            )
            continue
        
        seq = ""
        residue_list = []
        
        for atom in ca_atoms:
            res_name = atom.res_name
            res_id = atom.res_id
            icode = getattr(atom, 'ins_code', '') or ''
            
            if res_name in AA_CODE:
                seq += AA_CODE[res_name]
            else:
                seq += 'X'
            
            residue_list.append((res_id, icode))
        
        chain_infos[chain_id] = ChainSequenceInfo(
            chain_id=chain_id,
            sequence=seq,
            residue_list=residue_list,
            is_protein=True,
            length=len(seq)
        )
    
    protein_count = sum(1 for c in chain_infos.values() if c.is_protein)
    logger.info(f"提取到 {len(chain_infos)} 条链序列 ({protein_count} 蛋白, {len(chain_infos) - protein_count} 非蛋白)")
    
    return chain_infos


# =============================================================================
# mmCIF Residue Mapping
# =============================================================================

def extract_residue_mapping_from_mmcif(
    cif_file: CIFFile,
    logger: logging.Logger
) -> dict[str, list[ResidueMapping]]:
    """
    从 mmCIF 的 pdbx_poly_seq_scheme 类别提取 residue 映射
    """
    logger.info("从 mmCIF 提取 residue 映射（pdbx_poly_seq_scheme）...")
    
    block = cif_file.block
    
    if 'pdbx_poly_seq_scheme' not in block:
        logger.warning("mmCIF 缺少 pdbx_poly_seq_scheme 类别")
        return {}
    
    scheme_category = block['pdbx_poly_seq_scheme']
    chain_mappings = {}
    
    asym_ids = list(scheme_category['asym_id'].as_array())
    pdb_strand_ids = list(scheme_category['pdb_strand_id'].as_array()) if 'pdb_strand_id' in scheme_category else asym_ids
    seq_ids = list(scheme_category['seq_id'].as_array())
    auth_seq_nums = list(scheme_category['auth_seq_num'].as_array()) if 'auth_seq_num' in scheme_category else seq_ids
    pdb_ins_codes = list(scheme_category['pdb_ins_code'].as_array()) if 'pdb_ins_code' in scheme_category else ['?'] * len(asym_ids)
    mon_ids = list(scheme_category['mon_id'].as_array())
    auth_mon_ids = list(scheme_category['auth_mon_id'].as_array()) if 'auth_mon_id' in scheme_category else mon_ids
    
    for i in range(len(asym_ids)):
        asym_id = asym_ids[i]
        pdb_strand_id = pdb_strand_ids[i] if i < len(pdb_strand_ids) else asym_id
        seq_id = seq_ids[i]
        auth_seq_num = auth_seq_nums[i] if i < len(auth_seq_nums) else seq_id
        pdb_ins_code = pdb_ins_codes[i] if i < len(pdb_ins_codes) else '?'
        mon_id = mon_ids[i]
        auth_mon_id = auth_mon_ids[i] if i < len(auth_mon_ids) else mon_id
        
        if pdb_strand_id == '?' or pdb_strand_id is None:
            pdb_strand_id = asym_id
        if auth_seq_num == '?' or auth_seq_num is None:
            auth_seq_num = seq_id
        if pdb_ins_code == '?' or pdb_ins_code is None:
            pdb_ins_code = ''
        
        chain_id = pdb_strand_id
        
        if chain_id not in chain_mappings:
            chain_mappings[chain_id] = []
        
        mapping = ResidueMapping(
            abs_index=0,
            auth_chain_id=chain_id,
            auth_seq_id=int(auth_seq_num) if auth_seq_num and auth_seq_num != '?' else 0,
            pdbx_ins_code=pdb_ins_code.strip() if isinstance(pdb_ins_code, str) else '',
            label_seq_id=int(seq_id) if seq_id and seq_id != '?' else 0,
            label_comp_id=auth_mon_id if auth_mon_id else mon_id,
            label_asym_id=asym_id if asym_id else chain_id
        )
        
        chain_mappings[chain_id].append(mapping)
    
    logger.info(f"  提取到 {len(chain_mappings)} 条链的映射信息")
    return chain_mappings


def export_residue_mapping_table(
    chain_mappings: dict[str, list[ResidueMapping]],
    output_path: str,
    logger: logging.Logger
):
    """导出 residue 映射表为 TSV"""
    logger.info(f"导出 residue 映射表: {output_path}")
    
    with open(output_path, 'w') as f:
        f.write("abs_index\tauth_chain\tauth_resnum\ticode\tlabel_seq_id\t"
                "res_name\tlabel_asym_id\tchothia_num\tchothia_icode\tcdr_type\n")
        
        for chain_id in sorted(chain_mappings.keys()):
            mappings = chain_mappings[chain_id]
            for m in mappings:
                f.write(f"{m.abs_index}\t{m.auth_chain_id}\t{m.auth_seq_id}\t"
                       f"{m.pdbx_ins_code}\t{m.label_seq_id}\t{m.label_comp_id}\t"
                       f"{m.label_asym_id}\t{m.chothia_num or ''}\t"
                       f"{m.chothia_icode or ''}\t{m.cdr_type or ''}\n")
