"""
Crystal2HLT HLT 导出模块

导出 framework_HLT.pdb 和 target_HLT.pdb
"""

from __future__ import annotations

import os
import logging

import numpy as np
from biotite.structure import AtomArray, residue_iter
from biotite.structure.io.pdb import PDBFile

from ..config import (
    Config,
    RoleAssignment,
    CDRBoundary,
    ResidueTag,
    CHOTHIA_CDR_RANGES,
    PROTEIN_RESIDUES,
)
from ..analysis.anarci import get_cdr_boundaries


AA_CODE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def export_framework_hlt(
    structure: AtomArray,
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger
) -> tuple[str, list[CDRBoundary]]:
    """
    导出 framework_HLT.pdb
    - 重命名: heavy→H, light→L
    - 按H→L顺序写出
    - 残基重新编号为全局1-index
    """
    output_path = os.path.join(config.output_dir, "framework_HLT.pdb")
    logger.info(f"导出 framework_HLT.pdb...")
    
    # 确定链
    chains_to_export = []
    if assignment.heavy_chain:
        chains_to_export.append(('H', assignment.heavy_chain))
    if assignment.light_chain:
        chains_to_export.append(('L', assignment.light_chain))
    
    all_atoms = []
    all_cdr_boundaries = []
    global_res_idx = 0
    
    for new_chain_id, orig_chain_id in chains_to_export:
        chain_mask = structure.chain_id == orig_chain_id
        chain_atoms = structure[chain_mask]
        
        # 获取 ANARCI 结果
        anarci_result = assignment.anarci_results.get(orig_chain_id)
        
        # scFv 处理：只取对应域
        if assignment.mode == 'scfv' and orig_chain_id in assignment.anarci_results:
            # scFv: 需要根据 domain 边界过滤原子
            pass  # 简化版本，完整实现在原始文件中
        
        # 重新编号并收集原子
        visited_residues = set()
        for atom in chain_atoms:
            res_key = (atom.res_id, getattr(atom, 'ins_code', ''))
            if res_key not in visited_residues:
                visited_residues.add(res_key)
                global_res_idx += 1
            
            new_atom = atom.copy()
            new_atom.chain_id = new_chain_id
            new_atom.res_id = global_res_idx
            all_atoms.append(new_atom)
        
        # 提取 CDR 边界
        if anarci_result:
            chain_type = anarci_result.get('chain_type', 'H')
            if chain_type in ['K', 'L']:
                chain_type = 'L'
            
            cdr_ranges = CHOTHIA_CDR_RANGES.get(chain_type, {})
            for cdr_name, (start, end) in cdr_ranges.items():
                all_cdr_boundaries.append(CDRBoundary(
                    cdr_name=cdr_name,
                    chothia_start=start,
                    chothia_end=end,
                    abs_start=0,  # 简化
                    abs_end=0
                ))
    
    if not all_atoms:
        logger.error("没有找到要导出的原子")
        return output_path, []
    
    # 构建 AtomArray 并写入
    from biotite.structure import AtomArray as BA, array as ba_array
    new_structure = ba_array(all_atoms)
    
    pdb_file = PDBFile()
    pdb_file.set_structure(new_structure)
    pdb_file.write(output_path)
    
    logger.info(f"  写入 {len(new_structure)} 个原子到 {output_path}")
    
    return output_path, all_cdr_boundaries


def export_target_hlt(
    structure: AtomArray,
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger,
    resnum_conflict_info: dict | None = None
) -> tuple[str, dict]:
    """
    导出 target_HLT.pdb
    - 所有抗原链重命名为T
    - 如果有冲突，自动重编号
    """
    output_path = os.path.join(config.output_dir, "target_HLT.pdb")
    logger.info("导出 target_HLT.pdb...")
    
    resnum_mapping = {}
    
    # 过滤抗原链
    antigen_mask = np.isin(structure.chain_id, assignment.antigen_chains)
    antigen_atoms = structure[antigen_mask]
    
    if len(antigen_atoms) == 0:
        logger.warning("没有抗原原子")
        return output_path, resnum_mapping
    
    # 检测是否需要重编号
    need_renumber = resnum_conflict_info and resnum_conflict_info.get('has_conflict', False)
    
    all_atoms = []
    new_resnum = 0
    current_residue = None
    
    for atom in antigen_atoms:
        res_key = (atom.chain_id, atom.res_id, getattr(atom, 'ins_code', ''))
        
        if res_key != current_residue:
            current_residue = res_key
            new_resnum += 1
            resnum_mapping[res_key] = new_resnum
        
        new_atom = atom.copy()
        new_atom.chain_id = 'T'
        if need_renumber:
            new_atom.res_id = new_resnum
        else:
            new_atom.res_id = atom.res_id
        all_atoms.append(new_atom)
    
    if not all_atoms:
        return output_path, resnum_mapping
    
    from biotite.structure import array as ba_array
    new_structure = ba_array(all_atoms)
    
    pdb_file = PDBFile()
    pdb_file.set_structure(new_structure)
    pdb_file.write(output_path)
    
    logger.info(f"  写入 {len(new_structure)} 个原子到 {output_path}")
    
    return output_path, resnum_mapping


def write_sequences(
    structure: AtomArray,
    assignment: RoleAssignment,
    output_dir: str,
    logger: logging.Logger
):
    """写入序列FASTA文件"""
    logger.info("  写入序列文件")
    
    seq_dir = os.path.join(output_dir, "sequences")
    os.makedirs(seq_dir, exist_ok=True)
    
    def extract_seq(chain_id):
        mask = structure.chain_id == chain_id
        atoms = structure[mask]
        ca_mask = atoms.atom_name == 'CA'
        ca_atoms = atoms[ca_mask]
        seq = ""
        for atom in ca_atoms:
            seq += AA_CODE.get(atom.res_name, 'X')
        return seq
    
    if assignment.heavy_chain:
        seq = extract_seq(assignment.heavy_chain)
        with open(os.path.join(seq_dir, "H.fasta"), 'w') as f:
            f.write(f">H\n{seq}\n")
    
    if assignment.light_chain:
        seq = extract_seq(assignment.light_chain)
        with open(os.path.join(seq_dir, "L.fasta"), 'w') as f:
            f.write(f">L\n{seq}\n")
    
    if assignment.antigen_chains:
        seqs = [extract_seq(c) for c in assignment.antigen_chains]
        with open(os.path.join(seq_dir, "T.fasta"), 'w') as f:
            for i, seq in enumerate(seqs):
                f.write(f">T_{i+1}\n{seq}\n")
