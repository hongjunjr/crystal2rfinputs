#!/usr/bin/env python3
"""
Crystal2HLT-RFInputs: 共晶结构 → RFantibody HLT 输入自动构建管线

功能：
- 输入：RCSB下载的 complex.pdb 或 complex.cif
- 输出：
  - target_HLT.pdb: 抗原结构（链统一命名为T）
  - framework_HLT.pdb: 抗体框架（Heavy=H, Light=L）+ REMARK CDR标注
  - epitope_patch.txt: 抗原界面残基
  - hotspots_rfantibody.txt: ppi.hotspot_res=[T305,T456,...]
  - design_loops_rfantibody.txt: antibody.design_loops=[...]
  - qc_report.json / qc_report.md: 质量控制报告
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import logging
import hashlib
from dataclasses import dataclass, field, asdict
from typing import Literal, Optional, Any
from pathlib import Path
from datetime import datetime

import numpy as np

# Biotite for structure handling
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure import AtomArray, AtomArrayStack, array, residue_iter, get_chains
from biotite.structure import filter_amino_acids


# For neighbor search
from scipy.spatial import KDTree


# =============================================================================
# 数据结构定义
# =============================================================================

PROTEIN_RESIDUES = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]

# Chothia CDR定义
CHOTHIA_CDR_RANGES = {
    'H': {
        'H1': (26, 32),   # Heavy chain CDR1
        'H2': (52, 56),   # Heavy chain CDR2
        'H3': (95, 102),  # Heavy chain CDR3
    },
    'L': {
        'L1': (24, 34),   # Light chain CDR1
        'L2': (50, 56),   # Light chain CDR2
        'L3': (89, 97),   # Light chain CDR3
    },
}

# Fv region 范围 (用于裁剪)
CHOTHIA_FV_RANGES = {
    'H': (1, 113),  # Heavy chain Fv
    'L': (1, 107),  # Light chain Fv
}


@dataclass
class ResidueTag:
    """残基唯一标识（含插入码）"""
    chain_id: str
    res_num: int
    icode: str = ""
    res_name: str = ""
    
    def to_pdb_format(self) -> str:
        """T:305A"""
        return f"{self.chain_id}:{self.res_num}{self.icode}"
    
    def to_rfantibody_format(self) -> str:
        """用于 ppi.hotspot_res: T305 (不含icode，RFantibody格式)"""
        return f"{self.chain_id}{self.res_num}"
    
    def __hash__(self):
        return hash((self.chain_id, self.res_num, self.icode))
    
    def __eq__(self, other):
        return (self.chain_id, self.res_num, self.icode) == (other.chain_id, other.res_num, other.icode)


@dataclass
class AntibodyPairing:
    """单个抗体配对信息 (P1.3)"""
    heavy_chain: str               # VH 链 ID
    light_chain: str | None        # VL 链 ID，nanobody 为 None
    contact_score: float = 0.0     # VH-VL 接触分数
    index: int = 0                 # 配对索引 (1, 2, 3...)
    mode: Literal["antibody", "nanobody"] = "antibody"


@dataclass
class RoleAssignment:
    """链角色分配结果"""
    heavy_chain: str | None        # 原始chain_id (主抗体)
    light_chain: str | None        # nanobody时为None
    antigen_chains: list[str]      # 可多条
    mode: Literal["antibody", "nanobody"]
    vh_vl_contact_score: float = 0.0
    all_heavy_chains: list[str] = field(default_factory=list)
    all_light_chains: list[str] = field(default_factory=list)
    anarci_results: dict = field(default_factory=dict)
    # P1.3: 所有抗体配对
    all_pairings: list[AntibodyPairing] = field(default_factory=list)


@dataclass
class CDRBoundary:
    """CDR边界定义"""
    cdr_name: str                  # H1/H2/H3/L1/L2/L3
    chothia_start: int             # Chothia编号起点
    chothia_end: int               # Chothia编号终点
    abs_start: int                 # HLT全局索引起点
    abs_end: int                   # HLT全局索引终点
    length: int = 0
    
    def __post_init__(self):
        if self.length == 0:
            self.length = self.abs_end - self.abs_start + 1


@dataclass
class HotspotResidue:
    """热点残基评分"""
    tag: ResidueTag
    score: float
    contact_count: int = 0
    min_distance: float = 0.0
    avg_cb_dist_to_cdr5: float | None = None
    
    def to_dict(self) -> dict:
        return {
            'residue': self.tag.to_pdb_format(),
            'res_name': self.tag.res_name,
            'score': round(self.score, 4),
            'contact_count': self.contact_count,
            'min_distance': round(self.min_distance, 4),
            'avg_cb_dist_to_cdr5': round(self.avg_cb_dist_to_cdr5, 4) if self.avg_cb_dist_to_cdr5 else None
        }


@dataclass
class ResidueMapping:
    """
    完整的 residue 映射（基于 mmCIF pdbx_poly_seq_scheme）
    用于追溯从原始编号到输出编号的映射关系
    """
    abs_index: int              # 输出文件中的全局索引（1-indexed）
    auth_chain_id: str          # author 指定的链ID
    auth_seq_id: int            # author 序列编号
    pdbx_ins_code: str          # 插入码
    label_seq_id: int           # label 序列编号
    label_comp_id: str          # 残基名称（3字母码）
    label_asym_id: str          # label asymmetric unit ID
    chothia_num: int | None = None      # Chothia 编号（仅抗体链）
    chothia_icode: str | None = None    # Chothia 插入码
    cdr_type: str | None = None         # H1/H2/H3/L1/L2/L3（仅CDR）


# =============================================================================
# V2 新增数据结构
# =============================================================================

@dataclass
class FastaRecord:
    """FASTA 记录（带解析后的 metadata）"""
    id: str                     # FASTA ID (如 "1VFB_1")
    description: str            # 完整描述
    sequence: str               # 序列
    chain_ids: list[str]        # 提取的链ID (如 ["A"])
    role_hints: list[str]       # 角色提示 (如 ["heavy_chain"])
    entity_hints: list[str]     # entity 提示
    confidence: float = 0.0     # 解析置信度
    raw_header: str = ""        # 原始 header


@dataclass
class ChainSequenceInfo:
    """结构链序列信息"""
    chain_id: str
    sequence: str
    residue_list: list[tuple[int, str]] = field(default_factory=list)  # [(resnum, icode), ...]
    is_protein: bool = True
    length: int = 0
    filter_reason: str = ""     # 如果被过滤，原因
    
    def __post_init__(self):
        if self.length == 0:
            self.length = len(self.sequence)


@dataclass
class ChainMapping:
    """结构链 ↔ FASTA 映射结果"""
    struct_chain_id: str
    fasta_id: str | None
    identity: float
    coverage_struct: float
    coverage_fasta: float
    is_ambiguous: bool = False
    confidence: float = 0.0


@dataclass
class ChainRoleInfo:
    """链角色判定结果"""
    chain_id: str
    is_antibody: bool
    domain_type: str | None = None      # VH/VL/VHH/VNAR/None
    domain_range: tuple[int, int] | None = None  # 域在链上的残基范围
    anarci_score: float = 0.0
    fasta_hint: str | None = None
    final_confidence: float = 0.0
    is_crystal_packing: bool = False


@dataclass
class Config:
    """全局配置"""
    input_file: str
    outdir: str
    pdbid: str = ""
    mode: Literal["antibody", "nanobody", "auto"] = "auto"
    numbering: Literal["chothia", "imgt", "aho"] = "chothia"
    interface_cutoff: float = 4.5
    hotspots_k: int = 12
    hotspots_method: Literal["contact_score", "rfantibody_cbeta5"] = "contact_score"
    space_dedup_angstrom: float = 6.0
    use_bioassembly: Literal["auto", "true", "false"] = "auto"
    bioassembly_id: int = 1
    remove_waters: bool = True
    remove_hetero: bool = True
    keep_glycans: bool = False
    target_contact_only: bool = False
    crop_target: bool = False
    crop_margin_angstrom: float = 10.0
    h_crop: int = 115
    l_crop: int = 110
    use_anarci: bool = True
    fail_on_no_anarci: bool = False  # 改为False，不强制需要ANARCI
    verbose: bool = False
    
    # 手动链ID指定（优先于自动识别）
    heavy_chain_id: str = ""   # 手动指定重链ID
    light_chain_id: str = ""   # 手动指定轻链ID
    antigen_chain_ids: str = ""  # 手动指定抗原链ID（逗号分隔）
    
    # V2 新增: FASTA 输入支持
    fasta_file: str = ""       # FASTA 文件路径
    fasta_mode: Literal["auto", "strict", "off"] = "auto"  # FASTA 使用模式
    target_chain_ids: str = ""  # 手动指定 target 链ID（逗号分隔）
    prefer_cif_assembly: Literal["auto", "asu", "assembly"] = "auto"  # CIF assembly 偏好
    multi_antibody: Literal["off", "top1", "topN"] = "top1"  # 多抗体处理模式
    
    # 内部使用
    work_dir: str = ""
    output_dir: str = ""
    input_dir: str = ""
    reports_dir: str = ""



@dataclass
class QCReport:
    """质量控制报告 (V2 增强)"""
    meta: dict = field(default_factory=dict)
    structure_info: dict = field(default_factory=dict)
    chain_assignment: dict = field(default_factory=dict)
    cdr_info: dict = field(default_factory=dict)
    interface_info: dict = field(default_factory=dict)
    hotspots_info: dict = field(default_factory=dict)
    hlt_validation: dict = field(default_factory=dict)
    rfantibody_snippets: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    success: bool = True
    
    # V2 新增字段
    fasta_info: dict = field(default_factory=dict)  # FASTA 解析信息
    chain_mapping_info: dict = field(default_factory=dict)  # 链映射信息
    chain_role_table: list = field(default_factory=list)  # 链角色表
    target_selection_info: dict = field(default_factory=dict)  # target 选择依据
    resnum_conflict_info: dict = field(default_factory=dict)  # resnum 冲突信息


# =============================================================================
# 日志设置
# =============================================================================

def setup_logging(verbose: bool = False, log_file: str | None = None):
    """配置日志"""
    level = logging.DEBUG if verbose else logging.INFO
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)


# =============================================================================
# Module FASTA: FASTA 解析与角色提取 (V2 新增)
# =============================================================================

def parse_fasta(fasta_path: str, logger: logging.Logger) -> list[FastaRecord]:
    """
    解析 PDB FASTA 文件，返回 FastaRecord 列表
    使用 Biopython SeqIO
    """
    try:
        from Bio import SeqIO
    except ImportError:
        logger.error("Biopython 未安装。请运行: pip install biopython")
        raise ImportError("Biopython is required for FASTA parsing. Install with: pip install biopython")
    
    records = []
    logger.info(f"解析 FASTA 文件: {fasta_path}")
    
    for rec in SeqIO.parse(fasta_path, "fasta"):
        chain_ids = extract_chain_ids_from_header(rec.description)
        role_hints = extract_role_hints_from_header(rec.description)
        entity_hints = extract_entity_hints_from_header(rec.description)
        
        fasta_rec = FastaRecord(
            id=rec.id,
            description=rec.description,
            sequence=str(rec.seq),
            chain_ids=chain_ids,
            role_hints=role_hints,
            entity_hints=entity_hints,
            raw_header=rec.description
        )
        # 计算置信度
        fasta_rec.confidence = compute_fasta_parse_confidence(fasta_rec)
        records.append(fasta_rec)
        
        logger.debug(f"  FASTA 记录: {rec.id}, 链={chain_ids}, 角色={role_hints}")
    
    logger.info(f"  解析到 {len(records)} 条 FASTA 记录")
    return records


def extract_chain_ids_from_header(header: str) -> list[str]:
    """
    从 FASTA header 提取链ID
    支持格式:
    - "Chain A" / "Chains A, B"
    - "[auth H]" / "[auth L]"
    - "pdb_strand_id=A"
    
    示例:
    - "1VFB_1|Chain A|IGG1-KAPPA..." -> ["A"]
    - "7ZXK_3|Chains E[auth H], G[auth J]|..." -> ["H", "J"]
    """
    chain_ids = []
    
    # 模式1: "Chain A" 或 "Chains A, B, C"
    match = re.search(r'Chains?\s+([A-Za-z0-9,\s\[\]auth]+?)(?:\||$)', header)
    if match:
        chains_str = match.group(1)
        # 处理 "E[auth H], G[auth J]" 格式 - 优先使用 auth 编号
        auth_matches = re.findall(r'\[auth\s+([A-Za-z0-9]+)\]', chains_str)
        if auth_matches:
            chain_ids.extend(auth_matches)
        else:
            # 简单格式 "A, B"
            simple_chains = re.findall(r'([A-Za-z0-9])(?:\s*,|\s*$|(?=\[))', chains_str)
            if simple_chains:
                chain_ids.extend(simple_chains)
            else:
                # 最简单情况 "Chain A"
                single_match = re.search(r'Chain\s+([A-Za-z0-9])', header)
                if single_match:
                    chain_ids.append(single_match.group(1))
    
    return list(set(chain_ids))


def extract_role_hints_from_header(header: str) -> list[str]:
    """
    提取角色关键词: heavy, light, vhh, vnar, scfv, fab, nanobody
    返回标准化的角色提示列表
    """
    header_lower = header.lower()
    role_keywords = [
        ('heavy chain', 'heavy_chain'),
        ('light chain', 'light_chain'),
        ('vhh', 'vhh'),
        ('vnar', 'vnar'),
        ('scfv', 'scfv'),
        ('sc-fv', 'scfv'),
        ('single chain fv', 'scfv'),
        ('nanobody', 'nanobody'),
        ('fab fragment', 'fab'),
        ('fab', 'fab'),
        ('fv ', 'fv'),
        ('(fv)', 'fv'),
    ]
    
    hints = []
    for pattern, hint in role_keywords:
        if pattern in header_lower:
            hints.append(hint)
    
    return list(set(hints))


def extract_entity_hints_from_header(header: str) -> list[str]:
    """
    提取 entity/polymer/分子描述
    从 "|...|" 格式中提取描述文本
    """
    parts = header.split('|')
    hints = []
    
    if len(parts) >= 3:
        # 取第三部分作为分子描述
        description = parts[2].strip()
        if description:
            hints.append(description)
    
    return hints


def compute_fasta_parse_confidence(record: FastaRecord) -> float:
    """
    计算 FASTA 解析置信度
    基于: 是否提取到链ID、是否有角色提示、序列长度
    """
    score = 0.0
    
    # 有链ID
    if record.chain_ids:
        score += 0.4
    
    # 有角色提示
    if record.role_hints:
        score += 0.3
    
    # 有分子描述
    if record.entity_hints:
        score += 0.2
    
    # 序列长度合理（蛋白质）
    if 50 < len(record.sequence) < 1000:
        score += 0.1
    
    return min(score, 1.0)


def compute_sequence_alignment(seq1: str, seq2: str) -> tuple[float, float, float]:
    """
    计算两条序列的比对分数
    返回: (identity, coverage_seq1, coverage_seq2)
    使用 difflib 进行简单比对
    """
    from difflib import SequenceMatcher
    
    if not seq1 or not seq2:
        return 0.0, 0.0, 0.0
    
    matcher = SequenceMatcher(None, seq1, seq2)
    identity = matcher.ratio()
    
    # 计算覆盖度
    matching_blocks = matcher.get_matching_blocks()
    total_match = sum(block.size for block in matching_blocks)
    
    cov1 = total_match / len(seq1) if seq1 else 0
    cov2 = total_match / len(seq2) if seq2 else 0
    
    return identity, cov1, cov2


def map_structure_chains_to_fasta(
    chain_seqs: dict[str, ChainSequenceInfo],
    fasta_records: list[FastaRecord],
    mode: str,  # auto/strict
    logger: logging.Logger
) -> dict[str, ChainMapping]:
    """
    将结构链映射到 FASTA 记录
    使用两阶段策略: 快速筛选 + 序列比对
    """
    logger.info(f"映射结构链到 FASTA (mode={mode})...")
    mappings = {}
    
    for chain_id, chain_info in chain_seqs.items():
        if not chain_info.is_protein:
            logger.debug(f"  跳过非蛋白链: {chain_id}")
            continue
        
        # 阶段1: 快速候选筛选
        candidates = []
        for fasta_rec in fasta_records:
            priority = 0.0
            
            # 链ID 直接匹配（最高优先级）
            if chain_id in fasta_rec.chain_ids:
                priority = 1.0
            # 长度相近
            elif abs(len(fasta_rec.sequence) - chain_info.length) < 30:
                priority = 0.5
            
            if priority > 0:
                candidates.append((fasta_rec, priority))
        
        # 阶段2: 序列比对
        best_match = None
        best_score = 0.0
        is_ambiguous = False
        second_best_score = 0.0
        
        for fasta_rec, priority in candidates:
            identity, cov_s, cov_f = compute_sequence_alignment(
                chain_info.sequence, fasta_rec.sequence
            )
            score = identity * min(cov_s, cov_f) * (0.5 + 0.5 * priority)
            
            if score > best_score:
                second_best_score = best_score
                best_score = score
                best_match = ChainMapping(
                    struct_chain_id=chain_id,
                    fasta_id=fasta_rec.id,
                    identity=identity,
                    coverage_struct=cov_s,
                    coverage_fasta=cov_f,
                    is_ambiguous=False,
                    confidence=score
                )
            elif score > second_best_score:
                second_best_score = score
        
        # 检查歧义
        if best_match and second_best_score > 0 and best_score - second_best_score < 0.1:
            best_match.is_ambiguous = True
            logger.warning(f"  链 {chain_id} 映射存在歧义 (分数差 < 0.1)")
        
        if best_match and best_match.identity > 0.8:
            mappings[chain_id] = best_match
            logger.info(f"  链 {chain_id} -> {best_match.fasta_id} "
                       f"(identity={best_match.identity:.2f}, cov={best_match.coverage_struct:.2f})")
        elif mode == 'strict':
            raise ValueError(f"链 {chain_id} 无法映射到 FASTA (strict 模式要求所有链必须映射)")
        else:
            logger.warning(f"  链 {chain_id} 未找到高质量 FASTA 匹配")
            mappings[chain_id] = ChainMapping(
                struct_chain_id=chain_id,
                fasta_id=None,
                identity=0, coverage_struct=0, coverage_fasta=0,
                is_ambiguous=False, confidence=0
            )
    
    logger.info(f"  完成映射: {sum(1 for m in mappings.values() if m.fasta_id)} / {len(mappings)} 链成功映射")
    return mappings


def get_fasta_hint_for_chain(
    chain_id: str,
    chain_mappings: dict[str, ChainMapping],
    fasta_records: list[FastaRecord]
) -> str | None:
    """
    获取链对应的 FASTA 角色提示
    """
    if not chain_mappings or chain_id not in chain_mappings:
        return None
    
    mapping = chain_mappings[chain_id]
    if not mapping.fasta_id:
        return None
    
    for frec in fasta_records:
        if frec.id == mapping.fasta_id:
            if frec.role_hints:
                return ','.join(frec.role_hints)
    
    return None


# =============================================================================
# Module 0: 读取与规范化
# =============================================================================

def detect_file_format(filepath: str) -> str:
    """检测文件格式 (pdb/cif)"""
    ext = Path(filepath).suffix.lower()
    if ext in ['.pdb', '.ent']:
        return 'pdb'
    elif ext in ['.cif', '.mmcif']:
        return 'cif'
    else:
        # 尝试从内容判断
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('data_'):
                return 'cif'
            else:
                return 'pdb'


def normalize_structure(input_path: str, config: Config, logger: logging.Logger) -> AtomArray:
    """
    读取PDB/mmCIF，统一为Biotite AtomArray
    """
    logger.info(f"读取结构文件: {input_path}")
    
    file_format = detect_file_format(input_path)
    logger.info(f"检测到文件格式: {file_format}")
    
    if file_format == 'pdb':
        pdb_file = PDBFile.read(input_path)
        structure = pdb_file.get_structure(model=1)
    else:  # cif
        cif_file = CIFFile.read(input_path)
        structure = get_structure(cif_file, model=1)
    
    # 统计信息
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


# =============================================================================
# Module 0.5: mmCIF Metadata 与 Biological Assembly
# =============================================================================

def extract_mmcif_assemblies(cif_file: CIFFile, logger: logging.Logger) -> dict:
    """
    从 mmCIF 提取 biological assembly 信息
    使用 Biotite 的 list_assemblies()
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


def build_bioassembly(
    cif_file: CIFFile,
    assembly_id: str = "1",
    logger: logging.Logger = None
) -> AtomArray:
    """
    使用 Biotite 的 get_assembly() 构建 biological assembly
    
    参考 Biotite API:
    - biotite.structure.io.pdbx.get_assembly()
    - 使用 pdbx_struct_assembly_gen, pdbx_struct_oper_list, atom_site
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


def extract_residue_mapping_from_mmcif(
    cif_file: CIFFile,
    logger: logging.Logger
) -> dict[str, list[ResidueMapping]]:
    """
    从 mmCIF 的 pdbx_poly_seq_scheme 类别提取 residue 映射
    
    pdbx_poly_seq_scheme 包含:
    - asym_id: label asymmetric unit ID
    - entity_id: entity ID
    - seq_id: label sequence ID
    - mon_id: residue name
    - pdb_strand_id: author chain ID
    - auth_seq_num: author sequence number
    - pdb_ins_code: insertion code
    - auth_mon_id: author residue name
    """
    logger.info("从 mmCIF 提取 residue 映射（pdbx_poly_seq_scheme）...")
    
    # 访问 mmCIF block
    block = cif_file.block
    
    if 'pdbx_poly_seq_scheme' not in block:
        logger.warning("mmCIF 缺少 pdbx_poly_seq_scheme 类别，将使用结构数据构建简化映射")
        return {}
    
    scheme_category = block['pdbx_poly_seq_scheme']
    
    # 按链组织映射
    chain_mappings = {}
    
    # 获取各列数据（CIFCategory 像字典，用键访问后得到 CIFColumn，可转为数组）
    asym_ids = list(scheme_category['asym_id'].as_array())
    pdb_strand_ids = list(scheme_category['pdb_strand_id'].as_array()) if 'pdb_strand_id' in scheme_category else asym_ids
    seq_ids = list(scheme_category['seq_id'].as_array())
    auth_seq_nums = list(scheme_category['auth_seq_num'].as_array()) if 'auth_seq_num' in scheme_category else seq_ids
    pdb_ins_codes = list(scheme_category['pdb_ins_code'].as_array()) if 'pdb_ins_code' in scheme_category else ['?'] * len(asym_ids)
    mon_ids = list(scheme_category['mon_id'].as_array())
    auth_mon_ids = list(scheme_category['auth_mon_id'].as_array()) if 'auth_mon_id' in scheme_category else mon_ids
    
    # 遍历每一行
    for i in range(len(asym_ids)):
        asym_id = asym_ids[i] if i < len(asym_ids) else None
        pdb_strand_id = pdb_strand_ids[i] if i < len(pdb_strand_ids) else asym_id
        seq_id = seq_ids[i] if i < len(seq_ids) else None
        auth_seq_num = auth_seq_nums[i] if i < len(auth_seq_nums) else seq_id
        pdb_ins_code = pdb_ins_codes[i] if i < len(pdb_ins_codes) else '?'
        mon_id = mon_ids[i] if i < len(mon_ids) else None
        auth_mon_id = auth_mon_ids[i] if i < len(auth_mon_ids) else mon_id
        
        # 处理缺失值
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
            abs_index=0,  # 稍后填充
            auth_chain_id=chain_id,
            auth_seq_id=int(auth_seq_num) if auth_seq_num and auth_seq_num != '?' else 0,
            pdbx_ins_code=pdb_ins_code.strip() if isinstance(pdb_ins_code, str) else '',
            label_seq_id=int(seq_id) if seq_id and seq_id != '?' else 0,
            label_comp_id=auth_mon_id if auth_mon_id else mon_id,
            label_asym_id=asym_id if asym_id else chain_id
        )
        
        chain_mappings[chain_id].append(mapping)
    
    logger.info(f"  提取到 {len(chain_mappings)} 条链的映射信息")
    for chain_id, mappings in chain_mappings.items():
        logger.debug(f"    链 {chain_id}: {len(mappings)} 个残基")
    
    return chain_mappings


def integrate_anarci_to_mapping(
    chain_mappings: dict[str, list[ResidueMapping]],
    anarci_results: dict,
    cdr_ranges: dict,
    logger: logging.Logger
) -> dict[str, list[ResidueMapping]]:
    """
    将 ANARCI 的 Chothia 编号整合到 residue 映射表中
    """
    logger.info("整合 ANARCI Chothia 编号到映射表...")
    
    for chain_id, anarci_result in anarci_results.items():
        if chain_id not in chain_mappings:
            logger.warning(f"  链 {chain_id} 在映射表中不存在，跳过")
            continue
        
        numbering = anarci_result['numbering']  # [(pos, aa), ...]
        chain_type = anarci_result['chain_type']  # 'H' or 'L'
        
        mappings = chain_mappings[chain_id]
        
        # ANARCI numbering 与 structure 残基的对应
        seq_idx = 0
        for pos, aa in numbering:
            if aa == '-':  # gap
                continue
            
            if seq_idx >= len(mappings):
                break
            
            chothia_num, chothia_icode = pos
            mappings[seq_idx].chothia_num = chothia_num
            mappings[seq_idx].chothia_icode = chothia_icode if chothia_icode else ''
            
            # 判断是否在 CDR 中
            for cdr_name, (start, end) in cdr_ranges.get(chain_type, {}).items():
                if start <= chothia_num <= end:
                    mappings[seq_idx].cdr_type = cdr_name
                    break
            
            seq_idx += 1
        
        logger.debug(f"  链 {chain_id} ({chain_type}): 整合了 {seq_idx} 个 Chothia 编号")
    
    return chain_mappings


def export_residue_mapping_table(
    chain_mappings: dict[str, list[ResidueMapping]],
    output_path: str,
    logger: logging.Logger
):
    """导出 residue 映射表为 TSV（用于审计追溯）"""
    logger.info(f"导出 residue 映射表: {output_path}")
    
    with open(output_path, 'w') as f:
        # Header
        f.write("abs_index\tauth_chain\tauth_resnum\ticode\tlabel_seq_id\t"
                "res_name\tlabel_asym_id\tchothia_num\tchothia_icode\tcdr_type\n")
        
        for chain_id in sorted(chain_mappings.keys()):
            mappings = chain_mappings[chain_id]
            for m in mappings:
                f.write(f"{m.abs_index}\t{m.auth_chain_id}\t{m.auth_seq_id}\t"
                       f"{m.pdbx_ins_code}\t{m.label_seq_id}\t{m.label_comp_id}\t"
                       f"{m.label_asym_id}\t{m.chothia_num or ''}\t"
                       f"{m.chothia_icode or ''}\t{m.cdr_type or ''}\n")



# =============================================================================
# Module 1: 清理结构
# =============================================================================

def clean_structure(structure: AtomArray, config: Config, logger: logging.Logger) -> AtomArray:
    """
    清理结构：
    - 去除水分子 (HOH)
    - 去除非蛋白HETATM（可选保留糖链）
    - 处理altloc（保留occupancy最高或第一个）
    """
    logger.info("清理结构...")
    
    # 只保留蛋白
    cleaned = filter_protein_only(structure, logger)
    
    # 处理altloc - 保留第一个或occupancy最高的
    # Biotite的AtomArray没有直接的altloc字段，这里假设已经是单一构象
    
    return cleaned


# =============================================================================
# Module 2: 链角色识别
# =============================================================================

def extract_chain_sequences(structure: AtomArray) -> dict[str, str]:
    """从结构中提取每条链的序列"""
    aa_code = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    chains = {}
    for chain_id in np.unique(structure.chain_id):
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        # 获取CA原子以确定残基顺序
        ca_mask = chain_atoms.atom_name == 'CA'
        ca_atoms = chain_atoms[ca_mask]
        
        seq = ""
        for atom in ca_atoms:
            res_name = atom.res_name
            if res_name in aa_code:
                seq += aa_code[res_name]
            else:
                seq += 'X'
        
        chains[chain_id] = seq
    
    return chains


def extract_chain_sequences_v2(
    structure: AtomArray,
    logger: logging.Logger
) -> dict[str, ChainSequenceInfo]:
    """
    V2: 从结构中提取每条链的序列信息
    返回 ChainSequenceInfo 对象，包含残基列表和过滤信息
    """
    aa_code = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    chain_infos = {}
    
    for chain_id in np.unique(structure.chain_id):
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        # 获取CA原子以确定残基顺序
        ca_mask = chain_atoms.atom_name == 'CA'
        ca_atoms = chain_atoms[ca_mask]
        
        if len(ca_atoms) == 0:
            # 非蛋白链
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
            # 获取插入码（如果存在）
            icode = getattr(atom, 'ins_code', '') or ''
            
            if res_name in aa_code:
                seq += aa_code[res_name]
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
    
    # 统计
    protein_count = sum(1 for c in chain_infos.values() if c.is_protein)
    logger.info(f"提取到 {len(chain_infos)} 条链序列 ({protein_count} 蛋白, {len(chain_infos) - protein_count} 非蛋白)")
    
    return chain_infos


def run_anarci(sequence: str, scheme: str = 'chothia') -> dict | None:
    """
    运行ANARCI进行抗体编号
    返回：编号结果字典 或 None（如果不是抗体）
    """
    try:
        from anarci import anarci
        results = anarci([("seq", sequence)], scheme=scheme, output=False)
        
        if results is None or results[0] is None:
            return None
        
        numbering, alignment_details, hit_tables = results
        
        if numbering[0] is None:
            return None
        
        # 解析结果
        # ANARCI返回: (numbering, alignment_details, hit_tables)
        # numbering = [seq1_results, seq2_results, ...]  # 每个sequence一个
        # seq_results = [(numbered_list, start, end), ...]  # 每个hit一个
        
        # 获取第一个sequence的结果
        first_seq_results = numbering[0]
        if not first_seq_results or len(first_seq_results) == 0:
            return None
        
        # 获取第一个hit
        first_hit = first_seq_results[0]
        numbered_list, start_idx, end_idx = first_hit
        
        # 获取chain_type从 alignment_details（不是 hit_tables！）
        # alignment_details[0] = [{'chain_type': 'H/K/L', 'species': ..., ...}, ...]
        chain_type = 'H'  # default
        if alignment_details and len(alignment_details) > 0:
            first_seq_details = alignment_details[0]
            if first_seq_details and len(first_seq_details) > 0:
                # alignment_details[0][0] 是一个字典
                chain_type = first_seq_details[0].get('chain_type', 'H')
        
        return {
            'chain_type': chain_type,
            'numbering': numbered_list,
            'scheme': scheme
        }
    except ImportError:
        print("ANARCI ImportError - not installed")
        return None
    except Exception as e:
        # 调试：记录错误
        import traceback
        print(f"ANARCI Exception for sequence (len={len(sequence)}): {e}")
        traceback.print_exc()
        return None


def run_anarci_all_hits(sequence: str, scheme: str = 'chothia') -> list[dict] | None:
    """
    运行 ANARCI 并返回所有 hits（用于 scFv 检测）
    scFv 结构中，同一条链会有 VH 和 VL 两个域
    """
    try:
        from anarci import anarci
        results = anarci([("seq", sequence)], scheme=scheme, output=False)
        
        if results is None or results[0] is None:
            return None
        
        numbering, alignment_details, hit_tables = results
        
        if numbering[0] is None:
            return None
        
        first_seq_results = numbering[0]
        if not first_seq_results or len(first_seq_results) == 0:
            return None
        
        first_seq_details = alignment_details[0] if alignment_details else []
        
        # 收集所有 hits
        all_hits = []
        for i, hit in enumerate(first_seq_results):
            numbered_list, start_idx, end_idx = hit
            chain_type = 'H'
            if i < len(first_seq_details):
                chain_type = first_seq_details[i].get('chain_type', 'H')
            
            all_hits.append({
                'chain_type': chain_type,
                'numbering': numbered_list,
                'start': start_idx,
                'end': end_idx,
                'scheme': scheme
            })
        
        return all_hits if all_hits else None
        
    except Exception:
        return None


def detect_scfv(sequence: str, logger: logging.Logger) -> dict | None:
    """
    检测 scFv 结构（同一链上同时有 VH 和 VL 域）
    
    返回:
        {
            'is_scfv': True,
            'vh_domain': {'start': int, 'end': int, 'numbering': ...},
            'vl_domain': {'start': int, 'end': int, 'numbering': ...},
            'vh_first': bool  # VH 在 VL 之前
        }
    或 None 如果不是 scFv
    """
    hits = run_anarci_all_hits(sequence, 'chothia')
    
    if not hits or len(hits) < 2:
        return None
    
    # 找 VH 和 VL
    vh_hits = [h for h in hits if h['chain_type'] == 'H']
    vl_hits = [h for h in hits if h['chain_type'] in ['K', 'L']]
    
    if not vh_hits or not vl_hits:
        return None
    
    # 取第一个 VH 和 VL
    vh = vh_hits[0]
    vl = vl_hits[0]
    
    logger.info(f"  检测到 scFv: VH({vh['start']}-{vh['end']}) + VL({vl['start']}-{vl['end']})")
    
    return {
        'is_scfv': True,
        'vh_domain': vh,
        'vl_domain': vl,
        'vh_first': vh['start'] < vl['start']
    }


def identify_antibody_chains_simple(sequences: dict[str, str], logger: logging.Logger) -> tuple[list[str], list[str]]:
    """
    简单的抗体链识别（基于序列长度和保守残基）
    用于ANARCI不可用时的后备方案
    """
    heavy_chains = []
    light_chains = []
    
    for chain_id, seq in sequences.items():
        length = len(seq)
        
        # 粗略判断：Fv区域约110-130残基
        if 100 <= length <= 250:
            # 检查一些保守残基模式
            # Heavy chain VH通常在位置22有Cys
            # 这只是一个非常粗略的启发式
            if 'CXXW' in seq[20:30] if len(seq) > 30 else False:
                heavy_chains.append(chain_id)
            elif 'C' in seq[20:25] if len(seq) > 25 else False:
                light_chains.append(chain_id)
    
    return heavy_chains, light_chains


def parse_pdb_header_for_antibody_chains(input_file: str, logger: logging.Logger) -> tuple[list[str], list[str]]:
    """
    从PDB header的COMPND记录解析抗体链信息
    查找含有 "HEAVY CHAIN" 和 "LIGHT CHAIN" 关键词的分子
    """
    heavy_chains = []
    light_chains = []
    
    try:
        with open(input_file, 'r') as f:
            current_molecule = ""
            current_chains = []
            
            for line in f:
                if not line.startswith("COMPND"):
                    if line.startswith("SOURCE") or line.startswith("KEYWDS"):
                        break
                    continue
                
                content = line[10:].strip()
                
                if "MOLECULE:" in content:
                    # 处理上一个分子
                    if current_molecule and current_chains:
                        mol_upper = current_molecule.upper()
                        if "HEAVY CHAIN" in mol_upper or "HEAVY" in mol_upper:
                            heavy_chains.extend(current_chains)
                        elif "LIGHT CHAIN" in mol_upper or "LIGHT" in mol_upper:
                            light_chains.extend(current_chains)
                    
                    # 开始新分子
                    current_molecule = content.replace("MOLECULE:", "").strip().rstrip(";")
                    current_chains = []
                elif "CHAIN:" in content:
                    chains_str = content.replace("CHAIN:", "").strip().rstrip(";")
                    current_chains = [c.strip() for c in chains_str.split(",")]
            
            # 处理最后一个分子
            if current_molecule and current_chains:
                mol_upper = current_molecule.upper()
                if "HEAVY CHAIN" in mol_upper or "HEAVY" in mol_upper:
                    heavy_chains.extend(current_chains)
                elif "LIGHT CHAIN" in mol_upper or "LIGHT" in mol_upper:
                    light_chains.extend(current_chains)
        
        if heavy_chains or light_chains:
            logger.info(f"从PDB header解析到: 重链={heavy_chains}, 轻链={light_chains}")
    except Exception as e:
        logger.warning(f"解析PDB header失败: {e}")
    
    return heavy_chains, light_chains


def assign_roles(structure: AtomArray, config: Config, logger: logging.Logger) -> RoleAssignment:
    """
    识别并分配链角色：VH, VL (或VHH), 抗原
    
    优先级:
    1. 手动指定的链ID (config.heavy_chain_id, light_chain_id, antigen_chain_ids)
    2. ANARCI自动识别
    3. PDB header解析
    4. 简单启发式识别
    """
    logger.info("识别链角色...")
    
    # 提取每条链的序列
    sequences = extract_chain_sequences(structure)
    all_chain_ids = list(sequences.keys())
    logger.info(f"提取到 {len(sequences)} 条链的序列: {all_chain_ids}")
    
    heavy_chains = []
    light_chains = []
    antigen_chains = []
    anarci_results = {}
    
    # 方法1: 检查是否有手动指定
    if config.heavy_chain_id or config.light_chain_id or config.antigen_chain_ids:
        logger.info("使用手动指定的链ID...")
        
        if config.heavy_chain_id:
            heavy_chains = [c.strip() for c in config.heavy_chain_id.split(",")]
            logger.info(f"  手动指定重链: {heavy_chains}")
        
        if config.light_chain_id:
            light_chains = [c.strip() for c in config.light_chain_id.split(",")]
            logger.info(f"  手动指定轻链: {light_chains}")
        
        if config.antigen_chain_ids:
            antigen_chains = [c.strip() for c in config.antigen_chain_ids.split(",")]
            logger.info(f"  手动指定抗原链: {antigen_chains}")
        else:
            # 推断抗原链
            antigen_chains = [c for c in all_chain_ids if c not in heavy_chains and c not in light_chains]
            logger.info(f"  推断抗原链: {antigen_chains}")
        
        # ========== 关键修复：手动指定链时也要调用ANARCI获取编号 ==========
        logger.info("调用 ANARCI 获取 Chothia 编号...")
        for chain_id in heavy_chains + light_chains:
            if chain_id in sequences:
                result = run_anarci(sequences[chain_id], config.numbering)
                if result:
                    anarci_results[chain_id] = result
                    numbering = result.get('numbering', [])
                    num_positions = len(numbering) if isinstance(numbering, list) else 0
                    logger.info(f"  链 {chain_id}: 获取到 {result['chain_type']} 链 Chothia 编号 ({num_positions} 个位置)")
                else:
                    logger.warning(f"  链 {chain_id}: ANARCI 编号失败")
    
    # 方法2: ANARCI自动识别
    elif config.use_anarci:
        logger.info("尝试使用ANARCI进行抗体链识别...")
        anarci_available = False
        
        for chain_id, seq in sequences.items():
            result = run_anarci(seq, config.numbering)
            if result:
                anarci_available = True
                anarci_results[chain_id] = result
                chain_type = result['chain_type']
                if chain_type == 'H':
                    heavy_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 重链 (VH)")
                elif chain_type in ['K', 'L']:
                    light_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 轻链 (VL)")
                else:
                    antigen_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 非抗体 -> 抗原")
            else:
                if not anarci_available:
                    # ANARCI可能未安装，跳过尝试其他方法
                    pass
                else:
                    antigen_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 非抗体 -> 抗原")
        
        # 如果ANARCI不可用或未识别到抗体链，尝试其他方法
        if not heavy_chains and not light_chains:
            logger.info("ANARCI未识别到抗体链，尝试从PDB header解析...")
            heavy_chains, light_chains = parse_pdb_header_for_antibody_chains(config.input_file, logger)
            antigen_chains = [c for c in all_chain_ids if c not in heavy_chains and c not in light_chains]
    
    # 方法3: 如果仍然没有识别到，使用PDB header
    if not heavy_chains and not light_chains:
        logger.info("尝试从PDB header解析抗体链...")
        heavy_chains, light_chains = parse_pdb_header_for_antibody_chains(config.input_file, logger)
        antigen_chains = [c for c in all_chain_ids if c not in heavy_chains and c not in light_chains]
    
    # 方法4: 简单启发式（最后手段）
    if not heavy_chains and not light_chains:
        logger.warning("所有自动识别方法均失败，请使用 --heavy_chain 和 --light_chain 手动指定链ID")
        raise ValueError("无法识别抗体链。请使用 --heavy_chain 和 --light_chain 参数手动指定")
    
    # 确定模式
    if config.mode == "auto":
        if heavy_chains and not light_chains:
            mode = "nanobody"
            logger.info("自动检测模式: nanobody (仅有重链)")
        elif heavy_chains and light_chains:
            mode = "antibody"
            logger.info("自动检测模式: antibody")
        else:
            raise ValueError("无法自动确定模式：未识别到抗体链")
    else:
        mode = config.mode
        logger.info(f"使用指定模式: {mode}")
    
    # 选择最佳的VH-VL配对（如果有多个）
    selected_heavy = None
    selected_light = None
    best_contact_score = 0.0
    
    if mode == "antibody" and len(heavy_chains) > 0 and len(light_chains) > 0:
        if len(heavy_chains) == 1 and len(light_chains) == 1:
            selected_heavy = heavy_chains[0]
            selected_light = light_chains[0]
            best_contact_score = compute_chain_contact_score(structure, selected_heavy, selected_light)
        else:
            # 多个候选，计算每对H-L的接触
            logger.info("检测到多个抗体链，计算最佳VH-VL配对...")
            for h in heavy_chains:
                for l in light_chains:
                    score = compute_chain_contact_score(structure, h, l)
                    logger.debug(f"  {h}-{l}: 接触分数 = {score}")
                    if score > best_contact_score:
                        best_contact_score = score
                        selected_heavy = h
                        selected_light = l
        
        logger.info(f"选择VH-VL配对: {selected_heavy}-{selected_light} (接触分数: {best_contact_score})")
        
        # 更新抗原链列表（移除被选中的抗体链）
        antigen_chains = [c for c in all_chain_ids if c != selected_heavy and c != selected_light]
        
    elif mode == "nanobody" and heavy_chains:
        selected_heavy = heavy_chains[0]
        logger.info(f"Nanobody模式，选择重链: {selected_heavy}")
        antigen_chains = [c for c in all_chain_ids if c != selected_heavy]
    
    if mode == "antibody" and (selected_heavy is None or selected_light is None):
        raise ValueError("Antibody模式下无法确定VH-VL配对")
    if mode == "nanobody" and selected_heavy is None:
        raise ValueError("Nanobody模式下无法确定VH链")
    
    logger.info(f"最终分配: 重链={selected_heavy}, 轻链={selected_light}, 抗原链={antigen_chains}")
    
    return RoleAssignment(
        heavy_chain=selected_heavy,
        light_chain=selected_light,
        antigen_chains=antigen_chains,
        mode=mode,
        vh_vl_contact_score=best_contact_score,
        all_heavy_chains=heavy_chains,
        all_light_chains=light_chains,
        anarci_results=anarci_results
    )


def assign_roles_v2(
    structure: AtomArray, 
    config: Config, 
    fasta_records: list[FastaRecord] | None,
    chain_mappings: dict[str, ChainMapping] | None,
    logger: logging.Logger
) -> RoleAssignment:
    """
    V2 链角色识别：融合 ANARCI + FASTA hint
    
    策略:
    1. 如果手动指定链ID -> 直接使用
    2. ANARCI 扫描所有链
    3. FASTA hint 作为辅助提示（不覆盖 ANARCI 结果）
    4. 选择最佳 VH-VL 配对
    5. 标记非选中抗体链为 crystal_packing
    """
    logger.info("识别链角色 (V2: FASTA+ANARCI 融合)...")
    
    # 提取每条链的序列
    sequences = extract_chain_sequences(structure)
    all_chain_ids = list(sequences.keys())
    logger.info(f"提取到 {len(sequences)} 条链的序列: {all_chain_ids}")
    
    heavy_chains = []
    light_chains = []
    antigen_chains = []
    anarci_results = {}
    chain_roles = {}  # V2: 保存每条链的角色信息
    has_scfv = False  # P1.2: 标记是否检测到 scFv
    
    # ============ 方法1: 手动指定优先 ============
    if config.heavy_chain_id or config.light_chain_id or config.antigen_chain_ids:
        logger.info("使用手动指定的链ID...")
        
        if config.heavy_chain_id:
            heavy_chains = [c.strip() for c in config.heavy_chain_id.split(",")]
            logger.info(f"  手动指定重链: {heavy_chains}")
        
        if config.light_chain_id:
            light_chains = [c.strip() for c in config.light_chain_id.split(",")]
            logger.info(f"  手动指定轻链: {light_chains}")
        
        if config.antigen_chain_ids:
            antigen_chains = [c.strip() for c in config.antigen_chain_ids.split(",")]
            logger.info(f"  手动指定抗原链: {antigen_chains}")
        else:
            antigen_chains = [c for c in all_chain_ids if c not in heavy_chains and c not in light_chains]
            logger.info(f"  推断抗原链: {antigen_chains}")
        
        # ANARCI 编号
        logger.info("调用 ANARCI 获取 Chothia 编号...")
        for chain_id in heavy_chains + light_chains:
            if chain_id in sequences:
                result = run_anarci(sequences[chain_id], config.numbering)
                if result:
                    anarci_results[chain_id] = result
                    logger.info(f"  链 {chain_id}: 获取到 {result['chain_type']} 链编号")
                else:
                    logger.warning(f"  链 {chain_id}: ANARCI 编号失败")
    
    # ============ 方法2: FASTA+ANARCI 自动识别 ============
    else:
        logger.info("自动识别链角色 (ANARCI + FASTA hint)...")
        
        # Step 1: ANARCI 扫描所有链
        for chain_id, seq in sequences.items():
            result = run_anarci(seq, config.numbering)
            
            # 获取 FASTA hint
            fasta_hint = None
            if chain_mappings and chain_id in chain_mappings:
                fasta_hint = get_fasta_hint_for_chain(chain_id, chain_mappings, fasta_records or [])
            
            if result:
                anarci_results[chain_id] = result
                chain_type = result['chain_type']
                
                if chain_type == 'H':
                    heavy_chains.append(chain_id)
                    chain_roles[chain_id] = ChainRoleInfo(
                        chain_id=chain_id, is_antibody=True, domain_type='H',
                        anarci_score=1.0, fasta_hint=fasta_hint, final_confidence=0.95
                    )
                    logger.info(f"  链 {chain_id}: 重链 (VH) [ANARCI]" + (f" [FASTA: {fasta_hint}]" if fasta_hint else ""))
                elif chain_type in ['K', 'L']:
                    light_chains.append(chain_id)
                    chain_roles[chain_id] = ChainRoleInfo(
                        chain_id=chain_id, is_antibody=True, domain_type=chain_type,
                        anarci_score=1.0, fasta_hint=fasta_hint, final_confidence=0.95
                    )
                    logger.info(f"  链 {chain_id}: 轻链 (VL) [ANARCI]" + (f" [FASTA: {fasta_hint}]" if fasta_hint else ""))
                else:
                    antigen_chains.append(chain_id)
                    chain_roles[chain_id] = ChainRoleInfo(
                        chain_id=chain_id, is_antibody=False, anarci_score=0.0,
                        fasta_hint=fasta_hint, final_confidence=0.9
                    )
            else:
                # ANARCI 未识别 -> 使用 FASTA hint
                if fasta_hint:
                    if 'heavy_chain' in fasta_hint.lower():
                        heavy_chains.append(chain_id)
                        chain_roles[chain_id] = ChainRoleInfo(
                            chain_id=chain_id, is_antibody=True, domain_type='H',
                            anarci_score=0.0, fasta_hint=fasta_hint, final_confidence=0.6
                        )
                        logger.info(f"  链 {chain_id}: 重链 (VH) [FASTA hint, ANARCI失败]")
                    elif 'light_chain' in fasta_hint.lower():
                        light_chains.append(chain_id)
                        chain_roles[chain_id] = ChainRoleInfo(
                            chain_id=chain_id, is_antibody=True, domain_type='L',
                            anarci_score=0.0, fasta_hint=fasta_hint, final_confidence=0.6
                        )
                        logger.info(f"  链 {chain_id}: 轻链 (VL) [FASTA hint, ANARCI失败]")
                    else:
                        antigen_chains.append(chain_id)
                        chain_roles[chain_id] = ChainRoleInfo(
                            chain_id=chain_id, is_antibody=False, fasta_hint=fasta_hint,
                            final_confidence=0.8
                        )
                        logger.info(f"  链 {chain_id}: 非抗体 -> 抗原 [FASTA: {fasta_hint}]")
                else:
                    antigen_chains.append(chain_id)
                    chain_roles[chain_id] = ChainRoleInfo(
                        chain_id=chain_id, is_antibody=False, final_confidence=0.7
                    )
                    logger.info(f"  链 {chain_id}: 非抗体 -> 抗原")
        
        # ============ P1.2: scFv 检测 ============
        # 如果只有轻链没有重链，可能是 scFv（ANARCI 可能只返回了第一个域）
        if not heavy_chains and light_chains:
            logger.info("未检测到重链，尝试检测 scFv 结构...")
            scfv_found = False
            for chain_id in list(light_chains):  # 复制列表以便修改
                seq = sequences.get(chain_id, "")
                scfv_info = detect_scfv(seq, logger)
                if scfv_info:
                    # 找到 scFv，使用其 VH 和 VL 信息
                    light_chains.remove(chain_id)
                    heavy_chains.append(chain_id)  # scFv 同时作为 H 和 L
                    
                    # 存储 ANARCI 结果（使用 VH 域的编号）
                    anarci_results[chain_id] = {
                        'chain_type': 'H',
                        'numbering': scfv_info['vh_domain']['numbering'],
                        'scheme': 'chothia',
                        'is_scfv': True,
                        'scfv_info': scfv_info
                    }
                    
                    chain_roles[chain_id] = ChainRoleInfo(
                        chain_id=chain_id, is_antibody=True, domain_type='scFv',
                        anarci_score=1.0, final_confidence=0.9
                    )
                    scfv_found = True
                    logger.info(f"  链 {chain_id}: scFv (VH+VL 单链抗体)")
            
            if scfv_found:
                # 更新抗原链
                antigen_chains = [c for c in all_chain_ids if c not in heavy_chains]
                # 标记有 scFv
                has_scfv = True
        
        # 如果仍未识别到抗体链，尝试 PDB header
        if not heavy_chains and not light_chains:
            logger.info("ANARCI+FASTA 未识别到抗体链，尝试 PDB header...")
            heavy_chains, light_chains = parse_pdb_header_for_antibody_chains(config.input_file, logger)
            antigen_chains = [c for c in all_chain_ids if c not in heavy_chains and c not in light_chains]
    
    # ============ 检查结果 ============
    if not heavy_chains and not light_chains:
        logger.warning("所有自动识别方法均失败")
        raise ValueError("无法识别抗体链。请使用 --heavy_chain 和 --light_chain 参数手动指定，或提供 --fasta 文件辅助识别")
    
    # ============ 确定模式 ============
    if config.mode == "auto":
        if has_scfv:
            # scFv 有 VH+VL 在同一条链上
            mode = "scfv"
            logger.info("自动检测模式: scfv (单链抗体 VH-linker-VL)")
        elif heavy_chains and not light_chains:
            mode = "nanobody"
            logger.info("自动检测模式: nanobody (仅有重链)")
        elif heavy_chains and light_chains:
            mode = "antibody"
            logger.info("自动检测模式: antibody")
        else:
            raise ValueError("无法自动确定模式：未识别到抗体链")
    else:
        mode = config.mode
        logger.info(f"使用指定模式: {mode}")
    
    # ============ VH-VL 配对选择 (P1.3: 构建所有配对) ============
    selected_heavy = None
    selected_light = None
    best_contact_score = 0.0
    all_pairings = []  # P1.3: 存储所有有效配对
    
    if mode == "antibody" and len(heavy_chains) > 0 and len(light_chains) > 0:
        # 计算所有 VH-VL 配对的接触分数
        pair_scores = []
        for h in heavy_chains:
            for l in light_chains:
                score = compute_chain_contact_score(structure, h, l)
                pair_scores.append((h, l, score))
                logger.debug(f"  {h}-{l}: 接触分数 = {score}")
        
        # 按分数降序排序
        pair_scores.sort(key=lambda x: -x[2])
        
        # P1.3: 贪婪选择不冲突的配对（每条链只能用一次）
        used_chains = set()
        pairing_index = 1
        
        for h, l, score in pair_scores:
            if h not in used_chains and l not in used_chains:
                if score > 0:  # 有效接触
                    all_pairings.append(AntibodyPairing(
                        heavy_chain=h,
                        light_chain=l,
                        contact_score=score,
                        index=pairing_index,
                        mode="antibody"
                    ))
                    used_chains.add(h)
                    used_chains.add(l)
                    pairing_index += 1
                    logger.info(f"  配对 {pairing_index-1}: {h}-{l} (分数: {score})")
        
        # 第一个配对作为主抗体
        if all_pairings:
            selected_heavy = all_pairings[0].heavy_chain
            selected_light = all_pairings[0].light_chain
            best_contact_score = all_pairings[0].contact_score
        
        if len(all_pairings) > 1:
            logger.info(f"识别到 {len(all_pairings)} 套抗体配对!")
        elif all_pairings:
            logger.info(f"选择VH-VL配对: {selected_heavy}-{selected_light} (接触分数: {best_contact_score})")
        
        # 更新抗原链 - 排除所有抗体链
        antigen_chains = [c for c in all_chain_ids 
                        if c not in heavy_chains and c not in light_chains]
        
    elif mode == "nanobody" and heavy_chains:
        # Nanobody: 每个 VH 都是一个独立配对
        for i, h in enumerate(heavy_chains, 1):
            score = 0.0  # nanobody 不需要接触分数
            all_pairings.append(AntibodyPairing(
                heavy_chain=h,
                light_chain=None,
                contact_score=score,
                index=i,
                mode="nanobody"
            ))
        
        selected_heavy = heavy_chains[0]
        logger.info(f"Nanobody模式，识别到 {len(all_pairings)} 个 VHH")
        antigen_chains = [c for c in all_chain_ids if c not in heavy_chains]
    
    elif mode == "scfv" and heavy_chains:
        # scFv: 每个 scFv 链同时包含 VH 和 VL
        for i, h in enumerate(heavy_chains, 1):
            score = 0.0  # scFv 不需要配对分数
            all_pairings.append(AntibodyPairing(
                heavy_chain=h,
                light_chain=None,  # scFv 模式下 light_chain 在 heavy_chain 同链内
                contact_score=score,
                index=i,
                mode="scfv"
            ))
        
        selected_heavy = heavy_chains[0]
        logger.info(f"scFv模式，识别到 {len(all_pairings)} 个 scFv")
        antigen_chains = [c for c in all_chain_ids if c not in heavy_chains]
    
    if mode == "antibody" and (selected_heavy is None or selected_light is None):
        raise ValueError("Antibody模式下无法确定VH-VL配对")
    if mode == "nanobody" and selected_heavy is None:
        raise ValueError("Nanobody模式下无法确定VH链")
    if mode == "scfv" and selected_heavy is None:
        raise ValueError("scFv模式下无法确定scFv链")
    
    logger.info(f"主抗体: 重链={selected_heavy}, 轻链={selected_light}")
    logger.info(f"抗原链: {antigen_chains}")
    
    return RoleAssignment(
        heavy_chain=selected_heavy,
        light_chain=selected_light,
        antigen_chains=antigen_chains,
        mode=mode,
        vh_vl_contact_score=best_contact_score,
        all_heavy_chains=heavy_chains,
        all_light_chains=light_chains,
        anarci_results=anarci_results,
        all_pairings=all_pairings  # P1.3
    )


def compute_chain_contact_score(structure: AtomArray, chain1: str, chain2: str, cutoff: float = 8.0) -> float:
    """计算两条链之间的接触分数"""
    mask1 = structure.chain_id == chain1
    mask2 = structure.chain_id == chain2
    
    coords1 = structure[mask1].coord
    coords2 = structure[mask2].coord
    
    if len(coords1) == 0 or len(coords2) == 0:
        return 0.0
    
    # 使用KDTree快速计算接触
    tree = KDTree(coords2)
    contacts = tree.query_ball_point(coords1, r=cutoff)
    
    total_contacts = sum(len(c) for c in contacts)
    return float(total_contacts)


# =============================================================================
# Module 3: 导出 framework_HLT.pdb
# =============================================================================

def get_chain_residue_mapping(structure: AtomArray, chain_id: str, anarci_result: dict | None = None) -> list[tuple[int, str, str]]:
    """
    获取链的残基映射：[(原始res_num, icode, Chothia编号), ...]
    如果没有ANARCI结果，使用原始编号
    """
    chain_mask = structure.chain_id == chain_id
    chain_atoms = structure[chain_mask]
    
    # 获取唯一残基
    seen = set()
    residues = []
    for atom in chain_atoms:
        key = (atom.res_id, atom.ins_code if hasattr(atom, 'ins_code') else '')
        if key not in seen:
            seen.add(key)
            icode = atom.ins_code if hasattr(atom, 'ins_code') else ''
            residues.append((atom.res_id, icode, atom.res_name))
    
    return residues


def identify_cdr_residues(anarci_result: dict, cdr_ranges: dict) -> dict[str, list[int]]:
    """
    根据ANARCI编号结果识别CDR残基的位置
    返回：{CDR名: [Chothia编号列表]}
    """
    if not anarci_result:
        return {}
    
    chain_type = anarci_result['chain_type']
    numbering = anarci_result['numbering']
    
    cdr_residues = {}
    ranges = cdr_ranges.get('H' if chain_type == 'H' else 'L', {})
    
    for cdr_name, (start, end) in ranges.items():
        cdr_residues[cdr_name] = []
        for pos, aa in numbering:
            chothia_num = pos[0]  # (num, icode)
            if aa != '-' and start <= chothia_num <= end:
                cdr_residues[cdr_name].append(pos)
    
    return cdr_residues


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
    - 末尾写REMARK PDBinfo-LABEL
    """
    logger.info("导出 framework_HLT.pdb...")
    
    output_path = os.path.join(config.output_dir, "framework_HLT.pdb")
    
    atom_list = []
    cdr_boundaries = []
    residue_counter = 1
    cdr_abs_indices = {}  # {CDR名: [abs_index列表]}
    
    chain_order = ['H']
    chain_mapping = {assignment.heavy_chain: 'H'}
    is_scfv = False
    scfv_info = None
    
    if assignment.mode == "scfv":
        # scFv: 同一条链包含 VH 和 VL，需要拆分
        chain_order = ['H', 'L']
        chain_mapping = {assignment.heavy_chain: 'H', assignment.heavy_chain: 'L'}
        is_scfv = True
        # 获取 scFv 域信息
        anarci_result = assignment.anarci_results.get(assignment.heavy_chain)
        if anarci_result and anarci_result.get('is_scfv'):
            scfv_info = anarci_result.get('scfv_info')
        logger.info(f"  scFv 模式: 将拆分链 {assignment.heavy_chain} 为 H (VH) + L (VL)")
    elif assignment.mode == "antibody" and assignment.light_chain:
        chain_order.append('L')
        chain_mapping[assignment.light_chain] = 'L'
    
    # 获取CDR范围
    cdr_ranges = CHOTHIA_CDR_RANGES
    
    for new_chain_id in chain_order:
        # 找到原始链ID
        orig_chain = assignment.heavy_chain  # scFv 模式下都是同一条链
        
        if not is_scfv:
            # 非 scFv，按原来的逻辑
            orig_chain = None
            for orig, new in chain_mapping.items():
                if new == new_chain_id:
                    orig_chain = orig
                    break
            if orig_chain is None:
                continue
        
        chain_mask = structure.chain_id == orig_chain
        chain_atoms = structure[chain_mask]
        
        # 获取ANARCI结果
        anarci_result = assignment.anarci_results.get(orig_chain)
        
        # ============ scFv 特殊处理 ============
        if is_scfv and scfv_info:
            if new_chain_id == 'H':
                # 使用 VH 域的 numbering
                domain_info = scfv_info['vh_domain']
                numbering = domain_info.get('numbering', [])
                domain_start = domain_info['start']
                domain_end = domain_info['end']
                logger.info(f"  导出 VH 域: 残基 {domain_start}-{domain_end}")
            else:
                # 使用 VL 域的 numbering
                domain_info = scfv_info['vl_domain']
                numbering = domain_info.get('numbering', [])
                domain_start = domain_info['start']
                domain_end = domain_info['end']
                logger.info(f"  导出 VL 域: 残基 {domain_start}-{domain_end}")
            
            # 重新设置 anarci_result 以使用正确域的 numbering
            anarci_result = {'numbering': numbering, 'chain_type': 'H' if new_chain_id == 'H' else 'L'}
        
        # 确定当前链的CDR范围
        if new_chain_id == 'H':
            curr_cdr_ranges = cdr_ranges.get('H', {})
            chain_type = 'H'
        else:
            curr_cdr_ranges = cdr_ranges.get('L', {})
            chain_type = 'L'
        
        # 初始化CDR索引跟踪
        for cdr_name in curr_cdr_ranges.keys():
            cdr_abs_indices[cdr_name] = []
        
        # ============ scFv: 过滤残基到目标域 ============
        if is_scfv and scfv_info:
            # 获取域边界（0-indexed 残基位置）
            domain_start = domain_info['start']
            domain_end = domain_info['end']
            
            # 过滤 chain_atoms 到目标域
            # 首先获取所有唯一的残基 ID
            residue_ids = []
            seen_res = set()
            for atom in chain_atoms:
                res_key = (atom.res_id, atom.res_name)
                if res_key not in seen_res:
                    seen_res.add(res_key)
                    residue_ids.append(atom.res_id)
            
            # 按序列位置过滤（domain_start/end 是序列中的位置，不是 res_id）
            # 所以我们需要选择第 domain_start 到 domain_end 个残基
            target_res_ids = set()
            for i, res_id in enumerate(residue_ids):
                if domain_start <= i <= domain_end:
                    target_res_ids.add(res_id)
            
            # 过滤原子
            domain_mask = np.array([atom.res_id in target_res_ids for atom in chain_atoms])
            chain_atoms = chain_atoms[domain_mask]
            logger.info(f"  scFv 域过滤: {len(target_res_ids)} 残基, {sum(domain_mask)} 原子")
        
        # ============ 关键修复：使用 ANARCI Chothia 编号 ============
        if anarci_result:
            # 有ANARCI结果 - 使用Chothia编号
            numbering = anarci_result['numbering']  # [(pos, aa), ...]
            
            # 创建迭代器，跳过gap
            numbering_idx = 0
            
            # 处理每个残基
            for residue in residue_iter(chain_atoms):
                # 找到对应的Chothia位置（跳过gap）
                chothia_num, chothia_icode = None, None
                
                while numbering_idx < len(numbering):
                    pos, aa = numbering[numbering_idx]
                    numbering_idx += 1
                    
                    if aa != '-':  # 跳过对齐gap
                        chothia_num, chothia_icode = pos
                        break
                
                if chothia_num is None:
                    logger.warning(f"链 {orig_chain} 残基超出ANARCI编号范围，跳过")
                    continue
                
                # 使用Chothia编号进行裁剪
                crop_limit = config.h_crop if new_chain_id == 'H' else config.l_crop
                if chothia_num > crop_limit:
                    continue  # 超出Fv区域
                
                # 使用Chothia编号判断CDR
                for cdr_name, (start, end) in curr_cdr_ranges.items():
                    if start <= chothia_num <= end:
                        cdr_abs_indices[cdr_name].append(residue_counter)
                        break
                
                # 重新编号并写入
                residue.res_id = np.full(len(residue), residue_counter)
                residue.chain_id = np.full(len(residue), new_chain_id)
                # 清除插入码
                if hasattr(residue, 'ins_code'):
                    residue.ins_code = np.full(len(residue), '')
                
                atom_list.append(residue)
                residue_counter += 1
        
        else:
            # 没有ANARCI结果 - fail-fast
            logger.error(f"链 {orig_chain} 缺少ANARCI编号，无法可靠识别CDR")
            logger.error("建议：")
            logger.error("  1. 安装ANARCI: conda install -c bioconda anarci")
            logger.error("  2. 或使用 --heavy_chain 和 --light_chain 手动指定链ID")
            raise ValueError(
                f"链 {orig_chain} 缺少ANARCI编号结果。"
                f"无法可靠识别CDR和Fv边界。请安装ANARCI或手动指定链ID。"
            )
    
    # 构建CDR边界信息
    for cdr_name, indices in cdr_abs_indices.items():
        if indices:
            chain_type = 'H' if cdr_name.startswith('H') else 'L'
            cdr_range = curr_cdr_ranges.get(cdr_name, (0, 0)) if chain_type == new_chain_id else CHOTHIA_CDR_RANGES.get(chain_type, {}).get(cdr_name, (0, 0))
            
            cdr_boundaries.append(CDRBoundary(
                cdr_name=cdr_name,
                chothia_start=cdr_range[0] if cdr_range else 0,
                chothia_end=cdr_range[1] if cdr_range else 0,
                abs_start=min(indices),
                abs_end=max(indices),
                length=len(indices)
            ))
    
    # 合并原子数组
    if atom_list:
        # 将所有残基合并
        all_atoms = []
        for residue in atom_list:
            for atom in residue:
                all_atoms.append(atom)
        hlt_structure = array(all_atoms)
        
        # 写入PDB
        pdb_file = PDBFile()
        pdb_file.set_structure(hlt_structure)
        
        with open(output_path, 'w') as f:
            pdb_file.write(f)
            
            # 写入CDR REMARK
            for cdr in sorted(cdr_boundaries, key=lambda x: x.cdr_name):
                for abs_idx in range(cdr.abs_start, cdr.abs_end + 1):
                    f.write(f"REMARK PDBinfo-LABEL: {abs_idx:4d} {cdr.cdr_name}\n")
        
        logger.info(f"  写入 {len(hlt_structure)} 个原子到 {output_path}")
        logger.info(f"  CDR边界: {[(c.cdr_name, c.abs_start, c.abs_end) for c in cdr_boundaries]}")
    
    return output_path, cdr_boundaries


# =============================================================================
# Module 4: 导出 target_HLT.pdb
# =============================================================================

def detect_target_resnum_conflicts(
    structure: AtomArray,
    antigen_chains: list[str],
    logger: logging.Logger
) -> dict:
    """
    V2 P0.9: 检测多链 target 合并为 T 后的 resnum 冲突
    
    当多条抗原链合并为单一 T 链时，如果不同链有相同的 (resnum, icode)，
    会导致 hotspot 残基无法唯一定位。
    
    返回:
        {
            'has_conflict': bool,
            'conflict_count': int,
            'conflicts': [{'resnum': int, 'icode': str, 'chains': [str, ...]}, ...],
            'chain_ranges': {chain_id: (min_resnum, max_resnum), ...}
        }
    """
    logger.info("检测 target resnum 冲突...")
    
    # 收集每条链的 (resnum, icode) 集合
    chain_resnums = {}  # {chain_id: set of (resnum, icode)}
    chain_ranges = {}   # {chain_id: (min, max)}
    
    for chain_id in antigen_chains:
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        # 获取 CA 原子的残基编号
        ca_mask = chain_atoms.atom_name == 'CA'
        ca_atoms = chain_atoms[ca_mask]
        
        resnums = set()
        min_res = float('inf')
        max_res = float('-inf')
        
        for atom in ca_atoms:
            res_id = atom.res_id
            icode = getattr(atom, 'ins_code', '') or ''
            resnums.add((res_id, icode))
            min_res = min(min_res, res_id)
            max_res = max(max_res, res_id)
        
        chain_resnums[chain_id] = resnums
        if resnums:
            chain_ranges[chain_id] = (int(min_res), int(max_res))
    
    # 检测冲突
    all_resnums = {}  # {(resnum, icode): [chain_ids]}
    for chain_id, resnums in chain_resnums.items():
        for key in resnums:
            if key not in all_resnums:
                all_resnums[key] = []
            all_resnums[key].append(chain_id)
    
    # 找出冲突
    conflicts = []
    for (resnum, icode), chains in all_resnums.items():
        if len(chains) > 1:
            conflicts.append({
                'resnum': resnum,
                'icode': icode,
                'chains': chains
            })
    
    has_conflict = len(conflicts) > 0
    
    if has_conflict:
        logger.warning(f"  ⚠ 检测到 {len(conflicts)} 个 resnum 冲突!")
        for c in conflicts[:5]:  # 只显示前5个
            logger.warning(f"    resnum {c['resnum']}{c['icode']}: 链 {c['chains']}")
        if len(conflicts) > 5:
            logger.warning(f"    ... 还有 {len(conflicts) - 5} 个冲突")
    else:
        logger.info("  ✓ 无 resnum 冲突")
    
    return {
        'has_conflict': has_conflict,
        'conflict_count': len(conflicts),
        'conflicts': conflicts,
        'chain_ranges': chain_ranges
    }


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
    - P1.1: 如果有冲突，自动重编号为连续 1..N
    
    返回: (output_path, resnum_mapping)
        resnum_mapping = {(orig_chain, orig_resnum, icode): new_resnum, ...}
    """
    logger.info("导出 target_HLT.pdb...")
    
    output_path = os.path.join(config.output_dir, "target_HLT.pdb")
    
    # 判断是否需要重编号
    need_renumber = (
        resnum_conflict_info is not None and 
        resnum_conflict_info.get('has_conflict', False)
    )
    
    if need_renumber:
        logger.info("  检测到冲突，执行连续重编号 (P1.1)...")
    
    atom_list = []
    resnum_mapping = {}  # {(orig_chain, orig_resnum, icode): new_resnum}
    new_resnum = 1
    current_residue_key = None
    
    for orig_chain in assignment.antigen_chains:
        chain_mask = structure.chain_id == orig_chain
        chain_atoms = structure[chain_mask]
        
        # 按残基顺序处理
        chain_atoms = chain_atoms.copy()
        chain_atoms.chain_id = np.full(len(chain_atoms), 'T')
        
        for atom in chain_atoms:
            orig_resnum = atom.res_id
            icode = getattr(atom, 'ins_code', '') or ''
            key = (str(orig_chain), int(orig_resnum), icode)
            
            if need_renumber:
                # 检查是否是新残基
                if key != current_residue_key:
                    current_residue_key = key
                    resnum_mapping[key] = new_resnum
                    new_resnum += 1
                
                # 更新残基编号
                atom.res_id = resnum_mapping[key]
            else:
                # 不重编号，直接记录映射
                if key not in resnum_mapping:
                    resnum_mapping[key] = int(orig_resnum)
            
            atom_list.append(atom)
    
    if atom_list:
        target_structure = array(atom_list)
        
        pdb_file = PDBFile()
        pdb_file.set_structure(target_structure)
        pdb_file.write(output_path)
        
        logger.info(f"  写入 {len(target_structure)} 个原子到 {output_path}")
        
        if need_renumber:
            logger.info(f"  重编号范围: T:1 - T:{new_resnum - 1}")
            
            # 保存映射文件
            mapping_path = os.path.join(config.output_dir, "target_resmap.tsv")
            with open(mapping_path, 'w') as f:
                f.write("orig_chain\torig_resnum\ticode\tnew_resnum\n")
                for (oc, orn, ic), nrn in sorted(resnum_mapping.items(), key=lambda x: x[1]):
                    f.write(f"{oc}\t{orn}\t{ic}\t{nrn}\n")
            logger.info(f"  映射表: {mapping_path}")
    
    return output_path, resnum_mapping


def update_residue_numbers(
    residue_tags: list[ResidueTag],
    resnum_mapping: dict
) -> list[ResidueTag]:
    """
    P1.1: 更新残基标签的编号（基于重编号映射）
    
    残基标签的 chain_id 已经是 'T'，但 res_num 是原始编号。
    需要根据 resnum_mapping 更新为新编号。
    """
    updated = []
    
    for tag in residue_tags:
        # 在 mapping 中查找：遍历所有原始 key 找到匹配的 res_num
        new_num = tag.res_num
        found = False
        
        for (orig_chain, orig_resnum, icode), new_resnum in resnum_mapping.items():
            if orig_resnum == tag.res_num and icode == (tag.icode or ''):
                new_num = new_resnum
                found = True
                break
        
        updated.append(ResidueTag(
            chain_id=tag.chain_id,
            res_num=new_num,
            icode=tag.icode,
            res_name=tag.res_name
        ))
    
    return updated


# =============================================================================
# Module 5: 界面计算
# =============================================================================

def compute_epitope_patch(
    structure: AtomArray,
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger
) -> list[ResidueTag]:
    """
    计算抗原界面残基（epitope）
    """
    logger.info(f"计算界面残基 (cutoff={config.interface_cutoff}Å)...")
    
    # 获取抗体原子坐标
    ab_chains = [assignment.heavy_chain]
    if assignment.light_chain:
        ab_chains.append(assignment.light_chain)
    
    ab_mask = np.isin(structure.chain_id, ab_chains)
    ab_coords = structure[ab_mask].coord
    
    if len(ab_coords) == 0:
        logger.warning("未找到抗体原子")
        return []
    
    # 构建KDTree
    ab_tree = KDTree(ab_coords)
    
    # 查找抗原界面残基
    epitope_residues = set()
    
    for ag_chain in assignment.antigen_chains:
        chain_mask = structure.chain_id == ag_chain
        chain_atoms = structure[chain_mask]
        
        for residue in residue_iter(chain_atoms):
            # 检查是否有任何原子在cutoff内
            for atom in residue:
                distances, _ = ab_tree.query(atom.coord.reshape(1, -1), k=1)
                if distances[0] <= config.interface_cutoff:
                    icode = atom.ins_code if hasattr(atom, 'ins_code') else ''
                    tag = ResidueTag(
                        chain_id='T',  # 统一为T
                        res_num=atom.res_id,
                        icode=icode,
                        res_name=atom.res_name
                    )
                    epitope_residues.add(tag)
                    break
    
    epitope_list = sorted(epitope_residues, key=lambda x: (x.chain_id, x.res_num))
    logger.info(f"  找到 {len(epitope_list)} 个界面残基")
    
    return epitope_list


def detect_hetatm_interface(
    input_file: str,
    ab_chains: list[str],
    cutoff: float,
    logger: logging.Logger
) -> dict:
    """
    P2.1: 检测 HETATM 配体界面（非蛋白抗原）
    
    用于诊断当蛋白界面为空时，是否有小分子/糖/核酸配体与抗体接触
    
    返回:
        {
            'has_hetatm_interface': bool,
            'hetatm_residues': [(chain, resnum, resname), ...],
            'hetatm_types': {'ligand': n, 'sugar': n, 'nucleic': n}
        }
    """
    try:
        from biotite.structure.io.pdbx import CIFFile, get_structure
        
        # 读取完整结构（不过滤 HETATM）
        cif = CIFFile.read(input_file)
        structure = get_structure(cif, model=1)
        
        # 获取抗体原子
        ab_mask = np.isin(structure.chain_id, ab_chains)
        ab_coords = structure[ab_mask].coord
        
        if len(ab_coords) == 0:
            return {'has_hetatm_interface': False, 'hetatm_residues': [], 'hetatm_types': {}}
        
        ab_tree = KDTree(ab_coords)
        
        # 查找接近抗体的 HETATM
        hetatm_residues = []
        hetatm_types = {'ligand': 0, 'sugar': 0, 'nucleic': 0, 'metal': 0}
        
        # 常见糖类残基
        sugar_names = {'NAG', 'MAN', 'BMA', 'FUC', 'GAL', 'GLC', 'SIA'}
        # 核酸
        nucleic_names = {'A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT'}
        # 金属
        metal_names = {'ZN', 'MG', 'CA', 'FE', 'MN', 'CU', 'NI', 'CO'}
        
        seen = set()
        for atom in structure:
            # 只检查 HETATM（hetero=True），排除标准蛋白残基
            if not atom.hetero:
                continue
            
            # 检查距离
            dist, _ = ab_tree.query(atom.coord.reshape(1, -1), k=1)
            if dist[0] <= cutoff:
                key = (atom.chain_id, atom.res_id, atom.res_name)
                if key not in seen:
                    seen.add(key)
                    hetatm_residues.append(key)
                    
                    res_name = atom.res_name.upper()
                    if res_name in sugar_names:
                        hetatm_types['sugar'] += 1
                    elif res_name in nucleic_names:
                        hetatm_types['nucleic'] += 1
                    elif res_name in metal_names:
                        hetatm_types['metal'] += 1
                    elif len(res_name) <= 3 and res_name not in {'HOH', 'WAT'}:
                        hetatm_types['ligand'] += 1
        
        # 转换为原生 Python 类型（避免 JSON 序列化问题）
        hetatm_residues_native = [
            (str(chain), int(resid), str(resname)) 
            for chain, resid, resname in hetatm_residues[:20]
        ]
        
        return {
            'has_hetatm_interface': len(hetatm_residues) > 0,
            'hetatm_residues': hetatm_residues_native,
            'hetatm_types': hetatm_types
        }
    except Exception as e:
        logger.debug(f"HETATM 检测失败: {e}")
        return {'has_hetatm_interface': False, 'hetatm_residues': [], 'hetatm_types': {}}


# =============================================================================
# Module 6: Hotspots 选择
# =============================================================================

def select_hotspots_contact_score(
    structure: AtomArray,
    epitope_patch: list[ResidueTag],
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger
) -> list[HotspotResidue]:
    """
    使用contact_score方法选择热点
    score = w1 * contact_count + w2 * (1/min_distance)
    """
    logger.info("使用 contact_score 方法选择热点...")
    
    # 获取抗体原子
    ab_chains = [assignment.heavy_chain]
    if assignment.light_chain:
        ab_chains.append(assignment.light_chain)
    
    ab_mask = np.isin(structure.chain_id, ab_chains)
    ab_atoms = structure[ab_mask]
    ab_coords = ab_atoms.coord
    ab_tree = KDTree(ab_coords)
    
    hotspots = []
    
    for tag in epitope_patch:
        # 找到抗原中对应的残基
        for ag_chain in assignment.antigen_chains:
            chain_mask = structure.chain_id == ag_chain
            res_mask = chain_mask & (structure.res_id == tag.res_num)
            
            res_atoms = structure[res_mask]
            if len(res_atoms) == 0:
                continue
            
            # 计算接触数和最小距离
            contact_count = 0
            min_dist = float('inf')
            
            for atom in res_atoms:
                distances, _ = ab_tree.query(atom.coord.reshape(1, -1), k=1)
                d = distances[0]
                if d <= config.interface_cutoff:
                    contact_count += 1
                if d < min_dist:
                    min_dist = d
            
            # 计算分数
            if min_dist > 0:
                score = contact_count + 1.0 / min_dist
                hotspots.append(HotspotResidue(
                    tag=tag,
                    score=score,
                    contact_count=contact_count,
                    min_distance=min_dist
                ))
            break
    
    # 排序并选择top-k
    hotspots.sort(key=lambda x: x.score, reverse=True)
    
    # 空间去重
    if config.space_dedup_angstrom > 0:
        hotspots = spatial_dedup_hotspots(structure, hotspots, config, assignment, logger)
    
    # 选择top-k
    selected = hotspots[:config.hotspots_k]
    logger.info(f"  选择了 {len(selected)} 个热点")
    
    return selected


def spatial_dedup_hotspots(
    structure: AtomArray,
    hotspots: list[HotspotResidue],
    config: Config,
    assignment: RoleAssignment,
    logger: logging.Logger
) -> list[HotspotResidue]:
    """空间去重：确保热点之间有足够距离"""
    if not hotspots:
        return hotspots
    
    # 获取每个热点的CA坐标
    def get_ca_coord(tag: ResidueTag) -> np.ndarray | None:
        for ag_chain in assignment.antigen_chains:
            chain_mask = structure.chain_id == ag_chain
            res_mask = chain_mask & (structure.res_id == tag.res_num)
            ca_mask = res_mask & (structure.atom_name == 'CA')
            ca_atoms = structure[ca_mask]
            if len(ca_atoms) > 0:
                return ca_atoms.coord[0]
        return None
    
    selected = []
    for hs in hotspots:
        coord = get_ca_coord(hs.tag)
        if coord is None:
            continue
        
        # 检查与已选热点的距离
        too_close = False
        for sel_hs in selected:
            sel_coord = get_ca_coord(sel_hs.tag)
            if sel_coord is not None:
                dist = np.linalg.norm(coord - sel_coord)
                if dist < config.space_dedup_angstrom:
                    too_close = True
                    break
        
        if not too_close:
            selected.append(hs)
    
    return selected


def select_hotspots_rfantibody_cbeta5(
    structure: AtomArray,
    epitope_patch: list[ResidueTag],
    cdr_boundaries: list[CDRBoundary],
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger
) -> list[HotspotResidue]:
    """
    RFantibody训练时使用的热点定义：
    计算每个target残基到最近5个CDR残基的平均Cβ距离
    选择 avg_cb_dist < 8.0 Å 的残基
    """
    logger.info("使用 rfantibody_cbeta5 方法选择热点...")
    
    # 1. 收集所有抗体链残基的Cβ/Cα坐标
    cdr_cb_coords = []
    
    for chain_id in [assignment.heavy_chain, assignment.light_chain]:
        if chain_id is None:
            continue
        
        chain_mask = structure.chain_id == chain_id
        chain_atoms = structure[chain_mask]
        
        # 获取每个残基的Cβ（GLY用Cα）
        for residue in residue_iter(chain_atoms):
            # Cβ优先
            cb_mask = residue.atom_name == 'CB'
            if any(cb_mask):
                cdr_cb_coords.append(residue[cb_mask].coord[0])
            else:
                # GLY 用 Cα
                ca_mask = residue.atom_name == 'CA'
                if any(ca_mask):
                    cdr_cb_coords.append(residue[ca_mask].coord[0])
    
    if len(cdr_cb_coords) < 5:
        logger.warning(f"抗体残基不足5个 ({len(cdr_cb_coords)})，无法使用 cbeta5 方法")
        logger.warning("回退到 contact_score 方法")
        return select_hotspots_contact_score(structure, epitope_patch, assignment, config, logger)
    
    cdr_cb_array = np.array(cdr_cb_coords)
    cdr_tree = KDTree(cdr_cb_array)
    
    logger.info(f"  收集了 {len(cdr_cb_coords)} 个抗体残基的Cβ/Cα坐标")
    
    # 2. 对每个epitope残基计算到最近5个的平均距离
    hotspots = []
    for tag in epitope_patch:
        for ag_chain in assignment.antigen_chains:
            res_mask = (structure.chain_id == ag_chain) & (structure.res_id == tag.res_num)
            res_atoms = structure[res_mask]
            
            if len(res_atoms) == 0:
                continue
            
            # 获取Cβ（GLY用Cα）
            cb_mask = res_atoms.atom_name == 'CB'
            if any(cb_mask):
                cb_coord = res_atoms[cb_mask].coord[0]
            else:
                ca_mask = res_atoms.atom_name == 'CA'
                if any(ca_mask):
                    cb_coord = res_atoms[ca_mask].coord[0]
                else:
                    continue
            
            # 查询最近5个抗体残基
            k = min(5, len(cdr_cb_coords))
            distances, _ = cdr_tree.query(cb_coord.reshape(1, -1), k=k)
            avg_cb_dist = np.mean(distances[0])
            min_dist = np.min(distances[0])
            
            # score = -avg_dist（越近分数越高）
            hotspots.append(HotspotResidue(
                tag=tag,
                score=-avg_cb_dist,
                avg_cb_dist_to_cdr5=avg_cb_dist,
                min_distance=min_dist,
                contact_count=0
            ))
            break
    
    # 3. 过滤 avg_cb_dist  8.0 Å
    hotspots = [h for h in hotspots if h.avg_cb_dist_to_cdr5 < 8.0]
    logger.info(f"  过滤 avg_cb_dist < 8Å: {len(hotspots)} 个残基")
    
    # 4. 排序（按score降序）
    hotspots.sort(key=lambda x: x.score, reverse=True)
    
    # 5. 空间去重
    if config.space_dedup_angstrom > 0:
        hotspots = spatial_dedup_hotspots(structure, hotspots, config, assignment, logger)
        logger.info(f"  空间去重后: {len(hotspots)} 个热点")
    
    # 6. 取top-k
    hotspots = hotspots[:config.hotspots_k]
    logger.info(f"  选择了 {len(hotspots)} 个热点 (cbeta5 方法)")
    
    return hotspots


def select_hotspots(
    structure: AtomArray,
    epitope_patch: list[ResidueTag],
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger,
    cdr_boundaries: list[CDRBoundary] = None
) -> list[HotspotResidue]:
    """选择热点（调度不同方法）"""
    if config.hotspots_method == "contact_score":
        return select_hotspots_contact_score(structure, epitope_patch, assignment, config, logger)
    elif config.hotspots_method == "rfantibody_cbeta5":
        if cdr_boundaries is None:
            logger.warning("cbeta5方法需要CDR boundaries，回退到contact_score")
            return select_hotspots_contact_score(structure, epitope_patch, assignment, config, logger)
        return select_hotspots_rfantibody_cbeta5(structure, epitope_patch, cdr_boundaries, assignment, config, logger)
    else:
        logger.warning(f"未知方法 {config.hotspots_method}，使用 contact_score")
        return select_hotspots_contact_score(structure, epitope_patch, assignment, config, logger)


# =============================================================================
# Module 7: Design Loops 生成
# =============================================================================

def build_design_loops(
    cdr_boundaries: list[CDRBoundary],
    config: Config,
    logger: logging.Logger
) -> str:
    """
    生成 antibody.design_loops 字符串
    默认：固定当前长度
    """
    logger.info("生成 design_loops...")
    
    loop_specs = []
    for cdr in sorted(cdr_boundaries, key=lambda x: x.cdr_name):
        # 默认固定长度
        loop_specs.append(f"{cdr.cdr_name}:{cdr.length}")
    
    design_loops_str = f"antibody.design_loops=[{','.join(loop_specs)}]"
    logger.info(f"  {design_loops_str}")
    
    return design_loops_str


# =============================================================================
# Module 8: 输出文件生成
# =============================================================================

def write_epitope_patch(epitope_patch: list[ResidueTag], output_dir: str, logger: logging.Logger):
    """写入 epitope_patch.txt"""
    epitopes_dir = os.path.join(output_dir, "epitopes")
    os.makedirs(epitopes_dir, exist_ok=True)
    
    output_path = os.path.join(epitopes_dir, "epitope_patch.txt")
    with open(output_path, 'w') as f:
        for tag in epitope_patch:
            f.write(f"{tag.to_pdb_format()}\n")
    
    logger.info(f"  写入 epitope_patch.txt ({len(epitope_patch)} 个残基)")


def write_hotspots(
    hotspots: list[HotspotResidue],
    output_dir: str,
    logger: logging.Logger
):
    """写入 hotspots 相关文件（兼容 insertion code）"""
    epitopes_dir = os.path.join(output_dir, "epitopes")
    os.makedirs(epitopes_dir, exist_ok=True)
    
    # 详细hotspots表
    hotspots_path = os.path.join(epitopes_dir, "hotspots.txt")
    with open(hotspots_path, 'w') as f:
        f.write("residue\tres_name\tscore\tcontact_count\tmin_distance\tavg_cb_dist_to_cdr5\n")
        for hs in hotspots:
            avg_cb = f"{hs.avg_cb_dist_to_cdr5:.4f}" if hs.avg_cb_dist_to_cdr5 else "N/A"
            f.write(f"{hs.tag.to_pdb_format()}\t{hs.tag.res_name}\t{hs.score:.4f}\t"
                   f"{hs.contact_count}\t{hs.min_distance:.4f}\t{avg_cb}\n")
    
    # 检查是否有insertion code
    has_icode = any(hs.tag.icode for hs in hotspots)
    
    # RFantibody格式 - digits_only版本（默认）
    rfab_digits_path = os.path.join(epitopes_dir, "hotspots_rfantibody.txt")
    hotspot_digits = [hs.tag.to_rfantibody_format() for hs in hotspots]
    rfab_digits_str = f"ppi.hotspot_res=[{','.join(hotspot_digits)}]"
    with open(rfab_digits_path, 'w') as f:
        f.write(rfab_digits_str + "\n")
    
    # 如果有insertion code，额外生成with_icode版本
    recommended = "digits_only"
    if has_icode:
        rfab_icode_path = os.path.join(epitopes_dir, "hotspots_rfantibody_with_icode.txt")
        hotspot_icode = [f"{hs.tag.chain_id}{hs.tag.res_num}{hs.tag.icode}" for hs in hotspots]
        rfab_icode_str = f"ppi.hotspot_res=[{','.join(hotspot_icode)}]"
        with open(rfab_icode_path, 'w') as f:
            f.write(rfab_icode_str + "\n")
        recommended = "with_icode"
        logger.info(f"  ⚠ 检测到 insertion code，生成两份 hotspots snippet")
        logger.info(f"    - hotspots_rfantibody.txt (digits_only)")
        logger.info(f"    - hotspots_rfantibody_with_icode.txt (推荐)")
    else:
        logger.info(f"  写入 hotspots.txt 和 hotspots_rfantibody.txt")
    
    return rfab_digits_str


def write_design_loops(design_loops_str: str, output_dir: str, logger: logging.Logger):
    """写入 design_loops_rfantibody.txt"""
    loops_dir = os.path.join(output_dir, "loops")
    os.makedirs(loops_dir, exist_ok=True)
    
    output_path = os.path.join(loops_dir, "design_loops_rfantibody.txt")
    with open(output_path, 'w') as f:
        f.write(design_loops_str + "\n")
    
    logger.info(f"  写入 design_loops_rfantibody.txt")


def write_cdr_boundaries(cdr_boundaries: list[CDRBoundary], output_dir: str, logger: logging.Logger):
    """写入 cdr_boundaries.json"""
    loops_dir = os.path.join(output_dir, "loops")
    os.makedirs(loops_dir, exist_ok=True)
    
    output_path = os.path.join(loops_dir, "cdr_boundaries.json")
    data = []
    for cdr in cdr_boundaries:
        data.append({
            'cdr_name': cdr.cdr_name,
            'chothia_start': cdr.chothia_start,
            'chothia_end': cdr.chothia_end,
            'abs_start': cdr.abs_start,
            'abs_end': cdr.abs_end,
            'length': cdr.length
        })
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    logger.info(f"  写入 cdr_boundaries.json")


def write_sequences(structure: AtomArray, assignment: RoleAssignment, output_dir: str, logger: logging.Logger):
    """写入序列FASTA文件"""
    seq_dir = os.path.join(output_dir, "sequences")
    os.makedirs(seq_dir, exist_ok=True)
    
    sequences = extract_chain_sequences(structure)
    
    # 重链
    if assignment.heavy_chain:
        with open(os.path.join(seq_dir, "H.fasta"), 'w') as f:
            f.write(f">H|{assignment.heavy_chain}\n")
            f.write(sequences.get(assignment.heavy_chain, "") + "\n")
    
    # 轻链
    if assignment.light_chain:
        with open(os.path.join(seq_dir, "L.fasta"), 'w') as f:
            f.write(f">L|{assignment.light_chain}\n")
            f.write(sequences.get(assignment.light_chain, "") + "\n")
    
    # 抗原
    with open(os.path.join(seq_dir, "target.fasta"), 'w') as f:
        for chain in assignment.antigen_chains:
            f.write(f">T|{chain}\n")
            f.write(sequences.get(chain, "") + "\n")
    
    logger.info(f"  写入序列文件")


# =============================================================================
# QC报告
# =============================================================================

# =============================================================================
# HLT Validation（真正的合规性检查）
# =============================================================================

def validate_hlt_compliance(
    framework_path: str,
    target_path: str,
    hotspots_path: str,
    assignment: RoleAssignment,
    logger: logging.Logger
) -> dict:
    """
    严格的 HLT 合规性检查 - fail-fast
    
    检查：
    1. framework链顺序必须 H → L（或只有H）
    2. REMARK格式正确且abs_index在有效范围内
    3. hotspots每个残基都能在target中定位
    """
    errors = []
    warnings = []
    
    logger.info("执行 HLT 合规性验证...")
    
    # 1. 检查framework链顺序
    try:
        pdb_file = PDBFile.read(framework_path)
        framework = pdb_file.get_structure()
        
        # 提取唯一链ID（保持出现顺序）
        seen = set()
        chains_in_order = []
        for cid in framework.chain_id:
            cid_str = str(cid)
            if cid_str not in seen:
                chains_in_order.append(cid_str)
                seen.add(cid_str)
        
        # 期望顺序
        if assignment.mode in ['antibody', 'scfv']:
            expected = ['H', 'L']
        else:
            expected = ['H']
        
        if chains_in_order != expected:
            errors.append(f"链顺序错误: {chains_in_order}，期望: {expected}")
            logger.error(f"  ❌ 链顺序: {chains_in_order} != {expected}")
        else:
            logger.info(f"  ✓ 链顺序正确: {chains_in_order}")
    
    except Exception as e:
        errors.append(f"读取framework_HLT.pdb失败: {e}")
    
    # 2. 解析并验证REMARK
    try:
        with open(framework_path, 'r') as f:
            content = f.read()
        
        remark_pattern = r'REMARK PDBinfo-LABEL:\s+(\d+)\s+([HL]\d)'
        remarks = re.findall(remark_pattern, content)
        
        if not remarks:
            errors.append("缺少 REMARK PDBinfo-LABEL CDR 标注")
            logger.error("  ❌ 缺少 CDR REMARK")
        else:
            # 获取残基总数
            n_residues = len(set(framework.res_id))
            
            # 检查每个REMARK的abs_index
            invalid_remarks = []
            for abs_idx_str, cdr_name in remarks:
                abs_idx = int(abs_idx_str)
                if abs_idx < 1 or abs_idx > n_residues:
                    invalid_remarks.append((abs_idx, cdr_name))
            
            if invalid_remarks:
                errors.append(f"REMARK abs_index 超出范围 [1, {n_residues}]: {invalid_remarks}")
                logger.error(f"  ❌ REMARK 索引超范围: {len(invalid_remarks)} 个")
            else:
                logger.info(f"  ✓ REMARK 格式正确: {len(remarks)} 个标注")
    
    except Exception as e:
        errors.append(f"解析REMARK失败: {e}")
    
    # 3. 验证hotspots可定位
    try:
        with open(hotspots_path, 'r') as f:
            hs_content = f.read()
        
        hs_match = re.search(r'ppi\.hotspot_res=\[(.*?)\]', hs_content)
        if hs_match:
            hotspots_str = hs_match.group(1)
            if hotspots_str.strip():
                hotspots_list = [h.strip() for h in hotspots_str.split(',')]
                
                target_pdb = PDBFile.read(target_path)
                target = target_pdb.get_structure()
                target_resids = set(target.res_id)
                
                unlocatable = []
                for hs in hotspots_list:
                    # 解析 T123 或 T123A
                    m = re.match(r'T(\d+)([A-Z])?', hs)
                    if not m:
                        unlocatable.append(f"{hs} (格式错误)")
                        continue
                    
                    resnum = int(m.group(1))
                    if resnum not in target_resids:
                        unlocatable.append(f"{hs} (不存在)")
                
                if unlocatable:
                    errors.append(f"hotspots 无法在 target 中定位: {unlocatable}")
                    logger.error(f"  ❌ hotspots 定位失败: {len(unlocatable)} 个")
                else:
                    logger.info(f"  ✓ hotspots 可定位: {len(hotspots_list)} 个")
            else:
                warnings.append("hotspots 列表为空")
                logger.warning("  ⚠ hotspots 列表为空")
    
    except Exception as e:
        errors.append(f"验证hotspots失败: {e}")
    
    # 构建结果
    validation = {
        'chain_order_ok': 'chain_order' not in str(errors).lower(),
        'cdr_remark_format_ok': 'remark' not in str(errors).lower(),
        'hotspots_locatable': 'hotspot' not in str(errors).lower() and 'unlocatable' not in str(errors).lower(),
        'errors': errors,
        'warnings': warnings,
        'passed': len(errors) == 0
    }
    
    # Fail-fast：如果有错误，立即失败
    if errors:
        logger.error("HLT 合规性验证失败！")
        for err in errors:
            logger.error(f"  - {err}")
        raise ValueError(f"HLT 合规性验证失败: {len(errors)} 个错误")
    
    if warnings:
        for warn in warnings:
            logger.warning(f"  - {warn}")
    
    logger.info("✓ HLT 合规性验证通过")
    return validation


# =============================================================================
# QC Report
# =============================================================================

def generate_qc_report(
    config: Config,
    assignment: RoleAssignment,
    cdr_boundaries: list[CDRBoundary],
    epitope_patch: list[ResidueTag],
    hotspots: list[HotspotResidue],
    hotspots_str: str,
    design_loops_str: str,
    logger: logging.Logger,
    # V2 新增参数
    fasta_records: list[FastaRecord] | None = None,
    chain_mappings: dict[str, ChainMapping] | None = None,
    resnum_conflict_info: dict | None = None
) -> QCReport:
    """生成QC报告 (V2 增强)"""
    
    report = QCReport()
    
    # 获取抗体类型描述
    mode_descriptions = {
        'antibody': 'Fab/IgG (配对 VH-VL)',
        'nanobody': 'VHH/Nanobody (仅重链)',
        'scfv': 'scFv (单链抗体 VH-linker-VL)'
    }
    detected_mode = assignment.mode
    antibody_type_desc = mode_descriptions.get(detected_mode, detected_mode)
    
    # Meta
    report.meta = {
        'pdbid': config.pdbid,
        'input_file': config.input_file,
        'numbering': config.numbering,
        'mode': config.mode,
        'detected_mode': detected_mode,  # V2: 实际检测到的模式
        'antibody_type': antibody_type_desc,  # V2: 人类可读描述
        'timestamp': datetime.now().isoformat(),
        'hotspots_method': config.hotspots_method,
        'interface_cutoff': config.interface_cutoff,
        # V2
        'version': 'v2',
        'fasta_file': config.fasta_file if config.fasta_file else None,
        'fasta_mode': config.fasta_mode
    }
    
    # Chain assignment
    report.chain_assignment = {
        'heavy_chain': assignment.heavy_chain,
        'light_chain': assignment.light_chain,
        'antigen_chains': assignment.antigen_chains,
        'mode': assignment.mode,
        'detected_mode': detected_mode,  # V2: 重复便于单独解析
        'antibody_type': antibody_type_desc,
        'vh_vl_contact_score': assignment.vh_vl_contact_score,
        'all_heavy_chains': assignment.all_heavy_chains,
        'all_light_chains': assignment.all_light_chains
    }
    
    # CDR info
    report.cdr_info = {
        'boundaries': [
            {
                'cdr_name': c.cdr_name,
                'abs_start': c.abs_start,
                'abs_end': c.abs_end,
                'length': c.length
            }
            for c in cdr_boundaries
        ]
    }
    
    # Interface info
    report.interface_info = {
        'epitope_count': len(epitope_patch),
        'epitope_residues': [t.to_pdb_format() for t in epitope_patch[:20]]  # 前20个示例
    }
    
    # Hotspots info
    report.hotspots_info = {
        'method': config.hotspots_method,
        'k': config.hotspots_k,
        'selected_count': len(hotspots),
        'hotspots': [hs.to_dict() for hs in hotspots]
    }
    
    # RFantibody snippets
    report.rfantibody_snippets = {
        'ppi.hotspot_res': hotspots_str,
        'antibody.design_loops': design_loops_str
    }
    
    # HLT validation
    report.hlt_validation = {
        'framework_chain_ids': ['H', 'L'] if assignment.mode == 'antibody' else ['H'],
        'target_chain_id': 'T',
        'chain_order_ok': True,
        'cdr_remark_format_ok': True
    }
    
    # ============ V2 新增字段 ============
    
    # FASTA 信息
    if fasta_records:
        report.fasta_info = {
            'fasta_file': config.fasta_file,
            'record_count': len(fasta_records),
            'records': [
                {
                    'id': r.id,
                    'chain_ids': r.chain_ids,
                    'role_hints': r.role_hints,
                    'confidence': r.confidence
                }
                for r in fasta_records
            ]
        }
    
    # 链映射信息
    if chain_mappings:
        report.chain_mapping_info = {
            'mode': config.fasta_mode,
            'mappings': [
                {
                    'chain_id': m.struct_chain_id,
                    'fasta_id': m.fasta_id,
                    'identity': round(m.identity, 3),
                    'coverage': round(m.coverage_struct, 3),
                    'is_ambiguous': m.is_ambiguous
                }
                for m in chain_mappings.values()
            ]
        }
    
    # resnum 冲突信息
    if resnum_conflict_info:
        report.resnum_conflict_info = resnum_conflict_info
    
    return report


def write_qc_report(report: QCReport, reports_dir: str, logger: logging.Logger):
    """写入QC报告"""
    os.makedirs(reports_dir, exist_ok=True)
    
    # JSON格式
    json_path = os.path.join(reports_dir, "qc_report.json")
    with open(json_path, 'w') as f:
        json.dump(asdict(report), f, indent=2, default=str)
    
    # Markdown格式
    md_path = os.path.join(reports_dir, "qc_report.md")
    with open(md_path, 'w') as f:
        f.write("# Crystal2HLT-RFInputs QC Report\n\n")
        
        # 检测模式 section (新增，放在最前面)
        detected_mode = report.meta.get('detected_mode', 'N/A')
        antibody_type = report.meta.get('antibody_type', 'N/A')
        f.write(f"## 检测到的抗体类型\n")
        f.write(f"- **模式**: `{detected_mode}`\n")
        f.write(f"- **类型描述**: {antibody_type}\n\n")
        
        f.write(f"## Meta\n")
        f.write(f"- **PDBID**: {report.meta.get('pdbid', 'N/A')}\n")
        f.write(f"- **Input Mode**: {report.meta.get('mode', 'N/A')}\n")
        f.write(f"- **Numbering**: {report.meta.get('numbering', 'N/A')}\n")
        f.write(f"- **Timestamp**: {report.meta.get('timestamp', 'N/A')}\n\n")
        
        f.write(f"## Chain Assignment\n")
        f.write(f"- **Heavy Chain**: {report.chain_assignment.get('heavy_chain', 'N/A')}\n")
        f.write(f"- **Light Chain**: {report.chain_assignment.get('light_chain', 'N/A')}\n")
        f.write(f"- **Antigen Chains**: {report.chain_assignment.get('antigen_chains', [])}\n\n")
        
        f.write(f"## CDR Boundaries\n")
        for cdr in report.cdr_info.get('boundaries', []):
            f.write(f"- **{cdr['cdr_name']}**: abs {cdr['abs_start']}-{cdr['abs_end']} (length={cdr['length']})\n")
        f.write("\n")
        
        f.write(f"## Interface\n")
        f.write(f"- **Epitope residues**: {report.interface_info.get('epitope_count', 0)}\n\n")
        
        f.write(f"## Hotspots\n")
        f.write(f"- **Method**: {report.hotspots_info.get('method', 'N/A')}\n")
        f.write(f"- **Selected**: {report.hotspots_info.get('selected_count', 0)}\n\n")
        
        f.write(f"## RFantibody Snippets\n")
        f.write(f"```\n{report.rfantibody_snippets.get('ppi.hotspot_res', '')}\n```\n\n")
        f.write(f"```\n{report.rfantibody_snippets.get('antibody.design_loops', '')}\n```\n\n")
        
        if report.errors:
            f.write(f"## Errors\n")
            for err in report.errors:
                f.write(f"- {err}\n")
        
        if report.warnings:
            f.write(f"## Warnings\n")
            for warn in report.warnings:
                f.write(f"- {warn}\n")
    
    logger.info(f"  写入 qc_report.json 和 qc_report.md")


def write_error_report(
    config: Config, 
    error_message: str, 
    error_type: str,
    assignment: Optional[RoleAssignment],
    hetatm_info: Optional[Dict] = None,
    logger: Optional[logging.Logger] = None
):
    """当管线失败时写入错误报告，提供诊断信息"""
    from datetime import datetime
    
    reports_dir = os.path.join(config.outdir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    
    # 构建错误报告
    error_report = {
        "status": "FAILED",
        "error_type": error_type,
        "error_message": error_message,
        "timestamp": datetime.now().isoformat(),
        "input_file": config.input_file,
        "pdbid": config.pdbid,
    }
    
    # 添加链分配信息（如果有）
    if assignment:
        error_report["chain_assignment"] = {
            "heavy_chain": assignment.heavy_chain,
            "light_chain": assignment.light_chain,
            "antigen_chains": assignment.antigen_chains,
        }
    
    # 添加 HETATM 信息（如果检测到）
    if hetatm_info and hetatm_info.get('has_hetatm_interface'):
        error_report["hetatm_interface"] = {
            "detected": True,
            "types": hetatm_info.get('hetatm_types', {}),
            "residues_sample": hetatm_info.get('hetatm_residues', [])[:10],
            "suggestion": "Crystal2HLT 目前仅支持蛋白-蛋白界面。非蛋白抗原(小分子/糖/核酸)暂不支持。"
        }
    
    # 写入 JSON
    json_path = os.path.join(reports_dir, "error_report.json")
    with open(json_path, 'w') as f:
        json.dump(error_report, f, indent=2, ensure_ascii=False)
    
    # 写入 Markdown
    md_path = os.path.join(reports_dir, "error_report.md")
    with open(md_path, 'w') as f:
        f.write("# Crystal2HLT-RFInputs Error Report\n\n")
        f.write("> [!CAUTION]\n")
        f.write(f"> 管线执行失败: {error_type}\n\n")
        
        f.write("## 基本信息\n")
        f.write(f"- **PDB ID**: {config.pdbid}\n")
        f.write(f"- **输入文件**: `{config.input_file}`\n")
        f.write(f"- **时间戳**: {error_report['timestamp']}\n\n")
        
        f.write("## 错误详情\n")
        f.write(f"```\n{error_message}\n```\n\n")
        
        if assignment:
            f.write("## 链分配 (执行到此步骤)\n")
            f.write(f"- **重链 (H)**: {assignment.heavy_chain}\n")
            f.write(f"- **轻链 (L)**: {assignment.light_chain or 'N/A'}\n")
            f.write(f"- **抗原链 (T)**: {assignment.antigen_chains}\n\n")
        
        if hetatm_info and hetatm_info.get('has_hetatm_interface'):
            f.write("## 检测到的非蛋白抗原\n")
            types = hetatm_info.get('hetatm_types', {})
            f.write(f"- **小分子配体**: {types.get('ligand', 0)}\n")
            f.write(f"- **糖类残基**: {types.get('sugar', 0)}\n")
            f.write(f"- **核酸残基**: {types.get('nucleic', 0)}\n")
            f.write(f"- **金属离子**: {types.get('metal', 0)}\n\n")
            
            residues = hetatm_info.get('hetatm_residues', [])[:5]
            if residues:
                f.write("**示例残基**:\n")
                for res in residues:
                    f.write(f"- `{res}`\n")
            f.write("\n")
            
            f.write("> [!IMPORTANT]\n")
            f.write("> Crystal2HLT 目前仅支持蛋白-蛋白界面。\n")
            f.write("> 对于抗半抗原(小分子)、抗糖或抗核酸的抗体，暂无法处理。\n\n")
        
        f.write("## 建议\n")
        if hetatm_info and hetatm_info.get('has_hetatm_interface'):
            f.write("1. 确认此抗体的抗原确实是非蛋白分子\n")
            f.write("2. 如需设计针对非蛋白抗原的抗体，可能需要手动准备输入文件\n")
            f.write("3. 或考虑使用其他工具处理小分子/糖/核酸抗原\n")
        else:
            f.write("1. 检查输入结构是否正确包含抗体-抗原复合物\n")
            f.write("2. 尝试使用 `--bioassembly True` 加载生物学组装体\n")
            f.write("3. 手动指定链角色: `--heavy_chain` / `--light_chain` / `--target_chains`\n")
    
    if logger:
        logger.info(f"  写入错误报告: {json_path}")
        logger.info(f"  写入错误报告: {md_path}")
    
    return error_report


# =============================================================================
# 主流程
# =============================================================================

def setup_directories(config: Config):
    """创建输出目录结构"""
    config.input_dir = os.path.join(config.outdir, "input")
    config.work_dir = os.path.join(config.outdir, "work")
    config.logs_dir = os.path.join(config.outdir, "logs")
    # 注意: 不再创建主目录的 outputs/ 和 reports/
    # 所有输出现在都在 antibody_N/ 子目录中
    
    for d in [config.input_dir, config.work_dir, config.logs_dir]:
        os.makedirs(d, exist_ok=True)


def copy_input_file(config: Config, logger: logging.Logger):
    """复制输入文件到input目录"""
    input_basename = os.path.basename(config.input_file)
    dest = os.path.join(config.input_dir, input_basename)
    shutil.copy2(config.input_file, dest)
    logger.info(f"复制输入文件到 {dest}")


def process_antibody_output(
    structure: AtomArray,
    pairing: AntibodyPairing,
    all_assignment: RoleAssignment,
    config: Config,
    fasta_records,
    chain_mappings,
    logger: logging.Logger,
    is_symmetric_copy: bool = False
) -> dict:
    """
    为单个抗体生成完整的独立输出目录
    
    返回:
        {
            'index': pairing.index,
            'mode': pairing.mode,
            'heavy_chain': pairing.heavy_chain,
            'light_chain': pairing.light_chain,
            'epitope_count': int,
            'hotspots_count': int,
            'is_symmetric_copy': bool,
            'output_dir': str
        }
    """
    # 创建抗体目录结构
    antibody_dir = os.path.join(config.outdir, f"antibody_{pairing.index}")
    output_dir = os.path.join(antibody_dir, "outputs")
    reports_dir = os.path.join(antibody_dir, "reports")
    epitopes_dir = os.path.join(output_dir, "epitopes")
    loops_dir = os.path.join(output_dir, "loops")
    sequences_dir = os.path.join(output_dir, "sequences")
    
    for d in [output_dir, reports_dir, epitopes_dir, loops_dir, sequences_dir]:
        os.makedirs(d, exist_ok=True)
    
    logger.info(f"\n--- 抗体 {pairing.index}: {pairing.heavy_chain}-{pairing.light_chain or 'None'} ({pairing.mode}) ---")
    
    # === 确定该抗体结合的抗原链 ===
    # 计算该抗体与每个抗原链的接触分数，只选择有显著接触的链
    ab_chains = [pairing.heavy_chain]
    if pairing.light_chain:
        ab_chains.append(pairing.light_chain)
    
    # 计算与每个抗原链的接触
    antigen_contacts = []
    for ag_chain in all_assignment.antigen_chains:
        # 计算接触分数
        contact_score = compute_chain_contact_score(structure, pairing.heavy_chain, ag_chain, cutoff=8.0)
        if pairing.light_chain:
            contact_score += compute_chain_contact_score(structure, pairing.light_chain, ag_chain, cutoff=8.0)
        antigen_contacts.append((ag_chain, contact_score))
    
    # 按接触分数排序，选择有显著接触的抗原链（至少有一些接触）
    antigen_contacts.sort(key=lambda x: -x[1])
    
    # 筛选有显著接触的抗原链（接触分数 > 0）
    bound_antigen_chains = [ag for ag, score in antigen_contacts if score > 0]
    
    if bound_antigen_chains:
        logger.info(f"  该抗体结合的抗原链: {bound_antigen_chains} (从 {len(all_assignment.antigen_chains)} 条中选出)")
    else:
        # 如果没有找到任何接触，使用所有抗原链
        bound_antigen_chains = all_assignment.antigen_chains
        logger.warning(f"  未检测到抗原接触，使用所有抗原链: {bound_antigen_chains}")
    
    # 创建此抗体的 RoleAssignment（使用筛选后的抗原链）
    sub_assignment = RoleAssignment(
        heavy_chain=pairing.heavy_chain,
        light_chain=pairing.light_chain,
        antigen_chains=bound_antigen_chains,  # 使用筛选后的抗原链
        mode=pairing.mode,
        vh_vl_contact_score=pairing.contact_score,
        all_heavy_chains=[pairing.heavy_chain],
        all_light_chains=[pairing.light_chain] if pairing.light_chain else [],
        anarci_results={
            k: v for k, v in all_assignment.anarci_results.items()
            if k in [pairing.heavy_chain, pairing.light_chain]
        }
    )
    
    # 创建此抗体的 Config
    sub_config = Config(
        input_file=config.input_file,
        outdir=antibody_dir,
        pdbid=config.pdbid,
        mode=pairing.mode,
        numbering=config.numbering,
        hotspots_method=config.hotspots_method,
        hotspots_k=config.hotspots_k,
        interface_cutoff=config.interface_cutoff,
        h_crop=config.h_crop,
        l_crop=config.l_crop
    )
    sub_config.output_dir = output_dir
    sub_config.reports_dir = reports_dir
    
    # === 导出 framework_HLT.pdb ===
    framework_path, cdr_bounds = export_framework_hlt(
        structure, sub_assignment, sub_config, logger
    )
    
    # === resnum 冲突检测与处理 ===
    resnum_conflict_info = None
    resnum_mapping = {}
    if len(sub_assignment.antigen_chains) > 1:
        resnum_conflict_info = detect_target_resnum_conflicts(
            structure, sub_assignment.antigen_chains, logger
        )
    
    # === 导出 target_HLT.pdb ===
    target_path, resnum_mapping = export_target_hlt(
        structure, sub_assignment, sub_config, logger, resnum_conflict_info
    )
    
    # === 计算 epitope ===
    sub_epitope = compute_epitope_patch(structure, sub_assignment, sub_config, logger)
    
    # 如果有重编号，更新 epitope
    if resnum_conflict_info and resnum_conflict_info.get('has_conflict'):
        sub_epitope = update_residue_numbers(sub_epitope, resnum_mapping)
    
    result = {
        'index': pairing.index,
        'mode': pairing.mode,
        'heavy_chain': pairing.heavy_chain,
        'light_chain': pairing.light_chain,
        'epitope_count': len(sub_epitope),
        'hotspots_count': 0,
        'is_symmetric_copy': is_symmetric_copy,
        'output_dir': antibody_dir
    }
    
    if not sub_epitope:
        logger.warning(f"  抗体 {pairing.index}: 未找到界面残基")
        return result
    
    # === 选择 hotspots ===
    sub_hotspots = select_hotspots(
        structure, sub_epitope, sub_assignment, sub_config, logger, cdr_bounds
    )
    result['hotspots_count'] = len(sub_hotspots)
    
    # === 生成 design_loops ===
    design_loops_str = build_design_loops(cdr_bounds, sub_config, logger)
    
    # === 写入输出文件 ===
    write_epitope_patch(sub_epitope, output_dir, logger)
    hotspots_str = write_hotspots(sub_hotspots, output_dir, logger)
    write_design_loops(design_loops_str, output_dir, logger)
    write_cdr_boundaries(cdr_bounds, output_dir, logger)
    write_sequences(structure, sub_assignment, output_dir, logger)
    
    # === HLT 验证 ===
    framework_hlt_path = os.path.join(output_dir, "framework_HLT.pdb")
    target_hlt_path = os.path.join(output_dir, "target_HLT.pdb")
    hotspots_rfa_path = os.path.join(epitopes_dir, "hotspots_rfantibody.txt")
    
    hlt_validation = validate_hlt_compliance(
        framework_hlt_path, target_hlt_path, hotspots_rfa_path,
        sub_assignment, logger
    )
    
    # === 生成 QC 报告 ===
    qc_report = generate_qc_report(
        sub_config, sub_assignment, cdr_bounds, sub_epitope,
        sub_hotspots, hotspots_str, design_loops_str, logger,
        fasta_records=fasta_records,
        chain_mappings=chain_mappings,
        resnum_conflict_info=resnum_conflict_info
    )
    qc_report.hlt_validation = hlt_validation
    
    # 添加对称拷贝标记
    if is_symmetric_copy:
        qc_report.meta['is_symmetric_copy'] = True
        qc_report.meta['symmetric_note'] = '此抗体是晶体对称相关拷贝，与其他抗体具有相同序列和/或epitope'
    
    write_qc_report(qc_report, reports_dir, logger)
    
    logger.info(f"  ✓ 抗体 {pairing.index}: {len(sub_epitope)} epitope, {len(sub_hotspots)} hotspots")
    
    return result


def run_pipeline(config: Config) -> QCReport:
    """运行完整管线"""
    
    # 设置日志
    log_file = os.path.join(config.outdir, "logs", "run.log") if config.verbose else None
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = setup_logging(config.verbose, log_file)
    
    logger.info("=" * 60)
    logger.info("Crystal2HLT-RFInputs Pipeline")
    logger.info("=" * 60)
    
    # 设置目录
    setup_directories(config)
    
    # 复制输入
    copy_input_file(config, logger)
    
    # Module 0: 读取结构（先用ASU）
    input_path = config.input_file
    file_format = detect_file_format(input_path)
    structure_asu = normalize_structure(input_path, config, logger)
    
    # Module 1: 清理ASU
    structure_asu_cleaned = clean_structure(structure_asu, config, logger)
    
    # ============ V2 FASTA 解析与映射 ============
    fasta_records = None
    chain_mappings = None
    chain_seq_info = None
    
    if config.fasta_file and config.fasta_mode != 'off':
        logger.info("======== V2 FASTA 辅助识别 ========")
        
        # 解析 FASTA
        fasta_records = parse_fasta(config.fasta_file, logger)
        
        # 提取结构链序列信息
        chain_seq_info = extract_chain_sequences_v2(structure_asu_cleaned, logger)
        
        # 映射结构链到 FASTA
        chain_mappings = map_structure_chains_to_fasta(
            chain_seq_info, fasta_records, config.fasta_mode, logger
        )
        
        logger.info("====================================")
    
    # Module 2: 链角色识别（ASU）- 使用 FASTA 辅助
    assignment_asu = assign_roles_v2(
        structure_asu_cleaned, config, 
        fasta_records, chain_mappings, logger
    )
    
    # Module 5: 先用ASU计算界面
    epitope_patch_asu = compute_epitope_patch(structure_asu_cleaned, assignment_asu, config, logger)
    
    # ============ B) Biological Assembly 自动切换 ============
    used_assembly = False
    assembly_id = None
    
    if len(epitope_patch_asu) == 0 and config.use_bioassembly in ["auto", "true"]:
        logger.warning("ASU中未找到界面残基，尝试使用 biological assembly...")
        
        if file_format == 'cif':
            try:
                # 读取CIF并提取assemblies
                from biotite.structure.io.pdbx import CIFFile
                cif_file = CIFFile.read(input_path)
                assemblies = extract_mmcif_assemblies(cif_file, logger)
                
                if assemblies:
                    # 使用指定的assembly ID
                    assembly_id = str(config.bioassembly_id)
                    if assembly_id not in assemblies:
                        logger.warning(f"指定的 assembly {assembly_id} 不存在，使用第一个可用的")
                        assembly_id = list(assemblies.keys())[0]
                    
                    logger.info(f"构建 biological assembly {assembly_id}...")
                    structure = build_bioassembly(cif_file, assembly_id, logger)
                    
                    # 重新清理和识别
                    structure = clean_structure(structure, config, logger)
                    assignment = assign_roles(structure, config, logger)
                    
                    # 重新计算epitope
                    epitope_patch = compute_epitope_patch(structure, assignment, config, logger)
                    
                    if len(epitope_patch) > 0:
                        logger.info(f"✓ Assembly {assembly_id} 找到 {len(epitope_patch)} 个界面残基")
                        used_assembly = True
                    else:
                        logger.warning(f"Assembly {assembly_id} 仍无界面残基")
                        
                        # 检测是否有抗原链
                        if not assignment.antigen_chains:
                            # P2.1: 检测 HETATM 配体界面
                            ab_chains = [assignment.heavy_chain]
                            if assignment.light_chain:
                                ab_chains.append(assignment.light_chain)
                            
                            hetatm_info = detect_hetatm_interface(
                                config.input_file, ab_chains, config.interface_cutoff, logger
                            )
                            
                            if hetatm_info['has_hetatm_interface']:
                                types = hetatm_info['hetatm_types']
                                error_msg = (f"无蛋白抗原链，但检测到非蛋白抗原 (HETATM): {types}. "
                                           "Crystal2HLT 仅支持蛋白-蛋白界面。")
                                write_error_report(
                                    config=config,
                                    error_message=error_msg,
                                    error_type="NON_PROTEIN_ANTIGEN",
                                    assignment=assignment,
                                    hetatm_info=hetatm_info,
                                    logger=logger
                                )
                                raise ValueError(error_msg)
                            else:
                                error_msg = "无抗原链，且未检测到 HETATM 配体界面"
                                write_error_report(
                                    config=config,
                                    error_message=error_msg,
                                    error_type="NO_ANTIGEN",
                                    assignment=assignment,
                                    hetatm_info=None,
                                    logger=logger
                                )
                                raise ValueError(error_msg)
                        else:
                            error_msg = "即使使用 biological assembly 仍未找到界面残基"
                            write_error_report(
                                config=config,
                                error_message=error_msg,
                                error_type="NO_INTERFACE",
                                assignment=assignment,
                                hetatm_info=None,
                                logger=logger
                            )
                            raise ValueError(error_msg)
                else:
                    error_msg = "ASU无接触且无可用assembly"
                    write_error_report(
                        config=config,
                        error_message=error_msg,
                        error_type="NO_ASSEMBLY",
                        assignment=assignment_asu,
                        hetatm_info=None,
                        logger=logger
                    )
                    logger.warning("CIF 文件中无可用的 biological assembly")
                    raise ValueError(error_msg)
            
            except Exception as e:
                logger.error(f"构建 biological assembly 失败: {e}")
                # 检查是否已生成错误报告（如果是我们抛出的 ValueError）
                reports_dir = os.path.join(config.outdir, "reports")
                error_report_path = os.path.join(reports_dir, "error_report.json")
                if not os.path.exists(error_report_path):
                    error_msg = f"ASU无接触，且assembly构建失败: {e}"
                    write_error_report(
                        config=config,
                        error_message=error_msg,
                        error_type="ASSEMBLY_FAILED",
                        assignment=assignment_asu,
                        hetatm_info=None,
                        logger=logger
                    )
                raise ValueError(f"ASU无接触，且assembly构建失败: {e}")
        else:
            error_msg = "ASU无接触，PDB格式不支持assembly"
            write_error_report(
                config=config,
                error_message=error_msg,
                error_type="PDB_NO_ASSEMBLY",
                assignment=assignment_asu,
                hetatm_info=None,
                logger=logger
            )
            logger.error("PDB格式不支持自动 assembly 构建，请使用 mmCIF 格式")
            raise ValueError(error_msg)
    
    else:
        # 使用ASU
        structure = structure_asu_cleaned
        assignment = assignment_asu
        epitope_patch = epitope_patch_asu
        logger.info(f"使用 ASU，找到 {len(epitope_patch)} 个界面残基")
    
    # 最终检查
    if len(epitope_patch) == 0:
        logger.error("未找到界面残基！")
        
        # P2.1: 检测 HETATM 配体界面
        ab_chains = [assignment.heavy_chain]
        if assignment.light_chain:
            ab_chains.append(assignment.light_chain)
        
        hetatm_info = detect_hetatm_interface(
            config.input_file, ab_chains, config.interface_cutoff, logger
        )
        
        if hetatm_info['has_hetatm_interface']:
            types = hetatm_info['hetatm_types']
            residues = hetatm_info['hetatm_residues'][:5]
            
            logger.warning("=" * 60)
            logger.warning("检测到 非蛋白抗原 (HETATM) 界面!")
            logger.warning(f"  配体类型统计: 小分子={types.get('ligand', 0)} 糖={types.get('sugar', 0)} "
                          f"核酸={types.get('nucleic', 0)} 金属={types.get('metal', 0)}")
            logger.warning(f"  示例残基: {residues}")
            logger.warning("")
            logger.warning("Crystal2HLT 目前仅支持蛋白-蛋白界面。")
            logger.warning("非蛋白抗原 (半抗原/糖/核酸) 暂不支持。")
            logger.warning("=" * 60)
            
            # 生成错误报告
            error_msg = f"检测到非蛋白抗原: {types}. Crystal2HLT 仅支持蛋白-蛋白界面。"
            write_error_report(
                config=config,
                error_message=error_msg,
                error_type="NON_PROTEIN_ANTIGEN",
                assignment=assignment,
                hetatm_info=hetatm_info,
                logger=logger
            )
            
            raise ValueError(error_msg)
        else:
            # 生成错误报告（无界面）
            error_msg = "Epitope patch为空，请检查输入结构或使用 biological assembly"
            write_error_report(
                config=config,
                error_message=error_msg,
                error_type="NO_INTERFACE",
                assignment=assignment,
                hetatm_info=None,
                logger=logger
            )
            raise ValueError(error_msg)
    
    # 保存cleaned结构
    cleaned_path = os.path.join(config.work_dir, "cleaned.pdb")
    pdb_file = PDBFile()
    pdb_file.set_structure(structure)
    pdb_file.write(cleaned_path)
    
    # ============ 抗体中心化输出 ============
    # 所有抗体都输出到独立的 antibody_N 目录
    logger.info("=" * 60)
    logger.info(f"为 {len(assignment.all_pairings)} 套抗体生成完整输出...")
    logger.info("=" * 60)
    
    # 检测对称拷贝 (相同序列的抗体)
    sequences_seen = {}
    symmetric_copies = set()
    for pairing in assignment.all_pairings:
        h_seq = extract_chain_sequences(structure).get(pairing.heavy_chain, "")
        l_seq = extract_chain_sequences(structure).get(pairing.light_chain, "") if pairing.light_chain else ""
        seq_key = f"{h_seq}|{l_seq}"
        
        if seq_key in sequences_seen:
            symmetric_copies.add(pairing.index)
            logger.info(f"  抗体 {pairing.index} 检测为对称相关拷贝 (与抗体 {sequences_seen[seq_key]} 相同序列)")
        else:
            sequences_seen[seq_key] = pairing.index
    
    # 处理每个抗体
    antibody_results = []
    for pairing in assignment.all_pairings:
        is_symmetric = pairing.index in symmetric_copies
        result = process_antibody_output(
            structure=structure,
            pairing=pairing,
            all_assignment=assignment,
            config=config,
            fasta_records=fasta_records,
            chain_mappings=chain_mappings,
            logger=logger,
            is_symmetric_copy=is_symmetric
        )
        antibody_results.append(result)
    
    # ============ 生成全局摘要报告 ============
    logger.info("\n生成全局摘要报告...")
    summary = {
        'input_file': config.input_file,
        'pdbid': config.pdbid,
        'timestamp': datetime.now().isoformat(),
        'total_antibodies': len(antibody_results),
        'antibodies': antibody_results,
        'symmetric_copies_detected': len(symmetric_copies),
        'symmetric_copy_indices': list(symmetric_copies)
    }
    
    # 写入 summary_report.json
    summary_json_path = os.path.join(config.outdir, "summary_report.json")
    with open(summary_json_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    # 写入 summary_report.md
    summary_md_path = os.path.join(config.outdir, "summary_report.md")
    with open(summary_md_path, 'w') as f:
        f.write("# Crystal2HLT 处理摘要\n\n")
        f.write(f"- **输入文件**: {config.input_file}\n")
        f.write(f"- **PDB ID**: {config.pdbid or 'N/A'}\n")
        f.write(f"- **处理时间**: {summary['timestamp']}\n")
        f.write(f"- **检测到抗体数量**: {summary['total_antibodies']}\n\n")
        
        if symmetric_copies:
            f.write("## ⚠ 对称相关拷贝\n\n")
            f.write(f"检测到 {len(symmetric_copies)} 个对称相关拷贝: {list(symmetric_copies)}\n\n")
            f.write("这些抗体与其他抗体具有相同的序列，可能是晶体堆积产生的对称相关分子。\n\n")
        
        f.write("## 抗体列表\n\n")
        f.write("| # | 模式 | 重链 | 轻链 | Epitope | Hotspots | 对称拷贝 |\n")
        f.write("|---|------|------|------|---------|----------|----------|\n")
        for ab in antibody_results:
            sym_mark = "✓" if ab['is_symmetric_copy'] else ""
            f.write(f"| {ab['index']} | {ab['mode']} | {ab['heavy_chain']} | {ab['light_chain'] or '-'} | {ab['epitope_count']} | {ab['hotspots_count']} | {sym_mark} |\n")
        f.write("\n")
        
        f.write("## 输出目录\n\n")
        for ab in antibody_results:
            f.write(f"- `antibody_{ab['index']}/` - {ab['mode']}\n")
    
    logger.info(f"  写入 {summary_json_path}")
    logger.info(f"  写入 {summary_md_path}")
    
    logger.info("")
    logger.info("=" * 60)
    logger.info("Pipeline 完成！")
    logger.info(f"输出目录: {config.outdir}")
    logger.info(f"抗体目录: antibody_1, antibody_2, ... (共 {len(antibody_results)} 套)")
    logger.info("=" * 60)
    
    # 返回第一个抗体的 QC 报告作为主报告
    first_ab_qc_path = os.path.join(config.outdir, "antibody_1", "reports", "qc_report.json")
    if os.path.exists(first_ab_qc_path):
        with open(first_ab_qc_path, 'r') as f:
            qc_data = json.load(f)
        return QCReport(**{k: v for k, v in qc_data.items() if k in QCReport.__dataclass_fields__})
    
    return QCReport()


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Crystal2HLT-RFInputs: 共晶结构 → RFantibody HLT 输入自动构建',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 自动识别链角色
  python crystal2hlt_rfinputs.py -i complex.pdb -o output_dir

  # 手动指定链ID
  python crystal2hlt_rfinputs.py -i complex.pdb -o output_dir --heavy_chain H --light_chain L

  # 多个重链/轻链（选择接触最多的一对）
  python crystal2hlt_rfinputs.py -i complex.pdb -o output_dir --heavy_chain H,J --light_chain I,L
        """
    )
    
    # 必需参数
    parser.add_argument('--input', '-i', required=True, help='输入PDB/CIF文件')
    parser.add_argument('--outdir', '-o', required=True, help='输出目录')
    
    # 链ID指定（可选，但强烈建议）
    parser.add_argument('--heavy_chain', '-H', default='',
                        help='重链ID（逗号分隔多个）')
    parser.add_argument('--light_chain', '-L', default='',
                        help='轻链ID（逗号分隔多个）')
    parser.add_argument('--antigen_chains', '-A', default='',
                        help='抗原链ID（逗号分隔，留空则自动推断）')
    
    # 其他参数
    parser.add_argument('--pdbid', default='', help='PDB ID')
    parser.add_argument('--mode', choices=['antibody', 'nanobody', 'auto'], default='auto',
                        help='模式: antibody/nanobody/auto')
    parser.add_argument('--numbering', choices=['chothia', 'imgt', 'aho'], default='chothia',
                        help='编号方案')
    parser.add_argument('--interface_cutoff', type=float, default=4.5,
                        help='界面cutoff距离(Å)')
    parser.add_argument('--hotspots_k', type=int, default=12,
                        help='选择的热点数量')
    parser.add_argument('--hotspots_method', choices=['contact_score', 'rfantibody_cbeta5'],
                        default='contact_score', help='热点选择方法')
    parser.add_argument('--space_dedup_angstrom', type=float, default=6.0,
                        help='热点空间去重距离(Å)')
    parser.add_argument('--use_bioassembly', choices=['auto', 'true', 'false'], default='auto',
                        help='是否使用biological assembly')
    parser.add_argument('--bioassembly_id', type=int, default=1,
                        help='Biological assembly ID')
    
    # Boolean 参数（修复 type=bool 问题）
    parser.add_argument('--remove-waters', dest='remove_waters',
                       action='store_true', default=True,
                       help='移除水分子（默认）')
    parser.add_argument('--no-remove-waters', dest='remove_waters',
                       action='store_false',
                       help='保留水分子')
    
    parser.add_argument('--remove-hetero', dest='remove_hetero',
                       action='store_true', default=True,
                       help='移除非蛋白HETATM（默认）')
    parser.add_argument('--no-remove-hetero', dest='remove_hetero',
                       action='store_false',
                       help='保留非蛋白HETATM')
    
    parser.add_argument('--keep-glycans', dest='keep_glycans',
                       action='store_true', default=False,
                       help='保留糖链')
    parser.add_argument('--no-keep-glycans', dest='keep_glycans',
                       action='store_false',
                       help='移除糖链（默认）')
    parser.add_argument('--h_crop', type=int, default=115,
                        help='重链裁剪位置')
    parser.add_argument('--l_crop', type=int, default=110,
                        help='轻链裁剪位置')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='详细输出')
    
    # V2 新增参数: FASTA 输入支持
    parser.add_argument('--fasta', '-f', default='',
                        help='PDB FASTA 文件 (.fasta/.fa/.fsa)')
    parser.add_argument('--fasta_mode', choices=['auto', 'strict', 'off'],
                        default='auto',
                        help='FASTA 使用模式: auto(辅助)/strict(必须映射)/off(禁用)')
    parser.add_argument('--target_chains', '-T', default='',
                        help='手动指定 target 链ID（逗号分隔）')
    parser.add_argument('--prefer_cif_assembly', 
                        choices=['auto', 'asu', 'assembly'],
                        default='auto',
                        help='CIF assembly 偏好: auto/asu/assembly')
    parser.add_argument('--multi_antibody', 
                        choices=['off', 'top1', 'topN'],
                        default='top1',
                        help='多抗体处理模式: off/top1/topN')
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    config = Config(
        input_file=args.input,
        outdir=args.outdir,
        pdbid=args.pdbid,
        mode=args.mode,
        numbering=args.numbering,
        interface_cutoff=args.interface_cutoff,
        hotspots_k=args.hotspots_k,
        hotspots_method=args.hotspots_method,
        space_dedup_angstrom=args.space_dedup_angstrom,
        use_bioassembly=args.use_bioassembly,
        bioassembly_id=args.bioassembly_id,
        remove_waters=args.remove_waters,
        remove_hetero=args.remove_hetero,
        keep_glycans=args.keep_glycans,
        h_crop=args.h_crop,
        l_crop=args.l_crop,
        verbose=args.verbose,
        heavy_chain_id=args.heavy_chain,
        light_chain_id=args.light_chain,
        antigen_chain_ids=args.antigen_chains,
        # V2 新增参数
        fasta_file=args.fasta,
        fasta_mode=args.fasta_mode,
        target_chain_ids=args.target_chains,
        prefer_cif_assembly=args.prefer_cif_assembly,
        multi_antibody=args.multi_antibody
    )
    
    try:
        run_pipeline(config)
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if config.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
