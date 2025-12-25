"""
Crystal2HLT 配置与数据结构定义

包含所有核心 dataclass 和常量配置
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


# =============================================================================
# 常量定义
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

__version__ = "2.1.0"


# =============================================================================
# 核心数据结构
# =============================================================================

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
    """单个抗体配对信息"""
    heavy_chain: str               # VH 链 ID
    light_chain: str | None        # VL 链 ID，nanobody 为 None
    contact_score: float = 0.0     # VH-VL 接触分数
    index: int = 0                 # 配对索引 (1, 2, 3...)
    mode: Literal["antibody", "nanobody", "scfv"] = "antibody"


@dataclass
class RoleAssignment:
    """链角色分配结果"""
    heavy_chain: str | None        # 原始chain_id (主抗体)
    light_chain: str | None        # nanobody时为None
    antigen_chains: list[str]      # 可多条
    mode: Literal["antibody", "nanobody", "scfv"]
    vh_vl_contact_score: float = 0.0
    all_heavy_chains: list[str] = field(default_factory=list)
    all_light_chains: list[str] = field(default_factory=list)
    anarci_results: dict = field(default_factory=dict)
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
    """完整的 residue 映射"""
    abs_index: int              # 输出文件中的全局索引（1-indexed）
    auth_chain_id: str          # author 指定的链ID
    auth_seq_id: int            # author 序列编号
    pdbx_ins_code: str          # 插入码
    label_seq_id: int           # label 序列编号
    label_comp_id: str          # 残基名称（3字母码）
    label_asym_id: str          # label asymmetric unit ID
    chothia_num: int | None = None
    chothia_icode: str | None = None
    cdr_type: str | None = None


# =============================================================================
# FASTA 相关数据结构
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
    residue_list: list[tuple[int, str]] = field(default_factory=list)
    is_protein: bool = True
    length: int = 0
    filter_reason: str = ""
    
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
    domain_range: tuple[int, int] | None = None
    anarci_score: float = 0.0
    fasta_hint: str | None = None
    final_confidence: float = 0.0
    is_crystal_packing: bool = False


# =============================================================================
# 主配置
# =============================================================================

@dataclass
class Config:
    """全局配置"""
    input_file: str
    outdir: str
    pdbid: str = ""
    mode: Literal["antibody", "nanobody", "scfv", "auto"] = "auto"
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
    fail_on_no_anarci: bool = False
    verbose: bool = False
    
    # 手动链ID指定
    heavy_chain_id: str = ""
    light_chain_id: str = ""
    antigen_chain_ids: str = ""
    
    # FASTA 输入支持
    fasta_file: str = ""
    fasta_mode: Literal["auto", "strict", "off"] = "auto"
    target_chain_ids: str = ""
    prefer_cif_assembly: Literal["auto", "asu", "assembly"] = "auto"
    multi_antibody: Literal["off", "top1", "topN"] = "top1"
    
    # 内部使用
    work_dir: str = ""
    output_dir: str = ""
    input_dir: str = ""
    reports_dir: str = ""
    logs_dir: str = ""


@dataclass
class QCReport:
    """质量控制报告"""
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
    fasta_info: dict = field(default_factory=dict)
    chain_mapping_info: dict = field(default_factory=dict)
    chain_role_table: list = field(default_factory=list)
    target_selection_info: dict = field(default_factory=dict)
    resnum_conflict_info: dict = field(default_factory=dict)
