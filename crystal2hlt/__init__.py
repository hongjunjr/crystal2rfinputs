"""
Crystal2HLT: 共晶结构 → RFantibody HLT 输入自动构建管线

将 RCSB 下载的抗体-抗原共晶结构自动转换为 RFantibody 的 HLT 格式输入。
"""

from .config import (
    __version__,
    Config,
    RoleAssignment,
    AntibodyPairing,
    ResidueTag,
    CDRBoundary,
    HotspotResidue,
    QCReport,
    CHOTHIA_CDR_RANGES,
    CHOTHIA_FV_RANGES,
    PROTEIN_RESIDUES,
)

__all__ = [
    "__version__",
    "Config",
    "RoleAssignment",
    "AntibodyPairing",
    "ResidueTag",
    "CDRBoundary",
    "HotspotResidue",
    "QCReport",
    "CHOTHIA_CDR_RANGES",
    "CHOTHIA_FV_RANGES",
    "PROTEIN_RESIDUES",
]
