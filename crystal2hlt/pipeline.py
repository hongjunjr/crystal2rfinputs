"""
Crystal2HLT Pipeline 模块

直接使用原始 crystal2hlt_rfinputs.py 的完整功能，
通过导入提供模块化接口。
"""

from __future__ import annotations

import sys
import os
import logging

# 确保可以导入原始模块
_script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

# 从原始文件导入所有需要的函数
from crystal2hlt_rfinputs import (
    # 配置和数据结构
    Config,
    RoleAssignment,
    AntibodyPairing,
    CDRBoundary,
    ResidueTag,
    HotspotResidue,
    QCReport,
    FastaRecord,
    ChainSequenceInfo,
    ChainMapping,
    ChainRoleInfo,
    ResidueMapping,
    PROTEIN_RESIDUES,
    CHOTHIA_CDR_RANGES,
    CHOTHIA_FV_RANGES,
    
    # 日志
    setup_logging,
    
    # FASTA
    parse_fasta,
    extract_chain_ids_from_header,
    extract_role_hints_from_header,
    map_structure_chains_to_fasta,
    get_fasta_hint_for_chain,
    
    # 结构
    detect_file_format,
    normalize_structure,
    filter_protein_only,
    clean_structure,
    extract_chain_sequences,
    extract_chain_sequences_v2,
    
    # Assembly
    extract_mmcif_assemblies,
    build_bioassembly,
    extract_residue_mapping_from_mmcif,
    
    # ANARCI
    run_anarci,
    run_anarci_all_hits,
    detect_scfv,
    identify_antibody_chains_simple,
    parse_pdb_header_for_antibody_chains,
    
    # 角色分配
    assign_roles,
    assign_roles_v2,
    compute_chain_contact_score,
    
    # 导出
    get_chain_residue_mapping,
    identify_cdr_residues,
    export_framework_hlt,
    detect_target_resnum_conflicts,
    export_target_hlt,
    update_residue_numbers,
    
    # 界面
    compute_epitope_patch,
    detect_hetatm_interface,
    
    # Hotspots
    select_hotspots_contact_score,
    spatial_dedup_hotspots,
    select_hotspots_rfantibody_cbeta5,
    select_hotspots,
    
    # 输出
    build_design_loops,
    write_epitope_patch,
    write_hotspots,
    write_design_loops,
    write_cdr_boundaries,
    write_sequences,
    
    # 验证
    validate_hlt_compliance,
    
    # 报告
    generate_qc_report,
    write_qc_report,
    generate_error_report,
    write_error_report,
    
    # 主流程
    setup_directories,
    copy_input_file,
    process_antibody_output,
    run_pipeline,
)

# 版本信息
__version__ = "2.2.0"

# 公开接口
__all__ = [
    # 核心函数
    "run_pipeline",
    "setup_logging",
    
    # 配置
    "Config",
    "RoleAssignment",
    "AntibodyPairing",
    
    # 分析
    "assign_roles_v2",
    "compute_epitope_patch",
    "select_hotspots",
    
    # 导出
    "export_framework_hlt",
    "export_target_hlt",
    
    # 报告
    "generate_qc_report",
    "write_qc_report",
]
