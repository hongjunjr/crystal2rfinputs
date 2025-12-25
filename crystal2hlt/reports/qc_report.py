"""
Crystal2HLT QC 报告模块
"""

from __future__ import annotations

import os
import json
import logging
from datetime import datetime
from dataclasses import asdict

from ..config import (
    Config,
    RoleAssignment,
    CDRBoundary,
    ResidueTag,
    HotspotResidue,
    QCReport,
    FastaRecord,
    ChainMapping,
)


def generate_qc_report(
    config: Config,
    assignment: RoleAssignment,
    cdr_boundaries: list[CDRBoundary],
    epitope_patch: list[ResidueTag],
    hotspots: list[HotspotResidue],
    hotspots_str: str,
    design_loops_str: str,
    logger: logging.Logger,
    fasta_records: list[FastaRecord] | None = None,
    chain_mappings: dict[str, ChainMapping] | None = None,
    resnum_conflict_info: dict | None = None
) -> QCReport:
    """生成QC报告"""
    
    # 模式描述
    mode_descriptions = {
        'scfv': 'scFv (单链抗体)',
        'antibody': 'Fab/IgG (配对 VH-VL)',
        'nanobody': 'VHH/Nanobody (单域抗体)'
    }
    
    report = QCReport(
        meta={
            'pdbid': config.pdbid,
            'input_mode': config.mode,
            'detected_mode': assignment.mode,
            'antibody_type': mode_descriptions.get(assignment.mode, assignment.mode),
            'numbering': config.numbering,
            'timestamp': datetime.now().isoformat()
        },
        chain_assignment={
            'heavy_chain': assignment.heavy_chain,
            'light_chain': assignment.light_chain,
            'antigen_chains': assignment.antigen_chains,
            'detected_mode': assignment.mode,
            'antibody_type': mode_descriptions.get(assignment.mode, assignment.mode)
        },
        cdr_info={
            'boundaries': [
                {
                    'name': cdr.cdr_name,
                    'chothia_start': cdr.chothia_start,
                    'chothia_end': cdr.chothia_end,
                    'abs_start': cdr.abs_start,
                    'abs_end': cdr.abs_end,
                    'length': cdr.length
                }
                for cdr in cdr_boundaries
            ]
        },
        interface_info={
            'epitope_count': len(epitope_patch),
            'epitope_residues': [tag.to_pdb_format() for tag in epitope_patch]
        },
        hotspots_info={
            'count': len(hotspots),
            'method': config.hotspots_method,
            'hotspots': [h.to_dict() for h in hotspots]
        },
        rfantibody_snippets={
            'hotspot_res': hotspots_str,
            'design_loops': design_loops_str
        },
        success=True
    )
    
    if resnum_conflict_info:
        report.resnum_conflict_info = resnum_conflict_info
    
    return report


def write_qc_report(report: QCReport, reports_dir: str, logger: logging.Logger):
    """写入QC报告"""
    os.makedirs(reports_dir, exist_ok=True)
    
    # JSON
    json_path = os.path.join(reports_dir, "qc_report.json")
    with open(json_path, 'w') as f:
        json.dump(asdict(report), f, indent=2, ensure_ascii=False)
    
    # Markdown
    md_path = os.path.join(reports_dir, "qc_report.md")
    with open(md_path, 'w') as f:
        f.write("# Crystal2HLT-RFInputs QC Report\n\n")
        
        f.write("## 检测到的抗体类型\n")
        f.write(f"- **模式**: `{report.meta.get('detected_mode', 'unknown')}`\n")
        f.write(f"- **类型描述**: {report.meta.get('antibody_type', '')}\n\n")
        
        f.write("## Meta\n")
        for k, v in report.meta.items():
            f.write(f"- **{k.replace('_', ' ').title()}**: {v}\n")
        
        f.write("\n## Chain Assignment\n")
        f.write(f"- **Heavy Chain**: {report.chain_assignment.get('heavy_chain')}\n")
        f.write(f"- **Light Chain**: {report.chain_assignment.get('light_chain')}\n")
        f.write(f"- **Antigen Chains**: {report.chain_assignment.get('antigen_chains')}\n")
        
        f.write("\n## CDR Boundaries\n")
        for cdr in report.cdr_info.get('boundaries', []):
            f.write(f"- **{cdr['name']}**: abs {cdr['abs_start']}-{cdr['abs_end']} (length={cdr['length']})\n")
        
        f.write("\n## Interface\n")
        f.write(f"- **Epitope Count**: {report.interface_info.get('epitope_count', 0)}\n")
        f.write(f"- **Hotspots Count**: {report.hotspots_info.get('count', 0)}\n")
        
        f.write("\n## RFantibody Snippets\n")
        f.write("```\n")
        f.write(f"{report.rfantibody_snippets.get('hotspot_res', '')}\n")
        f.write(f"{report.rfantibody_snippets.get('design_loops', '')}\n")
        f.write("```\n")
    
    logger.info(f"  写入 qc_report.json 和 qc_report.md")
