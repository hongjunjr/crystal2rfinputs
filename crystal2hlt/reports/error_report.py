"""
Crystal2HLT 错误报告模块
"""

from __future__ import annotations

import os
import json
import logging
from datetime import datetime


def generate_error_report(
    config,
    error_type: str,
    error_message: str,
    assignment=None,
    hetatm_info: dict | None = None,
    logger: logging.Logger = None
) -> dict:
    """生成错误报告"""
    
    report = {
        'error_type': error_type,
        'error_message': error_message,
        'input_file': config.input_file,
        'pdbid': config.pdbid,
        'timestamp': datetime.now().isoformat(),
    }
    
    if assignment:
        report['chain_assignment'] = {
            'heavy_chain': assignment.heavy_chain,
            'light_chain': assignment.light_chain,
            'antigen_chains': assignment.antigen_chains
        }
    
    if hetatm_info:
        report['hetatm_info'] = hetatm_info
    
    return report


def write_error_report(report: dict, reports_dir: str, logger: logging.Logger):
    """写入错误报告"""
    os.makedirs(reports_dir, exist_ok=True)
    
    # JSON
    json_path = os.path.join(reports_dir, "error_report.json")
    with open(json_path, 'w') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    logger.info(f"  写入错误报告: {json_path}")
    
    # Markdown
    md_path = os.path.join(reports_dir, "error_report.md")
    with open(md_path, 'w') as f:
        f.write("# Crystal2HLT-RFInputs Error Report\n\n")
        f.write(f"> [!CAUTION]\n")
        f.write(f"> 管线执行失败: {report.get('error_type', 'UNKNOWN')}\n\n")
        
        f.write("## 基本信息\n")
        f.write(f"- **PDB ID**: {report.get('pdbid', '')}\n")
        f.write(f"- **输入文件**: `{report.get('input_file', '')}`\n")
        f.write(f"- **时间戳**: {report.get('timestamp', '')}\n\n")
        
        f.write("## 错误详情\n")
        f.write("```\n")
        f.write(f"{report.get('error_message', '')}\n")
        f.write("```\n\n")
        
        if 'chain_assignment' in report:
            f.write("## 链分配 (执行到此步骤)\n")
            ca = report['chain_assignment']
            f.write(f"- **重链 (H)**: {ca.get('heavy_chain')}\n")
            f.write(f"- **轻链 (L)**: {ca.get('light_chain')}\n")
            f.write(f"- **抗原链 (T)**: {ca.get('antigen_chains')}\n\n")
        
        if 'hetatm_info' in report:
            f.write("## 检测到的非蛋白抗原\n")
            hi = report['hetatm_info']
            f.write(f"- **小分子配体**: {hi.get('ligand', 0)}\n")
            f.write(f"- **糖类残基**: {hi.get('sugar', 0)}\n")
            f.write(f"- **核酸残基**: {hi.get('nucleic', 0)}\n")
            f.write(f"- **金属离子**: {hi.get('metal', 0)}\n")
    
    logger.info(f"  写入错误报告: {md_path}")
