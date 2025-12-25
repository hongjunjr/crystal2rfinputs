"""
Crystal2HLT 摘要报告模块
"""

from __future__ import annotations

import os
import json
import logging
from datetime import datetime


def generate_summary_report(
    config,
    all_antibody_infos: list[dict],
    logger: logging.Logger
) -> dict:
    """生成全局摘要报告"""
    
    return {
        'input_file': config.input_file,
        'pdbid': config.pdbid,
        'timestamp': datetime.now().isoformat(),
        'antibody_count': len(all_antibody_infos),
        'antibodies': all_antibody_infos
    }


def write_summary_report(
    summary: dict,
    outdir: str,
    logger: logging.Logger
):
    """写入摘要报告"""
    
    # JSON
    json_path = os.path.join(outdir, "summary_report.json")
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    logger.info(f"  写入 {json_path}")
    
    # Markdown
    md_path = os.path.join(outdir, "summary_report.md")
    with open(md_path, 'w') as f:
        f.write("# Crystal2HLT 处理摘要\n\n")
        f.write(f"- **输入文件**: {summary.get('input_file', 'N/A')}\n")
        f.write(f"- **PDB ID**: {summary.get('pdbid', 'N/A')}\n")
        f.write(f"- **处理时间**: {summary.get('timestamp', '')}\n")
        f.write(f"- **检测到抗体数量**: {summary.get('antibody_count', 0)}\n")
        
        # 检测对称拷贝
        symmetric_copies = [
            ab for ab in summary.get('antibodies', [])
            if ab.get('is_symmetric_copy', False)
        ]
        if symmetric_copies:
            f.write(f"\n## ⚠ 对称相关拷贝\n\n")
            f.write(f"检测到 {len(symmetric_copies)} 个对称相关拷贝: ")
            f.write(str([ab.get('index', 0) for ab in symmetric_copies]))
            f.write("\n\n这些抗体与其他抗体具有相同的序列，可能是晶体堆积产生的对称相关分子。\n")
        
        f.write("\n## 抗体列表\n\n")
        f.write("| # | 模式 | 重链 | 轻链 | Epitope | Hotspots | 对称拷贝 |\n")
        f.write("|---|------|------|------|---------|----------|----------|\n")
        
        for ab in summary.get('antibodies', []):
            sym = "✓" if ab.get('is_symmetric_copy') else ""
            f.write(f"| {ab.get('index', '')} | {ab.get('mode', '')} | "
                   f"{ab.get('heavy_chain', '')} | {ab.get('light_chain', '-')} | "
                   f"{ab.get('epitope_count', 0)} | {ab.get('hotspots_count', 0)} | {sym} |\n")
        
        f.write("\n## 输出目录\n\n")
        for ab in summary.get('antibodies', []):
            f.write(f"- `antibody_{ab.get('index', '')}/` - {ab.get('mode', '')}\n")
    
    logger.info(f"  写入 {md_path}")
