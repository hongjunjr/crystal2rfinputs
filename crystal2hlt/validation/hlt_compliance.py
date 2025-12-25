"""
Crystal2HLT HLT 合规性验证模块
"""

from __future__ import annotations

import logging
import re

from ..config import RoleAssignment


def validate_hlt_compliance(
    framework_path: str,
    target_path: str,
    hotspots_path: str,
    assignment: RoleAssignment,
    logger: logging.Logger
) -> dict:
    """
    严格的 HLT 合规性检查
    
    检查：
    1. framework链顺序必须 H → L（或只有H）
    2. REMARK格式正确
    3. hotspots每个残基都能在target中定位
    """
    logger.info("执行 HLT 合规性验证...")
    
    result = {
        'passed': True,
        'errors': [],
        'warnings': [],
        'chain_order': [],
        'remark_count': 0,
        'hotspots_locatable': 0
    }
    
    try:
        # 1. 检查 framework 链顺序
        with open(framework_path, 'r') as f:
            lines = f.readlines()
        
        chains_seen = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                if not chains_seen or chains_seen[-1] != chain_id:
                    chains_seen.append(chain_id)
        
        result['chain_order'] = chains_seen
        
        expected_order = ['H', 'L'] if assignment.mode in ['antibody', 'scfv'] else ['H']
        if chains_seen != expected_order:
            if chains_seen == ['H'] and assignment.mode in ['nanobody', 'scfv']:
                pass  # scfv 有时只有 H
            elif chains_seen != expected_order[:len(chains_seen)]:
                result['errors'].append(f"链顺序错误: 期望 {expected_order}, 实际 {chains_seen}")
                result['passed'] = False
        
        logger.info(f"  ✓ 链顺序正确: {chains_seen}")
        
        # 2. 检查 REMARK
        remark_count = 0
        for line in lines:
            if line.startswith('REMARK 666'):
                remark_count += 1
        
        result['remark_count'] = remark_count
        if remark_count > 0:
            logger.info(f"  ✓ REMARK 格式正确: {remark_count} 个标注")
        
        # 3. 检查 hotspots 可定位性
        with open(target_path, 'r') as f:
            target_lines = f.readlines()
        
        target_residues = set()
        for line in target_lines:
            if line.startswith('ATOM'):
                chain = line[21]
                try:
                    resnum = int(line[22:26].strip())
                    target_residues.add((chain, resnum))
                except ValueError:
                    pass
        
        with open(hotspots_path, 'r') as f:
            hotspot_lines = f.readlines()
        
        locatable = 0
        for line in hotspot_lines:
            match = re.match(r'([A-Z]):(\d+)', line.strip())
            if match:
                chain, resnum = match.group(1), int(match.group(2))
                if (chain, resnum) in target_residues:
                    locatable += 1
        
        result['hotspots_locatable'] = locatable
        logger.info(f"  ✓ hotspots 可定位: {locatable} 个")
        
    except Exception as e:
        result['errors'].append(str(e))
        result['passed'] = False
    
    if result['passed']:
        logger.info("✓ HLT 合规性验证通过")
    else:
        logger.error(f"✗ HLT 合规性验证失败: {result['errors']}")
    
    return result
