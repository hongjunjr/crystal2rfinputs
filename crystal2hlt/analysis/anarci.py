"""
Crystal2HLT ANARCI 封装模块

提供 ANARCI 抗体编号、scFv 检测、CDR 边界识别等功能
"""

from __future__ import annotations

import logging

from ..config import CHOTHIA_CDR_RANGES


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
        
        first_seq_results = numbering[0]
        if not first_seq_results or len(first_seq_results) == 0:
            return None
        
        first_hit = first_seq_results[0]
        numbered_list, start_idx, end_idx = first_hit
        
        chain_type = 'H'
        if alignment_details and len(alignment_details) > 0:
            first_seq_details = alignment_details[0]
            if first_seq_details and len(first_seq_details) > 0:
                chain_type = first_seq_details[0].get('chain_type', 'H')
        
        return {
            'chain_type': chain_type,
            'numbering': numbered_list,
            'scheme': scheme,
            'start': start_idx,
            'end': end_idx
        }
    except ImportError:
        return None
    except Exception:
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
            'vh_first': bool
        }
    或 None 如果不是 scFv
    """
    hits = run_anarci_all_hits(sequence, 'chothia')
    
    if not hits or len(hits) < 2:
        return None
    
    vh_hits = [h for h in hits if h['chain_type'] == 'H']
    vl_hits = [h for h in hits if h['chain_type'] in ['K', 'L']]
    
    if not vh_hits or not vl_hits:
        return None
    
    vh = vh_hits[0]
    vl = vl_hits[0]
    
    logger.info(f"  检测到 scFv: VH({vh['start']}-{vh['end']}) + VL({vl['start']}-{vl['end']})")
    
    return {
        'is_scfv': True,
        'vh_domain': vh,
        'vl_domain': vl,
        'vh_first': vh['start'] < vl['start']
    }


def get_cdr_boundaries(
    anarci_numbering: list,
    chain_type: str,
    offset: int = 0
) -> list[tuple[str, int, int]]:
    """
    从 ANARCI 编号结果提取 CDR 边界
    
    Args:
        anarci_numbering: ANARCI numbering [(pos, aa), ...]
        chain_type: 'H' or 'L' (or 'K')
        offset: 残基索引偏移量（用于 H+L 合并后）
    
    Returns:
        [(cdr_name, abs_start, abs_end), ...]
    """
    if chain_type in ['K', 'L']:
        chain_type = 'L'
    
    cdr_ranges = CHOTHIA_CDR_RANGES.get(chain_type, {})
    if not cdr_ranges:
        return []
    
    boundaries = []
    
    for cdr_name, (chothia_start, chothia_end) in cdr_ranges.items():
        cdr_start = None
        cdr_end = None
        
        for abs_idx, (pos, aa) in enumerate(anarci_numbering):
            if aa == '-':
                continue
            
            chothia_num = pos[0] if isinstance(pos, tuple) else pos
            
            if chothia_num >= chothia_start and cdr_start is None:
                cdr_start = abs_idx + offset
            
            if chothia_num >= chothia_start and chothia_num <= chothia_end:
                cdr_end = abs_idx + offset
        
        if cdr_start is not None and cdr_end is not None:
            boundaries.append((cdr_name, cdr_start, cdr_end))
    
    return boundaries


def identify_antibody_chains_simple(
    sequences: dict[str, str],
    logger: logging.Logger
) -> tuple[list[str], list[str]]:
    """
    简单的抗体链识别（基于序列长度和保守残基）
    用于ANARCI不可用时的后备方案
    """
    heavy_chains = []
    light_chains = []
    
    for chain_id, seq in sequences.items():
        length = len(seq)
        
        if 100 <= length <= 250:
            if 'CXXW' in seq[20:30] if len(seq) > 30 else False:
                heavy_chains.append(chain_id)
            elif 'C' in seq[20:25] if len(seq) > 25 else False:
                light_chains.append(chain_id)
    
    return heavy_chains, light_chains


def parse_pdb_header_for_antibody_chains(
    input_file: str,
    logger: logging.Logger
) -> tuple[list[str], list[str]]:
    """
    从PDB header的COMPND记录解析抗体链信息
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
                    if current_molecule and current_chains:
                        mol_upper = current_molecule.upper()
                        if "HEAVY CHAIN" in mol_upper or "HEAVY" in mol_upper:
                            heavy_chains.extend(current_chains)
                        elif "LIGHT CHAIN" in mol_upper or "LIGHT" in mol_upper:
                            light_chains.extend(current_chains)
                    
                    current_molecule = content.replace("MOLECULE:", "").strip().rstrip(";")
                    current_chains = []
                elif "CHAIN:" in content:
                    chains_str = content.replace("CHAIN:", "").strip().rstrip(";")
                    current_chains = [c.strip() for c in chains_str.split(",")]
            
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
