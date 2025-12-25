"""
Crystal2HLT FASTA 解析与链映射模块

处理 PDB FASTA 文件的解析和结构链到 FASTA 记录的映射
"""

from __future__ import annotations

import re
import logging
from difflib import SequenceMatcher

from ..config import (
    FastaRecord,
    ChainSequenceInfo,
    ChainMapping,
)


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
    """
    chain_ids = []
    
    match = re.search(r'Chains?\s+([A-Za-z0-9,\s\[\]auth]+?)(?:\||$)', header)
    if match:
        chains_str = match.group(1)
        auth_matches = re.findall(r'\[auth\s+([A-Za-z0-9]+)\]', chains_str)
        if auth_matches:
            chain_ids.extend(auth_matches)
        else:
            simple_chains = re.findall(r'([A-Za-z0-9])(?:\s*,|\s*$|(?=\[))', chains_str)
            if simple_chains:
                chain_ids.extend(simple_chains)
            else:
                single_match = re.search(r'Chain\s+([A-Za-z0-9])', header)
                if single_match:
                    chain_ids.append(single_match.group(1))
    
    return list(set(chain_ids))


def extract_role_hints_from_header(header: str) -> list[str]:
    """
    提取角色关键词: heavy, light, vhh, vnar, scfv, fab, nanobody
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
    """
    parts = header.split('|')
    hints = []
    
    if len(parts) >= 3:
        description = parts[2].strip()
        if description:
            hints.append(description)
    
    return hints


def compute_fasta_parse_confidence(record: FastaRecord) -> float:
    """
    计算 FASTA 解析置信度
    """
    score = 0.0
    
    if record.chain_ids:
        score += 0.4
    if record.role_hints:
        score += 0.3
    if record.entity_hints:
        score += 0.2
    if 50 < len(record.sequence) < 1000:
        score += 0.1
    
    return min(score, 1.0)


def compute_sequence_alignment(seq1: str, seq2: str) -> tuple[float, float, float]:
    """
    计算两条序列的比对分数
    返回: (identity, coverage_seq1, coverage_seq2)
    """
    if not seq1 or not seq2:
        return 0.0, 0.0, 0.0
    
    matcher = SequenceMatcher(None, seq1, seq2)
    identity = matcher.ratio()
    
    matching_blocks = matcher.get_matching_blocks()
    total_match = sum(block.size for block in matching_blocks)
    
    cov1 = total_match / len(seq1) if seq1 else 0
    cov2 = total_match / len(seq2) if seq2 else 0
    
    return identity, cov1, cov2


def map_chains_to_fasta(
    chain_seqs: dict[str, ChainSequenceInfo],
    fasta_records: list[FastaRecord],
    mode: str,
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
        
        candidates = []
        for fasta_rec in fasta_records:
            priority = 0.0
            if chain_id in fasta_rec.chain_ids:
                priority = 1.0
            elif abs(len(fasta_rec.sequence) - chain_info.length) < 30:
                priority = 0.5
            
            if priority > 0:
                candidates.append((fasta_rec, priority))
        
        best_match = None
        best_score = 0.0
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
