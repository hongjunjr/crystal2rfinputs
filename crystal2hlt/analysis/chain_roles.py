"""
Crystal2HLT 链角色识别模块

提供链角色分配功能，融合 ANARCI + FASTA 提示
"""

from __future__ import annotations

import logging
from dataclasses import field

import numpy as np
from biotite.structure import AtomArray, get_chains

from ..config import (
    Config,
    RoleAssignment,
    AntibodyPairing,
    FastaRecord,
    ChainMapping,
)
from .anarci import run_anarci, detect_scfv
from ..io.structure import extract_chain_sequences


def compute_chain_contact_score(
    structure: AtomArray,
    chain1: str,
    chain2: str,
    cutoff: float = 8.0
) -> float:
    """计算两条链之间的接触分数"""
    from scipy.spatial import KDTree
    
    mask1 = structure.chain_id == chain1
    mask2 = structure.chain_id == chain2
    
    coords1 = structure[mask1].coord
    coords2 = structure[mask2].coord
    
    if len(coords1) == 0 or len(coords2) == 0:
        return 0.0
    
    tree = KDTree(coords2)
    distances, _ = tree.query(coords1, k=1)
    contact_count = np.sum(distances < cutoff)
    
    return float(contact_count)


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
    3. FASTA hint 作为辅助提示
    4. 选择最佳 VH-VL 配对
    """
    logger.info("识别链角色...")
    
    sequences = extract_chain_sequences(structure)
    all_chain_ids = list(sequences.keys())
    logger.info(f"提取到 {len(sequences)} 条链的序列: {all_chain_ids}")
    
    heavy_chains = []
    light_chains = []
    antigen_chains = []
    anarci_results = {}
    scfv_info = {}
    has_scfv = False
    
    # 检查手动指定
    if config.heavy_chain_id or config.light_chain_id:
        logger.info("使用手动指定的链ID...")
        if config.heavy_chain_id:
            heavy_chains = [c.strip() for c in config.heavy_chain_id.split(",")]
        if config.light_chain_id:
            light_chains = [c.strip() for c in config.light_chain_id.split(",")]
        if config.antigen_chain_ids:
            antigen_chains = [c.strip() for c in config.antigen_chain_ids.split(",")]
        else:
            antigen_chains = [c for c in all_chain_ids if c not in heavy_chains + light_chains]
        
        # 获取 ANARCI 编号
        for chain_id in heavy_chains + light_chains:
            if chain_id in sequences:
                result = run_anarci(sequences[chain_id], config.numbering)
                if result:
                    anarci_results[chain_id] = result
    else:
        # ANARCI 自动识别
        logger.info("尝试使用ANARCI进行抗体链识别...")
        
        for chain_id, seq in sequences.items():
            # 检测 scFv
            scfv_result = detect_scfv(seq, logger)
            if scfv_result and scfv_result['is_scfv']:
                logger.info(f"  链 {chain_id}: 检测到 scFv")
                has_scfv = True
                scfv_info[chain_id] = scfv_result
                heavy_chains.append(chain_id)
                continue
            
            result = run_anarci(seq, config.numbering)
            if result:
                chain_type = result['chain_type']
                anarci_results[chain_id] = result
                
                if chain_type == 'H':
                    heavy_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 重链 (VH) [ANARCI]")
                elif chain_type in ['K', 'L']:
                    light_chains.append(chain_id)
                    logger.info(f"  链 {chain_id}: 轻链 (VL) [ANARCI]")
        
        # 非抗体链为抗原
        antibody_chains = set(heavy_chains + light_chains)
        antigen_chains = [c for c in all_chain_ids if c not in antibody_chains]
    
    # 确定模式
    if has_scfv:
        mode = "scfv"
    elif heavy_chains and light_chains:
        mode = "antibody"
    elif heavy_chains and not light_chains:
        mode = "nanobody"
    else:
        mode = "antibody"
    
    logger.info(f"自动检测模式: {mode}")
    
    # 选择最佳配对
    all_pairings = []
    selected_heavy = None
    selected_light = None
    best_score = 0.0
    
    if mode == "scfv":
        for chain_id in heavy_chains:
            if chain_id in scfv_info:
                all_pairings.append(AntibodyPairing(
                    heavy_chain=chain_id,
                    light_chain=None,
                    contact_score=0.0,
                    index=len(all_pairings) + 1,
                    mode="scfv"
                ))
        if all_pairings:
            selected_heavy = all_pairings[0].heavy_chain
    elif mode == "antibody":
        for h in heavy_chains:
            for l in light_chains:
                score = compute_chain_contact_score(structure, h, l)
                all_pairings.append(AntibodyPairing(
                    heavy_chain=h,
                    light_chain=l,
                    contact_score=score,
                    index=len(all_pairings) + 1,
                    mode="antibody"
                ))
                if score > best_score:
                    best_score = score
                    selected_heavy = h
                    selected_light = l
        
        if all_pairings:
            all_pairings.sort(key=lambda x: x.contact_score, reverse=True)
            for i, p in enumerate(all_pairings):
                p.index = i + 1
    else:  # nanobody
        for i, h in enumerate(heavy_chains):
            all_pairings.append(AntibodyPairing(
                heavy_chain=h,
                light_chain=None,
                contact_score=0.0,
                index=i + 1,
                mode="nanobody"
            ))
        if all_pairings:
            selected_heavy = all_pairings[0].heavy_chain
    
    logger.info(f"主抗体: 重链={selected_heavy}, 轻链={selected_light}")
    logger.info(f"抗原链: {antigen_chains}")
    
    return RoleAssignment(
        heavy_chain=selected_heavy,
        light_chain=selected_light,
        antigen_chains=antigen_chains,
        mode=mode,
        vh_vl_contact_score=best_score,
        all_heavy_chains=heavy_chains,
        all_light_chains=light_chains,
        anarci_results=anarci_results,
        all_pairings=all_pairings
    )
