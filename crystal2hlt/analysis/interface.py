"""
Crystal2HLT 界面计算模块

提供 epitope 和 hotspot 计算功能
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.spatial import KDTree
from biotite.structure import AtomArray, residue_iter

from ..config import (
    Config,
    RoleAssignment,
    ResidueTag,
    HotspotResidue,
    CDRBoundary,
)


def compute_chain_contact_score(
    structure: AtomArray,
    chain1: str,
    chain2: str,
    cutoff: float = 8.0
) -> float:
    """计算两条链之间的接触分数"""
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


def compute_epitope_patch(
    structure: AtomArray,
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger
) -> list[ResidueTag]:
    """
    计算抗原界面残基（epitope）
    """
    logger.info(f"计算界面残基 (cutoff={config.interface_cutoff}Å)...")
    
    ab_chains = [assignment.heavy_chain]
    if assignment.light_chain:
        ab_chains.append(assignment.light_chain)
    
    ab_mask = np.isin(structure.chain_id, ab_chains)
    ab_coords = structure[ab_mask].coord
    
    if len(ab_coords) == 0:
        logger.warning("未找到抗体原子")
        return []
    
    ab_tree = KDTree(ab_coords)
    epitope_residues = set()
    
    for ag_chain in assignment.antigen_chains:
        chain_mask = structure.chain_id == ag_chain
        chain_atoms = structure[chain_mask]
        
        for residue in residue_iter(chain_atoms):
            for atom in residue:
                distances, _ = ab_tree.query(atom.coord.reshape(1, -1), k=1)
                if distances[0] <= config.interface_cutoff:
                    icode = atom.ins_code if hasattr(atom, 'ins_code') else ''
                    tag = ResidueTag(
                        chain_id='T',
                        res_num=atom.res_id,
                        icode=icode,
                        res_name=atom.res_name
                    )
                    epitope_residues.add(tag)
                    break
    
    epitope_list = sorted(epitope_residues, key=lambda x: (x.chain_id, x.res_num))
    logger.info(f"  找到 {len(epitope_list)} 个界面残基")
    
    return epitope_list


def select_hotspots(
    structure: AtomArray,
    epitope_patch: list[ResidueTag],
    assignment: RoleAssignment,
    config: Config,
    logger: logging.Logger,
    cdr_boundaries: list[CDRBoundary] = None
) -> list[HotspotResidue]:
    """
    选择热点残基
    使用 contact_score 方法
    """
    logger.info(f"使用 {config.hotspots_method} 方法选择热点...")
    
    if not epitope_patch:
        return []
    
    # 获取抗体原子
    ab_chains = [assignment.heavy_chain]
    if assignment.light_chain:
        ab_chains.append(assignment.light_chain)
    
    ab_mask = np.isin(structure.chain_id, ab_chains)
    ab_coords = structure[ab_mask].coord
    
    if len(ab_coords) == 0:
        return []
    
    ab_tree = KDTree(ab_coords)
    
    scored_hotspots = []
    
    for tag in epitope_patch:
        chain_mask = np.isin(structure.chain_id, assignment.antigen_chains)
        res_mask = structure.res_id == tag.res_num
        residue_mask = chain_mask & res_mask
        residue_atoms = structure[residue_mask]
        
        if len(residue_atoms) == 0:
            continue
        
        distances, _ = ab_tree.query(residue_atoms.coord, k=1)
        contact_count = np.sum(distances < config.interface_cutoff)
        min_dist = float(np.min(distances))
        
        score = contact_count + (1.0 / max(min_dist, 0.1))
        
        scored_hotspots.append(HotspotResidue(
            tag=tag,
            score=score,
            contact_count=int(contact_count),
            min_distance=min_dist
        ))
    
    # 排序并选择 top-k
    scored_hotspots.sort(key=lambda x: x.score, reverse=True)
    selected = scored_hotspots[:config.hotspots_k]
    
    logger.info(f"  选择了 {len(selected)} 个热点")
    
    return selected
