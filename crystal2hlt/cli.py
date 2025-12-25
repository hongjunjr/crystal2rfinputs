"""
Crystal2HLT 命令行接口

使用原始 crystal2hlt_rfinputs.py 的完整功能
"""

from __future__ import annotations

import argparse
import sys
import os

# 确保可以导入原始模块
_script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

from crystal2hlt_rfinputs import Config, run_pipeline, setup_logging


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="Crystal2HLT: 共晶结构 → RFantibody HLT 输入",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python -m crystal2hlt -i examples/1VFB.cif -o output_1vfb
  python -m crystal2hlt -i structure.pdb -o output --fasta structure.fasta -v
        """
    )
    
    # 必需参数
    parser.add_argument('-i', '--input', required=True,
                       help='输入结构文件 (.cif 或 .pdb)')
    parser.add_argument('-o', '--outdir', required=True,
                       help='输出目录')
    
    # FASTA 参数
    parser.add_argument('-f', '--fasta', default='',
                       help='PDB FASTA 文件（辅助链识别）')
    parser.add_argument('--fasta_mode', default='auto',
                       choices=['auto', 'strict', 'off'],
                       help='FASTA 使用模式')
    
    # 链指定参数
    parser.add_argument('-H', '--heavy_chain', default='',
                       help='手动指定重链 ID')
    parser.add_argument('-L', '--light_chain', default='',
                       help='手动指定轻链 ID')
    parser.add_argument('-A', '--antigen_chains', default='',
                       help='手动指定抗原链 ID（逗号分隔）')
    parser.add_argument('-T', '--target_chains', default='',
                       help='手动指定 target 链 ID（逗号分隔）')
    
    # 模式和编号
    parser.add_argument('--pdbid', default='',
                       help='PDB ID')
    parser.add_argument('--mode', default='auto',
                       choices=['antibody', 'nanobody', 'scfv', 'auto'],
                       help='抗体模式')
    parser.add_argument('--numbering', default='chothia',
                       choices=['chothia', 'imgt', 'aho'],
                       help='编号方案')
    
    # Hotspots
    parser.add_argument('--hotspots_method', default='contact_score',
                       choices=['contact_score', 'rfantibody_cbeta5'],
                       help='热点选择方法')
    parser.add_argument('--hotspots_k', type=int, default=12,
                       help='选择的热点数量')
    parser.add_argument('--interface_cutoff', type=float, default=4.5,
                       help='界面距离阈值 (Å)')
    
    # Assembly
    parser.add_argument('--use_bioassembly', default='auto',
                       choices=['auto', 'true', 'false'],
                       help='是否使用 biological assembly')
    parser.add_argument('--bioassembly_id', type=int, default=1,
                       help='Biological assembly ID')
    
    # 多抗体
    parser.add_argument('--multi_antibody', default='top1',
                       choices=['off', 'top1', 'topN'],
                       help='多抗体处理模式')
    
    # 其他
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='详细输出')
    
    return parser.parse_args()


def main():
    """主入口函数"""
    args = parse_args()
    
    # 构建配置
    config = Config(
        input_file=args.input,
        outdir=args.outdir,
        pdbid=args.pdbid,
        mode=args.mode,
        numbering=args.numbering,
        interface_cutoff=args.interface_cutoff,
        hotspots_k=args.hotspots_k,
        hotspots_method=args.hotspots_method,
        use_bioassembly=args.use_bioassembly,
        bioassembly_id=args.bioassembly_id,
        heavy_chain_id=args.heavy_chain,
        light_chain_id=args.light_chain,
        antigen_chain_ids=args.antigen_chains,
        target_chain_ids=args.target_chains,
        fasta_file=args.fasta,
        fasta_mode=args.fasta_mode,
        multi_antibody=args.multi_antibody,
        verbose=args.verbose
    )
    
    # 运行 pipeline
    try:
        run_pipeline(config)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
