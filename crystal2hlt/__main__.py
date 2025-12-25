"""
Crystal2HLT 命令行入口

运行: python -m crystal2hlt -i input.cif -o output_dir
"""

from .cli import main

if __name__ == "__main__":
    main()
