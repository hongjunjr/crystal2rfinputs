# Crystal2HLT-RFInputs v2.2

将 RCSB 下载的抗体-抗原共晶结构自动转换为 RFantibody 的 HLT 格式输入。

## 🆕 v2.2 新功能

| 功能 | 说明 |
|------|------|
| **模块化包结构** | `python -m crystal2hlt` 运行，可作为 Python 库导入 |
| **FASTA 辅助识别** | 使用 PDB FASTA 自动识别链角色（推荐） |
| **多 Fab 自动配对** | 自动选择最佳 VH-VL 配对，支持多抗体结构 |
| **scFv 支持** | 自动检测并正确拆分 scFv 单链抗体的 VH/VL 域 |
| **对称拷贝检测** | 识别晶体对称相关分子并标记 |
| **自动重编号** | 多链 target 合并后自动消除 resnum 冲突 |
| **非蛋白抗原检测** | 识别半抗原/糖/核酸等非蛋白配体并给出提示 |

## 安装

```bash
# 创建 conda 环境
conda env create -f environment.yml

# 激活环境
conda activate crystal2hlt

# 验证安装
python -c "import biotite, anarci, Bio; print('✓ 安装成功')"
```

## 快速开始

### 方式1: FASTA 辅助识别 (推荐)

使用原始脚本 + PDB FASTA 文件，可以最准确地识别链角色：

```bash
# 下载结构和 FASTA
curl -o 1VFB.cif https://files.rcsb.org/download/1VFB.cif
curl -o 1VFB.fasta https://www.rcsb.org/fasta/entry/1VFB

# 运行 pipeline (推荐方式)
python crystal2hlt_rfinputs.py \
    -i 1VFB.cif \
    -o output_1vfb \
    --fasta 1VFB.fasta \
    --pdbid 1VFB \
    -v
```

### 方式2: 自动识别 (简单快速)

无需 FASTA，直接使用 ANARCI 自动识别抗体链：

```bash
python crystal2hlt_rfinputs.py -i examples/1VFB.cif -o output_1vfb
```

### 方式3: 模块化运行 (可选)

也可以使用模块化包运行（功能与原始脚本完全相同）：

```bash
python -m crystal2hlt -i examples/1VFB.cif -o output_1vfb --fasta examples/1VFB.fasta -v
```

### 方式4: 手动指定链 ID

对于自动识别失败的结构，可以手动指定：

```bash
python crystal2hlt_rfinputs.py \
    -i 1VFB.cif \
    -o output_1vfb \
    --heavy_chain B \
    --light_chain A \
    --pdbid 1VFB \
    -v
```

## 参数说明

### 必需参数

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入结构文件 (.cif 或 .pdb) |
| `-o, --outdir` | 输出目录 |

### V2 FASTA 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--fasta, -f` | - | PDB FASTA 文件，用于辅助链角色识别 |
| `--fasta_mode` | `auto` | FASTA 使用模式: `auto`/`strict`/`off` |

### 链指定参数

| 参数 | 说明 |
|------|------|
| `--heavy_chain, -H` | 重链 ID (逗号分隔多个) |
| `--light_chain, -L` | 轻链 ID (逗号分隔多个) |
| `--antigen_chains, -A` | 抗原链 ID (留空则自动推断) |
| `--target_chains, -T` | 手动指定 target 链 (覆盖自动选择) |

### 其他参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--pdbid` | - | PDB ID |
| `--mode` | `auto` | 模式: `antibody`/`nanobody`/`scfv`/`auto` |
| `--hotspots_method` | `contact_score` | 热点方法: `contact_score`/`rfantibody_cbeta5` |
| `--hotspots_k` | `12` | 选择的热点数量 |
| `-v, --verbose` | - | 详细输出 |

## 输出结构

```
output/
├── antibody_1/                         # 第一套抗体完整输出
│   ├── outputs/
│   │   ├── framework_HLT.pdb           # 抗体 Fv (H+L链, 残基从1编号, 含CDR REMARK)
│   │   ├── target_HLT.pdb              # 该抗体结合的抗原 (T链, 自动重编号)
│   │   ├── target_resmap.tsv           # 残基映射表 (多链时生成)
│   │   ├── epitopes/
│   │   │   ├── epitope_patch.txt       # 界面残基列表
│   │   │   ├── hotspots.txt            # 热点残基
│   │   │   └── hotspots_rfantibody.txt # RFantibody 格式
│   │   ├── loops/
│   │   │   ├── cdr_boundaries.json     # CDR 边界
│   │   │   └── design_loops_rfantibody.txt
│   │   └── sequences/
│   │       ├── H.fasta, L.fasta, T.fasta
│   └── reports/
│       ├── qc_report.json              # QC 报告
│       └── qc_report.md
├── antibody_2/                         # 第二套抗体 (如有)
│   └── ...
├── input/                              # 原始输入副本
├── work/                               # 中间处理文件
├── logs/                               # 日志 (--verbose 时生成)
├── summary_report.json                 # 全局摘要
└── summary_report.md                   # 全局摘要 (人类可读)
```

### QC 报告包含

- **检测到的抗体类型**: scfv / antibody / nanobody
- **链分配**: 重链、轻链、抗原链
- **CDR 边界**: Chothia 编号
- **对称拷贝检测**: 标记晶体对称相关分子

## 用于 RFantibody

```bash
# 复制 hotspots 和 design_loops 到 RFantibody 命令
poetry run python /path/to/rfantibody/scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=output/antibody_1/outputs/target_HLT.pdb \
    antibody.framework_pdb=output/antibody_1/outputs/framework_HLT.pdb \
    inference.ckpt_override_path=/path/to/RFdiffusion_Ab.pt \
    'ppi.hotspot_res=[T125,T121,T116,T23,T27,T103,T129,T18]' \
    'antibody.design_loops=[H1:7,H2:5,H3:8,L1:11,L2:7,L3:9]' \
    inference.num_designs=20 \
    inference.output_prefix=output/design
```

## 项目结构

```
crystal2rfinputs/
├── crystal2hlt_rfinputs.py    # 完整实现 (3854行)
├── crystal2hlt/               # 模块化包
│   ├── __init__.py            # 公共接口导出
│   ├── __main__.py            # python -m crystal2hlt 入口
│   ├── cli.py                 # 命令行解析
│   ├── pipeline.py            # 所有函数导入
│   ├── config.py              # 数据结构定义
│   ├── io/                    # 结构加载/FASTA/导出
│   ├── analysis/              # ANARCI/链角色/界面
│   ├── validation/            # HLT 验证
│   └── reports/               # 报告生成
├── examples/                  # 测试结构
├── environment.yml
└── requirements.txt
```

## Python API

```python
# 导入使用
from crystal2hlt.pipeline import run_pipeline, Config

config = Config(
    input_file="1VFB.cif",
    outdir="output_1vfb",
    fasta_file="1VFB.fasta",
    pdbid="1VFB",
    verbose=True
)
run_pipeline(config)
```

## 测试结果 (v2.2)

| 场景 | 测试用例 | 状态 |
|------|----------|------|
| 基线 Fab | 1VFB | ✅ |
| VHH/Nanobody | 5M13, 7R74, 6DO1, 3P0G | ✅ |
| scFv | 2ZNW | ✅ |
| 多链抗原 | 6J15, 7V2A | ✅ |
| 多 Fab | 7WEA, 8ZHH, 8A1E, 9JT1, 3J8V | ✅ |
| GPCR | 8TH4 | ✅ |
| 对称拷贝检测 | 2ZNW, 7V2A, 7R74 | ✅ |
| 非蛋白抗原 | 1S3K (糖), 5E08 (RNA) | ⚠️ 预期失败 |

**蛋白抗原通过率: 14/14 (100%)**

## 依赖

- Python ≥3.9
- biotite ≥0.37.0 (结构解析)
- biopython ≥1.80 (FASTA 解析)
- numpy, scipy (计算)
- anarci (抗体编号, 通过 conda 安装)
- hmmer (ANARCI 依赖, 通过 conda 安装)

## 更新日志

### v2.2 (2024-12-25)

- 模块化包结构 (`crystal2hlt/` 包)
- 支持 `python -m crystal2hlt` 运行方式，与原始脚本 100% 功能一致
- 提供 Python API：`from crystal2hlt.pipeline import run_pipeline`
- 代码拆分为独立模块: io, analysis, validation, reports
- 推荐 FASTA 辅助识别作为默认使用方式

### v2.1 (2024-12)

- 抗体中心化输出结构 (`antibody_N/` 目录)
- 全局摘要报告 (`summary_report.json/md`)
- 对称拷贝检测与标记
- QC 报告显示检测到的模式 (scfv/antibody/nanobody)
- scFv VH/VL 域正确拆分为 H/L 链
- 每个抗体独立的抗原链选择

### v2.0 (2024-12)

- 新增 FASTA 辅助链角色识别
- 新增 scFv 自动检测和处理
- 新增多链 target 自动重编号 (消除冲突)
- 新增非蛋白抗原检测和优雅失败提示
- 增强 QC 报告 (fasta_info, chain_mapping_info, resnum_conflict_info)

### v1.0

- ANARCI Chothia 编号
- HLT 合规性验证
- Biological Assembly 支持
- rfantibody_cbeta5 热点方法

## License

MIT
