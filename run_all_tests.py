#!/usr/bin/env python3
"""
Crystal2HLT V2 批量测试脚本
运行所有测试用例并统计结果
"""

import subprocess
import os
import sys
import json
from pathlib import Path

# 配置
SCRIPT_DIR = Path(__file__).parent
SCRIPT_PATH = SCRIPT_DIR / "crystal2hlt_rfinputs.py"
EXAMPLES_DIR = SCRIPT_DIR / "examples"
OUTPUT_DIR = SCRIPT_DIR / "test_outputs"
PYTHON = sys.executable

# 测试用例定义
TEST_CASES = [
    # (PDB_ID, 测试目的, 特点)
    ("1VFB", "基线测试 (Fv + 蛋白抗原)", "D1.3 Fv–溶菌酶"),
    ("1JRH", "标准 Fab + 单链蛋白抗原", "IFN-γ 受体α链片段"),
    ("3HFM", "经典蛋白-蛋白界面", "HyHEL-10 Fab–溶菌酶"),
    ("2IGF", "短肽抗原", "Fab–19aa 合成肽"),
    ("2ZNW", "scFv 识别与拆分", "scFv–溶菌酶"),
    ("5M13", "VHH/合成纳米抗体", "synthetic nanobody–MBP"),
    ("7R74", "VHH + 蛋白抗原", "llama nanobody–HIV gp120"),
    ("6J15", "Fab + 单链受体 + 糖基化", "PD-1–Fab"),
    ("7V2A", "多链抗原三聚体重编号", "SARS-CoV-2 Spike trimer + Fab"),
    ("7WEA", "同一抗原两个Fab", "Omicron Spike + 两个Fab"),
    ("8ZHH", "三聚体 + 两个Fab + 表达标签", "Spike trimer(6P) + H18 Fabs"),
    ("8A1E", "三聚体 + 两种不同Fab", "狂犬病毒糖蛋白 + 两种Fab"),
    ("9JT1", "二聚体抗原 + 两个Fab", "HBsAg二聚体 + Fab"),
    ("9O2R", "多价多拷贝装配", "HIV Env trimer 交联"),
    ("3J8V", "ASU vs biological assembly", "HPV16 衣壳 + Fab"),
    ("6DO1", "膜蛋白 + nanobody", "AT1R (GPCR) + nanobody"),
    ("3P0G", "伴侣蛋白融合 + nanobody", "β2AR + T4L融合"),
    ("8TH4", "GPCR + nanobody + 小分子", "AT1R + nanobody + Losartan"),
    ("1IND", "非蛋白抗原 (金属螯合)", "Fab + indium chelate hapten"),
    ("1S3K", "糖类表位", "Fab + Lewis Y 四糖"),
    ("5E08", "核酸抗原 (RNA)", "Fab + 单链RNA"),
]

def run_test(pdb_id: str) -> dict:
    """运行单个测试用例"""
    cif_file = EXAMPLES_DIR / f"{pdb_id}.cif"
    fasta_file = EXAMPLES_DIR / f"rcsb_pdb_{pdb_id}.fasta"
    output_path = OUTPUT_DIR / f"test_{pdb_id}"
    
    # 检查文件存在
    if not cif_file.exists():
        return {"status": "SKIP", "error": f"CIF file not found: {cif_file}"}
    
    # 构建命令
    cmd = [
        PYTHON, str(SCRIPT_PATH),
        "-i", str(cif_file),
        "-o", str(output_path),
        "--pdbid", pdb_id,
    ]
    
    # 如果有 FASTA 文件则添加
    if fasta_file.exists():
        cmd.extend(["--fasta", str(fasta_file)])
    
    # 运行
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        # 检查输出文件
        framework_hlt = output_path / "outputs" / "framework_HLT.pdb"
        target_hlt = output_path / "outputs" / "target_HLT.pdb"
        hotspots_file = output_path / "outputs" / "epitopes" / "hotspots_rfantibody.txt"
        qc_report = output_path / "reports" / "qc_report.json"
        
        if result.returncode == 0:
            # 检查关键输出文件
            missing_files = []
            if not framework_hlt.exists():
                missing_files.append("framework_HLT.pdb")
            if not target_hlt.exists():
                missing_files.append("target_HLT.pdb")
            if not hotspots_file.exists():
                missing_files.append("hotspots_rfantibody.txt")
            
            if missing_files:
                return {
                    "status": "PARTIAL",
                    "error": f"Missing files: {missing_files}",
                    "stdout": result.stdout[-500:] if len(result.stdout) > 500 else result.stdout
                }
            
            # 读取 QC 报告获取详细信息
            qc_info = {}
            if qc_report.exists():
                with open(qc_report) as f:
                    qc_data = json.load(f)
                    qc_info = {
                        "heavy_chain": qc_data.get("chain_assignment", {}).get("heavy_chain"),
                        "light_chain": qc_data.get("chain_assignment", {}).get("light_chain"),
                        "antigen_chains": qc_data.get("chain_assignment", {}).get("antigen_chains"),
                        "epitope_count": qc_data.get("interface_info", {}).get("epitope_count"),
                        "hotspot_count": qc_data.get("hotspots_info", {}).get("selected_count"),
                        "has_conflict": qc_data.get("resnum_conflict_info", {}).get("has_conflict", False),
                    }
            
            return {"status": "SUCCESS", "qc": qc_info}
        else:
            # 分析错误类型
            error_msg = result.stderr[-500:] if len(result.stderr) > 500 else result.stderr
            if "无法识别抗体链" in error_msg or "无法识别抗体链" in result.stdout:
                return {"status": "FAIL_NO_AB", "error": "无法识别抗体链"}
            elif "界面残基" in error_msg or "无接触" in error_msg:
                return {"status": "FAIL_NO_INTERFACE", "error": "无界面残基"}
            else:
                return {"status": "FAIL", "error": error_msg[:200]}
    
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT", "error": "Timeout after 120s"}
    except Exception as e:
        return {"status": "ERROR", "error": str(e)}

def main():
    # 创建输出目录
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print("=" * 80)
    print("Crystal2HLT V2 批量测试")
    print("=" * 80)
    print()
    
    results = {}
    success_count = 0
    fail_count = 0
    partial_count = 0
    
    for pdb_id, purpose, description in TEST_CASES:
        print(f"[{pdb_id}] {purpose}")
        print(f"         特点: {description}")
        
        result = run_test(pdb_id)
        results[pdb_id] = result
        
        status = result["status"]
        if status == "SUCCESS":
            success_count += 1
            qc = result.get("qc", {})
            print(f"         ✅ 成功 | H={qc.get('heavy_chain')} L={qc.get('light_chain')} "
                  f"Ag={qc.get('antigen_chains')} epitope={qc.get('epitope_count')} "
                  f"hotspot={qc.get('hotspot_count')}")
            if qc.get("has_conflict"):
                print(f"         ⚠️  有 resnum 冲突")
        elif status == "PARTIAL":
            partial_count += 1
            print(f"         ⚠️  部分成功: {result.get('error')}")
        else:
            fail_count += 1
            print(f"         ❌ 失败 [{status}]: {result.get('error', '')[:80]}")
        
        print()
    
    # 汇总
    print("=" * 80)
    print("测试汇总")
    print("=" * 80)
    print(f"总计: {len(TEST_CASES)} 个测试用例")
    print(f"  ✅ 成功: {success_count}")
    print(f"  ⚠️  部分成功: {partial_count}")
    print(f"  ❌ 失败: {fail_count}")
    print()
    
    # 失败详情
    if fail_count > 0:
        print("失败用例详情:")
        for pdb_id, result in results.items():
            if result["status"] not in ["SUCCESS", "PARTIAL"]:
                print(f"  - {pdb_id}: [{result['status']}] {result.get('error', '')[:100]}")
    
    # 保存结果
    results_file = OUTPUT_DIR / "test_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n详细结果已保存到: {results_file}")
    
    return 0 if fail_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
