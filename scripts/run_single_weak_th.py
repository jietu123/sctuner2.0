"""运行单个weak_th值的完整流程并输出结果"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

def detect_root() -> Path:
    here = Path(__file__).resolve()
    try:
        return here.parents[1]
    except IndexError:
        return here.parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--missing_type", required=True)
    parser.add_argument("--weak_th", type=float, required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    
    root = detect_root()
    
    # 更新配置
    cfg_path = root / "configs" / "datasets" / f"{args.sample}.yaml"
    import yaml
    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    if "stage3" not in cfg:
        cfg["stage3"] = {}
    cfg["stage3"]["weak_th"] = args.weak_th
    cfg["stage3"]["strong_th"] = 0.80
    cfg["stage3"]["abstain_unknown_sc_only"] = True
    with cfg_path.open("w", encoding="utf-8") as f:
        yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True)
    
    # 运行Stage3
    print(f"Running Stage3 with weak_th={args.weak_th}...")
    cmd = ["E:\\ANACONDA\\envs\\cytospace_v1.1.0_py310\\python.exe", 
           "src/stages/stage3_type_plugin.py",
           "--sample", args.sample,
           "--dataset_config", str(cfg_path)]
    subprocess.run(cmd, cwd=root, check=True)
    
    # 运行Stage4 Route2
    print(f"Running Stage4 Route2...")
    cmd = ["E:\\ANACONDA\\envs\\cytospace_v1.1.0_py310\\python.exe",
           "src/stages/stage4_cytospace.py",
           "--sample", args.sample,
           "--filter_mode", "plugin_unknown",
           "--cell_type_column", "plugin_type",
           "--missing_type", args.missing_type,
           "--mean_cell_numbers", "5",
           "--solver_method", "lap_CSPR",
           "--seed", str(args.seed),
           "--sampling_sub_spots",
           "--n_subspots", "400",
           "--n_processors", "1"]
    subprocess.run(cmd, cwd=root, check=True)
    
    # 运行Stage5 Route2
    print(f"Running Stage5 Route2...")
    mt_slug = args.missing_type.lower().replace(" ", "_")
    sim_dir = root / "data" / "sim" / "real_brca" / "S0" / mt_slug / f"seed_{args.seed}"
    out_dir = root / "result" / args.sample / "stage5_route2_s0"
    stage4_dir = root / "result" / args.sample / "stage4_cytospace" / "cytospace_output"
    
    cmd = ["E:\\ANACONDA\\envs\\cytospace_v1.1.0_py310\\python.exe",
           "src/stages/stage5_route2_s0.py",
           "--sample", args.sample,
           "--run_tag", "route2",
           "--stage4_dir", str(stage4_dir),
           "--sim_dir", str(sim_dir),
           "--out_dir", str(out_dir)]
    subprocess.run(cmd, cwd=root, check=True)
    
    # 读取结果
    json_path = out_dir / "stage5_route2_s0__route2.json"
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    
    # 输出关键指标
    print(f"\n{'='*60}")
    print(f"Results for weak_th={args.weak_th}")
    print(f"{'='*60}")
    leakage = data.get("leakage", {})
    audit = data.get("filter_audit", {})
    
    print(f"missing_present_in_assignment: {leakage.get('missing_in_assignment', -1)}")
    if audit:
        print(f"missing_type_contrib_count: {audit.get('missing_type_contribution', {}).get('count', 0)}")
        print(f"total_filtered: {audit.get('total_filtered', 0)}")
        print(f"non_missing_filtered_fraction: {audit.get('non_missing_filtered', {}).get('fraction_of_total', 0.0)*100:.2f}%")
        print(f"top3_filtered_types:")
        for item in audit.get('top10_filtered_types', [])[:3]:
            print(f"  - {item['type']}: {item['count']} ({item['fraction']*100:.2f}%)")
    
    return leakage.get('missing_in_assignment', -1) == 0

if __name__ == "__main__":
    main()

