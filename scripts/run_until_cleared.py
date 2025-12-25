"""运行直到Route2清零"""
import sys
from scripts.run_single_weak_th import main as run_single
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--missing_type", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start_weak_th", type=float, default=0.50)
    parser.add_argument("--max_weak_th", type=float, default=0.70)
    parser.add_argument("--step", type=float, default=0.02)
    args = parser.parse_args()
    
    weak_th = args.start_weak_th
    
    while weak_th <= args.max_weak_th:
        print(f"\n{'='*80}")
        print(f"Testing weak_th = {weak_th:.2f}")
        print(f"{'='*80}\n")
        
        # 创建临时参数对象
        class TempArgs:
            def __init__(self):
                self.sample = args.sample
                self.missing_type = args.missing_type
                self.weak_th = weak_th
                self.seed = args.seed
        
        temp_args = TempArgs()
        
        try:
            cleared = run_single(temp_args)
            if cleared:
                print(f"\n{'='*80}")
                print(f"SUCCESS: Route2 cleared at weak_th = {weak_th:.2f}")
                print(f"{'='*80}")
                break
        except Exception as e:
            print(f"Error at weak_th={weak_th:.2f}: {e}")
            weak_th += args.step
            continue
        
        weak_th += args.step
    
    if weak_th > args.max_weak_th:
        print(f"\n{'='*80}")
        print(f"FAILED: Could not clear Route2 with weak_th up to {args.max_weak_th:.2f}")
        print(f"{'='*80}")

if __name__ == "__main__":
    main()

