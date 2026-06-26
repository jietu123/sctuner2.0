"""Command-line entrypoint for SVTuner."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from svtuner import __version__
from svtuner.bundle import create_bundle
from svtuner.run import run_pipeline


def _project_root(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).resolve()
    return Path(__file__).resolve().parents[2]


def cmd_run(args: argparse.Namespace) -> int:
    root = _project_root(args.project_root)
    os.chdir(root)
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    try:
        run_pipeline(
            root,
            args.sample,
            from_scratch=args.from_scratch,
            skip_simgen=args.skip_simgen,
            skip_stage1=args.skip_stage1,
            use_python_stage1=args.use_python_stage1,
            missing_type=args.missing_type,
            n_processors=args.n_processors,
            n_subspots=args.n_subspots,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    return 0


def cmd_envcheck(args: argparse.Namespace) -> int:
    root = _project_root(args.project_root)
    ret = subprocess.run([sys.executable, "-m", "src.stages.stage0_envcheck"], cwd=root)
    return int(ret.returncode)


def cmd_stage3b(args: argparse.Namespace) -> int:
    root = _project_root(args.project_root)
    cmd = [
        sys.executable,
        "-m",
        "src.stages.stage3b_st_unsupported",
        "--project_root",
        str(root),
        "--sample",
        args.sample,
        "--fdr",
        str(args.fdr),
        "--n_calibration",
        str(args.n_calibration),
        "--n_spatial_permutations",
        str(args.n_spatial_permutations),
        "--random_seed",
        str(args.random_seed),
        "--sc_expr_source",
        args.sc_expr_source,
        "--sc_profile_source",
        args.sc_profile_source,
        "--expression_scale",
        args.expression_scale,
        "--sc_profile_scale",
        args.sc_profile_scale,
        "--max_genes",
        str(args.max_genes),
    ]
    ret = subprocess.run(cmd, cwd=root)
    return int(ret.returncode)


def cmd_bundle(args: argparse.Namespace) -> int:
    root = _project_root(args.project_root)
    out_dir = (root / args.out_dir).resolve() if not Path(args.out_dir).is_absolute() else Path(args.out_dir).resolve()
    bundle_name = args.name or f"svtuner_bundle_{__version__}"
    try:
        info = create_bundle(
            root,
            out_dir,
            bundle_name,
            include_raw_data=args.include_raw_data,
            include_results=args.include_results,
            include_external=not args.no_external,
        )
    except Exception as e:
        print(f"[ERROR] bundle failed: {e}", file=sys.stderr)
        return 1

    print(json.dumps(info, ensure_ascii=False, indent=2))
    return 0


def cmd_version(_: argparse.Namespace) -> int:
    print(__version__)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="svtuner", description="SVTuner unified CLI")
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("run", help="run full pipeline")
    pr.add_argument("--project-root", default=None, help="repo root path")
    pr.add_argument("--sample", required=True, help="dataset id (configs/datasets/<sample>.yaml)")
    pr.add_argument("--from-scratch", action="store_true", help="run simgen preset before pipeline")
    pr.add_argument("--skip-simgen", action="store_true", help="skip simgen when using --from-scratch")
    pr.add_argument("--skip-stage1", action="store_true", help="skip Stage1 if already prepared")
    pr.add_argument("--use-python-stage1", action="store_true", help="use Python Stage1 preparation instead of R")
    pr.add_argument("--missing-type", default=None, help="override missing_type")
    pr.add_argument("--n-processors", type=int, default=None, help="Stage4 n_processors")
    pr.add_argument("--n-subspots", type=int, default=None, help="Stage4 n_subspots")
    pr.set_defaults(_fn=cmd_run)

    pe = sub.add_parser("envcheck", help="run Stage0 environment check")
    pe.add_argument("--project-root", default=None, help="repo root path")
    pe.set_defaults(_fn=cmd_envcheck)

    ps = sub.add_parser(
        "stage3b",
        help="detect ST regions unsupported by the SC reference",
    )
    ps.add_argument("--project-root", default=None, help="repo root path")
    ps.add_argument("--sample", required=True, help="dataset id")
    ps.add_argument("--fdr", type=float, default=0.05, help="BH FDR level")
    ps.add_argument(
        "--n-calibration",
        type=int,
        default=0,
        help="pseudo-ST count; 0 uses observed spot count",
    )
    ps.add_argument(
        "--n-spatial-permutations",
        type=int,
        default=200,
        help="spatial region permutation count",
    )
    ps.add_argument("--random-seed", type=int, default=42)
    ps.add_argument(
        "--sc-expr-source",
        choices=["normalized", "data", "counts", "auto"],
        default="normalized",
    )
    ps.add_argument(
        "--sc-profile-source",
        choices=["normalized", "data", "counts", "auto"],
        default="normalized",
    )
    ps.add_argument(
        "--expression-scale",
        choices=["log1p", "linear"],
        default="log1p",
        help="scale of Stage1 expression values",
    )
    ps.add_argument(
        "--sc-profile-scale",
        choices=["log1p", "linear"],
        default="log1p",
    )
    ps.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="compute-budget cap; 0 uses all common genes",
    )
    ps.set_defaults(_fn=cmd_stage3b)

    pb = sub.add_parser("bundle", help="create distributable project bundle zip")
    pb.add_argument("--project-root", default=None, help="repo root path")
    pb.add_argument("--out-dir", default="dist", help="output directory for bundle")
    pb.add_argument("--name", default=None, help="bundle base name (without .zip)")
    pb.add_argument("--include-raw-data", action="store_true", help="include data/raw in bundle")
    pb.add_argument("--include-results", action="store_true", help="include result/ in bundle")
    pb.add_argument("--no-external", action="store_true", help="exclude external/cytospace code")
    pb.set_defaults(_fn=cmd_bundle)

    pv = sub.add_parser("version", help="print version")
    pv.set_defaults(_fn=cmd_version)
    return p


def main() -> int:
    parser = build_parser()
    ns = parser.parse_args()
    fn = getattr(ns, "_fn", None)
    if fn is None:
        parser.print_help()
        return 2
    return int(fn(ns))


if __name__ == "__main__":
    raise SystemExit(main())

