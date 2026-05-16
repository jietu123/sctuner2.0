from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OFFICIAL_DIR = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "cytospace_fig2c_melanoma"
    / "official_example"
    / "CytoSPACE_example_melanoma"
)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run CytoSPACE Fig.2c melanoma using official example inputs.")
    p.add_argument("--project_root", default=str(PROJECT_ROOT))
    p.add_argument("--official_dir", default=str(OFFICIAL_DIR))
    p.add_argument("--out_dir", default=None)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--solver_method", default="lapjv", choices=["lapjv", "lapjv_compat", "lap_CSPR"])
    p.add_argument("--downsample_off", action="store_true")
    p.add_argument("--n_processors", type=int, default=1)
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    project_root = Path(args.project_root).resolve()
    official_dir = Path(args.official_dir).resolve()
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else project_root / "result" / "cytospace_fig2c_official_slide1" / "cytospace_output"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "scRNA_path": official_dir / "melanoma_scRNA_GEP.txt",
        "cell_type_path": official_dir / "melanoma_scRNA_celllabels.txt",
        "st_path": official_dir / "melanoma_STdata_slide1_GEP.txt",
        "coordinates_path": official_dir / "melanoma_STdata_slide1_coordinates.txt",
        "cell_type_fraction_estimation_path": official_dir / "melanoma_cell_fraction_estimates.txt",
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)

    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "external" / "cytospace"))
    from src.stages.stage4_cytospace import _patch_datatable_fread

    _patch_datatable_fread()
    from cytospace.cytospace import main_cytospace  # type: ignore

    _patch_datatable_fread()
    main_cytospace(
        scRNA_path=str(paths["scRNA_path"]),
        cell_type_path=str(paths["cell_type_path"]),
        n_cells_per_spot_path=None,
        st_cell_type_path=None,
        cell_type_fraction_estimation_path=str(paths["cell_type_fraction_estimation_path"]),
        st_path=str(paths["st_path"]),
        coordinates_path=str(paths["coordinates_path"]),
        output_folder=str(out_dir),
        output_prefix="",
        mean_cell_numbers=20,
        downsample_off=bool(args.downsample_off),
        solver_method=str(args.solver_method),
        distance_metric="Pearson_correlation",
        sampling_method="duplicates",
        single_cell=False,
        sampling_sub_spots=False,
        number_of_processors=int(args.n_processors),
        seed=int(args.seed),
        plot_off=True,
        geometry="square",
    )
    print(f"[OK] official CytoSPACE Fig.2c run wrote: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
