#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap


HIGHLIGHT_COLOR = "#00CFE8"
HIGHLIGHT_SIZE = 24
HIGHLIGHT_LW = 0.70
BASE_BG = "#E4E0D8"
MARKER_PANEL_TOP_N = 30
RAW_VIS_SPOT_RATIO_TRIGGER = 2.0
RAW_VIS_SPOT_DELTA_TRIGGER = 500


def _make_deep_purple_magma() -> LinearSegmentedColormap:
    # Drop the near-black lower tail so masked panels stay in a consistent deep-purple range.
    base = colormaps["magma"]
    colors = base(np.linspace(0.18, 1.0, 256))
    return LinearSegmentedColormap.from_list("magma_deep_purple", colors)


SIGNATURE_CMAP = _make_deep_purple_magma()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Render standardized real-data single-type profile-mask visualization "
            "(A/B/C/D panels with baseline/route2 overlays)."
        )
    )
    p.add_argument("--sample", required=True, help="profile-mask sample name")
    p.add_argument("--project_root", default=".", help="project root")
    p.add_argument(
        "--force_stage4",
        action="store_true",
        help="re-run stage4 baseline/route2 even if outputs already exist",
    )
    p.add_argument(
        "--mapping_cells_per_spot",
        type=int,
        default=None,
        help="optional manual override for stage4 mapping_cells_per_spot",
    )
    p.add_argument(
        "--n_subspots",
        type=int,
        default=None,
        help="optional manual override for stage4 n_subspots",
    )
    return p.parse_args()


def _slugify(text: str) -> str:
    out = []
    for ch in str(text).strip().lower():
        if ch.isalnum():
            out.append(ch)
        else:
            out.append("_")
    slug = "".join(out)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_") or "target"


def _load_info(raw_dir: Path) -> dict:
    info_path = raw_dir / "real_input_info.json"
    if not info_path.exists():
        raise FileNotFoundError(f"missing info json: {info_path}")
    return json.loads(info_path.read_text(encoding="utf-8"))


def _ensure_windows_junction(link_path: Path, target_path: Path) -> bool:
    if link_path.exists():
        return False
    link_path.parent.mkdir(parents=True, exist_ok=True)
    if os.name == "nt":
        subprocess.run(
            ["cmd", "/c", "mklink", "/J", str(link_path), str(target_path)],
            check=True,
            shell=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    else:
        os.symlink(target_path, link_path, target_is_directory=True)
    return True


def _remove_windows_junction(link_path: Path) -> None:
    if not link_path.exists():
        return
    if os.name == "nt":
        subprocess.run(
            ["cmd", "/c", "rmdir", str(link_path)],
            check=True,
            shell=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    else:
        link_path.unlink()


def _run_quiet_subprocess(cmd: list[str], cwd: Path, label: str) -> None:
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        capture_output=True,
    )
    if proc.returncode == 0:
        print(f"[OK] {label}")
        return

    merged = "\n".join(part for part in [proc.stdout, proc.stderr] if part).strip()
    tail = "\n".join(merged.splitlines()[-60:]) if merged else "(no subprocess output captured)"
    raise RuntimeError(f"{label} failed with exit code {proc.returncode}\n{tail}")


def _stage4_resource_profile(info: dict, mapping_cells_per_spot: int | None, n_subspots: int | None) -> tuple[int, int]:
    if mapping_cells_per_spot is not None and n_subspots is not None:
        return int(mapping_cells_per_spot), int(n_subspots)

    n_spots = int(info.get("n_spots") or 0)
    n_genes = int(info.get("n_genes") or 0)

    # Conservative defaults for large real datasets to avoid CytoSPACE OOM.
    if n_spots >= 7000 or n_genes >= 15000:
        cps, subspots = 1, 300
    elif n_spots >= 4000 or n_genes >= 10000:
        cps, subspots = 2, 400
    else:
        cps, subspots = 5, 800

    if mapping_cells_per_spot is not None:
        cps = int(mapping_cells_per_spot)
    if n_subspots is not None:
        subspots = int(n_subspots)
    return cps, subspots


def _run_stage4(
    project_root: Path,
    sample: str,
    target_type: str,
    force_stage4: bool,
    info: dict,
    mapping_cells_per_spot: int | None,
    n_subspots: int | None,
) -> None:
    result_root = project_root / "result" / sample
    baseline_csv = result_root / "stage4_cytospace_baseline" / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    route2_csv = result_root / "stage4_cytospace_route2" / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    if baseline_csv.exists() and route2_csv.exists() and not force_stage4:
        return

    grouped_processed = project_root / "data" / "processed" / "low_resolution_experiments" / sample
    if not grouped_processed.exists():
        raise FileNotFoundError(f"missing processed profile-mask sample: {grouped_processed}")

    cps, subspots = _stage4_resource_profile(info, mapping_cells_per_spot, n_subspots)
    print(f"[INFO] stage4 resource profile for {sample}: mapping_cells_per_spot={cps}, n_subspots={subspots}")

    top_level_processed = project_root / "data" / "processed" / sample
    created_link = _ensure_windows_junction(top_level_processed, grouped_processed)
    try:
        py = sys.executable
        base_cmd = [
            py,
            "-m",
            "src.stages.stage4_cytospace",
            "--sample",
            sample,
            "--project_root",
            str(project_root),
            "--missing_type",
            target_type,
            "--n_processors",
            "1",
            "--n_subspots",
            str(subspots),
            "--mapping_cells_per_spot",
            str(cps),
            "--sc_expr_source",
            "normalized",
        ]
        print(f"[RUN] stage4 baseline: {sample}")
        _run_quiet_subprocess(
            base_cmd
            + [
                "--filter_mode",
                "none",
                "--cell_type_column",
                "sc_meta",
                "--stage4_suffix",
                "_baseline",
            ],
            project_root,
            f"stage4 baseline for {sample}",
        )
        print(f"[RUN] stage4 route2: {sample}")
        _run_quiet_subprocess(
            base_cmd
            + [
                "--filter_mode",
                "plugin_unknown",
                "--filter_scope",
                "missing_only",
                "--cell_type_column",
                "plugin_type",
                "--stage4_suffix",
                "_route2",
            ],
            project_root,
            f"stage4 route2 for {sample}",
        )
    finally:
        if created_link:
            _remove_windows_junction(top_level_processed)


def _load_marker_panel(raw_dir: Path, orig_cols: list[str], mask_cols: list[str]) -> list[str]:
    panel_path = raw_dir / "st_profile_mask_gene_panel.csv"
    if not panel_path.exists():
        raise FileNotFoundError(f"missing marker panel: {panel_path}")
    panel = pd.read_csv(panel_path)
    if "gene" not in panel.columns:
        raise ValueError(f"marker panel missing gene column: {panel_path}")
    genes: list[str] = []
    allowed = set(orig_cols).intersection(mask_cols)
    for gene in panel["gene"].astype(str):
        g = gene.strip()
        if g and g in allowed and g not in genes:
            genes.append(g)
        if len(genes) >= MARKER_PANEL_TOP_N:
            break
    if not genes:
        raise ValueError(f"no shared panel genes found for {panel_path}")
    return genes


def _raw_st_paths(raw_sample_dir: Path) -> tuple[Path, Path]:
    exp_path = raw_sample_dir / "brca_STdata_GEP.txt"
    coord_path = raw_sample_dir / "brca_STdata_coordinates.txt"
    if not exp_path.exists():
        raise FileNotFoundError(exp_path)
    if not coord_path.exists():
        raise FileNotFoundError(coord_path)
    return exp_path, coord_path


def _raw_st_genes(raw_sample_dir: Path) -> list[str]:
    exp_path, _ = _raw_st_paths(raw_sample_dir)
    gene_col = pd.read_csv(exp_path, sep="\t", usecols=[0])
    return gene_col.iloc[:, 0].astype(str).tolist()


def _raw_st_spot_count(raw_sample_dir: Path) -> int:
    exp_path, coord_path = _raw_st_paths(raw_sample_dir)
    expr_spots = max(0, len(pd.read_csv(exp_path, sep="\t", nrows=0).columns) - 1)
    coord_spots = sum(1 for _ in coord_path.open("r", encoding="utf-8")) - 1
    return min(expr_spots, coord_spots)


def _should_use_raw_visual_st(raw_orig_dir: Path, raw_mask_dir: Path, processed_spots: int) -> bool:
    try:
        raw_spots = min(_raw_st_spot_count(raw_orig_dir), _raw_st_spot_count(raw_mask_dir))
    except FileNotFoundError:
        return False
    return (
        raw_spots >= processed_spots * RAW_VIS_SPOT_RATIO_TRIGGER
        and raw_spots >= processed_spots + RAW_VIS_SPOT_DELTA_TRIGGER
    )


def _load_raw_st_df(raw_sample_dir: Path, genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    exp_path, coord_path = _raw_st_paths(raw_sample_dir)
    gene_set = set(genes)
    wanted: list[pd.DataFrame] = []
    lib_size: pd.Series | None = None

    for chunk in pd.read_csv(exp_path, sep="\t", chunksize=512):
        gene_col = chunk.columns[0]
        expr = chunk.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sums = expr.sum(axis=0)
        lib_size = sums if lib_size is None else lib_size.add(sums, fill_value=0.0)

        hit = chunk[gene_col].astype(str).isin(gene_set)
        if hit.any():
            picked = chunk.loc[hit].copy()
            picked[gene_col] = picked[gene_col].astype(str)
            wanted.append(picked)

    if not wanted or lib_size is None:
        raise ValueError(f"no requested genes found in raw ST matrix: {exp_path}")

    mat = pd.concat(wanted, ignore_index=True).drop_duplicates(subset=[wanted[0].columns[0]], keep="first")
    mat = mat.set_index(mat.columns[0])
    use_genes = [g for g in genes if g in mat.index]
    mat = mat.loc[use_genes].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    lib_size = lib_size.reindex(mat.columns).replace(0, np.nan)
    norm = np.log1p(mat.divide(lib_size, axis=1) * 10000.0).fillna(0.0)
    expr = norm.T.reset_index().rename(columns={"index": "spot_id"})

    coords = pd.read_csv(coord_path, sep="\t", usecols=["spot_id", "row", "col"])
    coords["spot_id"] = coords["spot_id"].astype(str)
    df = coords.merge(expr, on="spot_id", how="inner")
    return df, use_genes


def _load_st_df(processed_sample_dir: Path, genes: list[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    export_dir = processed_sample_dir / "stage1_preprocess" / "exported"
    exp_path = export_dir / "st_expression_normalized.csv"
    coord_path = export_dir / "st_coordinates.csv"
    if not exp_path.exists():
        raise FileNotFoundError(exp_path)
    if not coord_path.exists():
        raise FileNotFoundError(coord_path)

    cols = pd.read_csv(exp_path, nrows=0).columns.tolist()
    if genes is None:
        return pd.read_csv(coord_path, usecols=["spot_id", "row", "col"]), cols

    use_genes = [g for g in genes if g in cols]
    expr = pd.read_csv(exp_path, usecols=["spot_id", *use_genes])
    coords = pd.read_csv(coord_path, usecols=["spot_id", "row", "col"])
    df = coords.merge(expr, on="spot_id", how="inner")
    return df, use_genes


def _resolve_target_column(df: pd.DataFrame, target_type: str) -> str:
    if target_type in df.columns:
        return target_type
    norm = {str(c).strip().casefold(): c for c in df.columns}
    hit = norm.get(str(target_type).strip().casefold())
    if hit is not None:
        return str(hit)
    raise KeyError(f"target type column not found in assignments: {target_type}")


def _load_assignment_overlay(project_root: Path, sample: str, target_type: str, tag: str) -> pd.DataFrame:
    assign_path = project_root / "result" / sample / f"stage4_cytospace_{tag}" / "cytospace_output" / "cell_type_assignments_by_spot.csv"
    if not assign_path.exists():
        raise FileNotFoundError(assign_path)
    df = pd.read_csv(assign_path).rename(columns=lambda c: str(c))
    df = df.rename(columns={df.columns[0]: "spot_id"})
    target_col = _resolve_target_column(df, target_type)
    df[target_col] = pd.to_numeric(df[target_col], errors="coerce").fillna(0)
    if "Total cells" in df.columns:
        total = pd.to_numeric(df["Total cells"], errors="coerce").fillna(0)
    else:
        ignore = {"spot_id", "Unknown_sc_only"}
        total = df[[c for c in df.columns if c not in ignore]].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    df["target_fraction"] = (df[target_col] / total.replace(0, np.nan)).fillna(0)
    df["target_dominant"] = df["target_fraction"] >= 0.5
    return df[["spot_id", target_col, "target_dominant"]].rename(columns={target_col: "target_cells"})


def _merge_overlay(mask: pd.DataFrame, overlay: pd.DataFrame) -> pd.DataFrame:
    merged = mask.merge(overlay, on="spot_id", how="left")
    merged["target_cells"] = pd.to_numeric(merged["target_cells"], errors="coerce").fillna(0)
    merged["target_dominant"] = [
        False if pd.isna(value) else bool(value)
        for value in merged["target_dominant"].tolist()
    ]
    return merged


def _draw_panel(ax, df: pd.DataFrame, title: str, subtitle: str, mask_col: str, show_ylabel: bool = False):
    ax.scatter(df["col"], -df["row"], s=6, c=BASE_BG, linewidths=0, alpha=0.72, rasterized=True)
    sc = ax.scatter(
        df["col"],
        -df["row"],
        s=10,
        c=df["target_signature"],
        cmap=SIGNATURE_CMAP,
        norm=_draw_panel.norm,
        linewidths=0,
        alpha=0.95,
        rasterized=True,
    )
    hi = df[df[mask_col].astype(bool)] if mask_col in df.columns else df.iloc[0:0]
    if not hi.empty:
        ax.scatter(
            hi["col"],
            -hi["row"],
            s=HIGHLIGHT_SIZE,
            facecolors="none",
            edgecolors=HIGHLIGHT_COLOR,
            linewidths=HIGHLIGHT_LW,
            alpha=1.0,
            antialiaseds=True,
            rasterized=False,
            zorder=10,
        )
    ax.set_title(f"{title}\n{subtitle}", loc="left", fontweight="bold", fontsize=10.2, pad=6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("spatial col", fontsize=8.3)
    ax.set_ylabel("spatial row" if show_ylabel else "", fontsize=12.0 if show_ylabel else 8.3)
    return sc


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    raw_dir = project_root / "data" / "raw" / "low_resolution_experiments" / args.sample
    processed_mask = project_root / "data" / "processed" / "low_resolution_experiments" / args.sample
    info = _load_info(raw_dir)
    source_sample = str(info.get("source_sample") or "").strip()
    target_type = str(((info.get("st_profile_mask") or {}).get("target_type")) or "").strip()
    if not source_sample or not target_type:
        raise ValueError(f"cannot infer source_sample/target_type from {raw_dir / 'real_input_info.json'}")

    _run_stage4(
        project_root,
        args.sample,
        target_type,
        args.force_stage4,
        info,
        args.mapping_cells_per_spot,
        args.n_subspots,
    )

    processed_orig = project_root / "data" / "processed" / "low_resolution_experiments" / source_sample
    raw_orig = project_root / "data" / "raw" / "low_resolution_experiments" / source_sample

    processed_orig_coords, orig_cols = _load_st_df(processed_orig, None)
    _, mask_cols = _load_st_df(processed_mask, None)
    use_raw_visual_st = _should_use_raw_visual_st(raw_orig, raw_dir, len(processed_orig_coords))
    if use_raw_visual_st:
        print(
            f"[INFO] using raw ST matrix for visualization: "
            f"processed_spots={len(processed_orig_coords)}, raw_spots={_raw_st_spot_count(raw_orig)}"
        )
        orig_cols = _raw_st_genes(raw_orig)
        mask_cols = _raw_st_genes(raw_dir)

    genes = _load_marker_panel(raw_dir, orig_cols, mask_cols)

    if use_raw_visual_st:
        orig, use_genes = _load_raw_st_df(raw_orig, genes)
        mask, _ = _load_raw_st_df(raw_dir, use_genes)
    else:
        orig, use_genes = _load_st_df(processed_orig, genes)
        mask, _ = _load_st_df(processed_mask, use_genes)

    orig["target_signature"] = orig[use_genes].mean(axis=1)
    mask["target_signature"] = mask[use_genes].mean(axis=1)
    cutoff = float(orig["target_signature"].quantile(0.90))
    orig["expr_high"] = orig["target_signature"] >= cutoff
    mask["expr_high"] = mask["target_signature"] >= cutoff
    combined = pd.concat([orig["target_signature"], mask["target_signature"]], ignore_index=True)
    _draw_panel.norm = Normalize(float(np.nanpercentile(combined, 1)), float(np.nanpercentile(combined, 99)))

    mask_base = _merge_overlay(mask, _load_assignment_overlay(project_root, args.sample, target_type, "baseline"))
    mask_route = _merge_overlay(mask, _load_assignment_overlay(project_root, args.sample, target_type, "route2"))

    target_short = target_type
    slug = _slugify(target_type)
    out_dir = project_root / "visualizations" / "masked_scenarios" / "real" / args.sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{slug}_mask_ABCD_expression_identical_cyan_outline.png"

    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(1, 4, figsize=(15.8, 5.0), dpi=220)
    ax_a, ax_b, ax_c, ax_d = axes

    sc = _draw_panel(ax_a, orig, "A  Original real ST", f"{target_short} marker signature score", "expr_high", show_ylabel=True)
    _draw_panel(ax_b, mask, "B  Profile-mask real ST", f"{target_short} signature after masking", "expr_high")
    _draw_panel(ax_c, mask_base, "C  Baseline mapping overlay", f"baseline {target_short} outline", "target_dominant")
    _draw_panel(ax_d, mask_route, "D  Route2 mapping overlay", f"route2 {target_short} outline", "target_dominant")

    fig.suptitle(
        f"Real-data single-type profile masking: {target_short} ST-expression signal with identical cyan outline overlays",
        x=0.03,
        ha="left",
        fontsize=14.0,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.024, right=0.928, top=0.80, bottom=0.16, wspace=-0.07)
    cax = fig.add_axes([0.941, 0.22, 0.013, 0.56])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(f"{target_short} marker signature\nmean normalized expression", fontsize=7.3, labelpad=5)
    cb.ax.tick_params(labelsize=7.3)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"[OK] wrote: {out_png}")
    print(f"[INFO] sample={args.sample} source_sample={source_sample} target_type={target_type}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
