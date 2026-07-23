#!/usr/bin/env python
"""Prepare the human lymph node cell2location input used by NCEM Figure 2."""

from __future__ import annotations

import argparse
import gc
import json
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata
import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from anndata import AnnData
from cell2location.models import Cell2location, RegressionModel
from cell2location.utils.filtering import filter_genes
from scipy import sparse


CELL_TYPE_MAPPING = {
    "B_Cycling": "B cells",
    "B_GC_DZ": "B cells",
    "B_GC_LZ": "B cells",
    "B_GC_prePB": "B cells",
    "B_IFN": "B cells",
    "B_activated": "B cells",
    "B_mem": "B cells",
    "B_naive": "B cells",
    "B_plasma": "B cells",
    "B_preGC": "B cells",
    "DC_CCR7+": "DC",
    "DC_cDC1": "DC",
    "DC_cDC2": "DC",
    "DC_pDC": "DC",
    "Endo": "Endo",
    "FDC": "FDC",
    "ILC": "ILC",
    "Macrophages_M1": "Macrophages",
    "Macrophages_M2": "Macrophages",
    "Mast": "Mast",
    "Monocytes": "Monocytes",
    "NK": "NK",
    "NKT": "NKT",
    "T_CD4+": "CD4 T cells",
    "T_CD4+_TfH": "CD4 T cells",
    "T_CD4+_TfH_GC": "CD4 T cells",
    "T_CD4+_naive": "CD4 T cells",
    "T_CD8+_CD161+": "CD8 T cells",
    "T_CD8+_cytotoxic": "CD8 T cells",
    "T_CD8+_naive": "CD8 T cells",
    "T_TIM3+": "T TIM",
    "T_TfR": "T TfR",
    "T_Treg": "T Treg",
    "VSMC": "VSMC",
}


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    raw_root = repo_root / "data" / "raw" / "ncem_fig2_human_lymph_node"
    output_root = repo_root / "data" / "processed" / "ncem_fig2_human_lymph_node"
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=raw_root)
    parser.add_argument("--output-root", type=Path, default=output_root)
    parser.add_argument("--reference-epochs", type=int, default=250)
    parser.add_argument("--spatial-epochs", type=int, default=30000)
    parser.add_argument("--reference-posterior-samples", type=int, default=1000)
    parser.add_argument("--spatial-posterior-samples", type=int, default=1000)
    parser.add_argument(
        "--spatial-posterior-mode",
        choices=("samples", "quantile"),
        default="samples",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--stop-after",
        choices=("reference", "spatial", "ncem"),
        default="ncem",
    )
    parser.add_argument("--force-reference", action="store_true")
    parser.add_argument("--force-spatial", action="store_true")
    return parser.parse_args()


def require_inputs(raw_root: Path) -> tuple[Path, Path]:
    sc_path = raw_root / "scrna_reference" / "sc.h5ad"
    visium_root = raw_root / "visium"
    required = [
        sc_path,
        visium_root / "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
        visium_root / "spatial" / "tissue_positions_list.csv",
        visium_root / "spatial" / "scalefactors_json.json",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    return sc_path, visium_root


def save_model_unicode_safe(model, destination: Path) -> None:
    """Save through an ASCII temporary path for PyTorch on Windows."""
    destination.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(tempfile.mkdtemp(prefix="ncem_model_"))
    temporary_model = temporary_root / "model"
    temporary_model.mkdir()
    try:
        model.save(temporary_model, overwrite=True, save_anndata=False)
        shutil.copytree(temporary_model, destination, dirs_exist_ok=True)
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def load_model_unicode_safe(model_class, source: Path, adata: AnnData):
    """Load through an ASCII temporary path for PyTorch on Windows."""
    temporary_root = Path(tempfile.mkdtemp(prefix="ncem_model_"))
    temporary_model = temporary_root / "model"
    try:
        shutil.copytree(source, temporary_model)
        return model_class.load(
            str(temporary_model),
            adata=adata,
            use_gpu=True,
        )
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def load_visium(visium_root: Path) -> AnnData:
    adata_vis = sc.read_visium(
        visium_root,
        count_file="V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
    )
    sample = next(iter(adata_vis.uns["spatial"]))
    adata_vis.obs["sample"] = sample
    adata_vis.var["SYMBOL"] = adata_vis.var_names.astype(str)
    adata_vis.var_names = adata_vis.var["gene_ids"].astype(str)
    adata_vis.var_names.name = None
    adata_vis.var_names_make_unique()

    mt_mask = adata_vis.var["SYMBOL"].astype(str).str.startswith("MT-").to_numpy()
    adata_vis.obsm["MT"] = adata_vis[:, mt_mask].X.toarray()
    return adata_vis[:, ~mt_mask].copy()


def load_reference(sc_path: Path) -> AnnData:
    adata_ref = sc.read_h5ad(sc_path)
    adata_ref.var["SYMBOL"] = adata_ref.var_names.astype(str)
    adata_ref.var_names = adata_ref.var["GeneID-2"].astype(str)
    adata_ref.var_names.name = None

    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_ref = adata_ref[:, selected].copy()
    adata_ref.obs["ncem_cluster"] = (
        adata_ref.obs["Subset"].map(CELL_TYPE_MAPPING).astype("category")
    )
    if adata_ref.obs["ncem_cluster"].isna().any():
        missing = sorted(
            adata_ref.obs.loc[
                adata_ref.obs["ncem_cluster"].isna(), "Subset"
            ].astype(str).unique()
        )
        raise ValueError(f"Unmapped reference cell types: {missing}")
    return adata_ref


def train_reference(
    adata_ref: AnnData,
    output_root: Path,
    epochs: int,
    posterior_samples: int,
) -> pd.DataFrame:
    signatures_path = output_root / "reference_signatures.csv.gz"
    if signatures_path.exists():
        print(f"Using cached reference signatures: {signatures_path}", flush=True)
        return pd.read_csv(signatures_path, index_col=0)

    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key="Sample",
        labels_key="ncem_cluster",
        categorical_covariate_keys=["Method"],
    )
    model = RegressionModel(adata_ref)
    print(f"Training reference regression model for {epochs} epochs", flush=True)
    model.train(
        max_epochs=epochs,
        batch_size=2500,
        train_size=1,
        lr=0.002,
        use_gpu=True,
    )
    reference_model_path = output_root / "reference_model"
    save_model_unicode_safe(model, reference_model_path)
    adata_ref = model.export_posterior(
        adata_ref,
        sample_kwargs={
            "num_samples": posterior_samples,
            "batch_size": 2500,
            "use_gpu": True,
        },
    )
    key = "means_per_cluster_mu_fg"
    if key in adata_ref.varm:
        factor_names = adata_ref.uns["mod"]["factor_names"]
        columns = [f"{key}_{name}" for name in factor_names]
        signatures = adata_ref.varm[key][columns].copy()
    else:
        factor_names = adata_ref.uns["mod"]["factor_names"]
        columns = [f"{key}_{name}" for name in factor_names]
        signatures = adata_ref.var[columns].copy()
    signatures.columns = factor_names
    signatures.to_csv(signatures_path, compression="gzip")
    pd.DataFrame(model.history_["elbo_train"]).to_csv(
        output_root / "reference_training_history.csv"
    )
    del model, adata_ref
    gc.collect()
    torch.cuda.empty_cache()
    return signatures


def train_spatial(
    adata_vis: AnnData,
    signatures: pd.DataFrame,
    output_root: Path,
    epochs: int,
    posterior_samples: int,
    posterior_mode: str,
    force: bool,
) -> tuple[AnnData, Cell2location]:
    common = np.intersect1d(adata_vis.var_names, signatures.index)
    if common.size == 0:
        raise ValueError("Reference and Visium data have no shared ENSEMBL gene IDs")
    adata_vis = adata_vis[:, common].copy()
    signatures = signatures.loc[common].copy()
    print(
        f"Spatial model input: {adata_vis.n_obs} spots, "
        f"{adata_vis.n_vars} genes, {signatures.shape[1]} cell types",
        flush=True,
    )

    Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    spatial_model_path = output_root / "spatial_model"
    if spatial_model_path.joinpath("model.pt").exists() and not force:
        print(f"Using cached spatial model: {spatial_model_path}", flush=True)
        model = load_model_unicode_safe(Cell2location, spatial_model_path, adata_vis)
    else:
        model = Cell2location(
            adata_vis,
            cell_state_df=signatures,
            N_cells_per_location=30,
            detection_alpha=200,
        )
        print(f"Training spatial cell2location model for {epochs} epochs", flush=True)
        model.train(
            max_epochs=epochs,
            batch_size=None,
            train_size=1,
            use_gpu=True,
        )
        save_model_unicode_safe(model, spatial_model_path)
        pd.DataFrame(model.history_["elbo_train"]).to_csv(
            output_root / "spatial_training_history.csv"
        )

    if posterior_mode == "samples":
        print(
            f"Sampling spatial posterior with {posterior_samples} samples",
            flush=True,
        )
        adata_vis = model.export_posterior(
            adata_vis,
            sample_kwargs={
                "num_samples": posterior_samples,
                "batch_size": model.adata.n_obs,
                "use_gpu": True,
            },
        )
    else:
        print("Computing spatial posterior q05 directly", flush=True)
        adata_vis = model.export_posterior(
            adata_vis,
            sample_kwargs={"batch_size": 512, "use_gpu": True},
            add_to_obsm=["q05"],
            use_quantiles=True,
        )
    factor_names = list(adata_vis.uns["mod"]["factor_names"])
    adata_vis.obs[factor_names] = adata_vis.obsm["q05_cell_abundance_w_sf"]
    adata_vis.write_h5ad(
        output_root / "cell2location_spatial_q05.h5ad",
        compression="gzip",
    )
    return adata_vis, model


def build_ncem_input(
    adata_vis: AnnData,
    model: Cell2location,
    output_root: Path,
) -> Path:
    expected = model.module.model.compute_expected_per_cell_type(
        adata_vis.uns["mod"]["post_sample_q05"],
        model.adata_manager,
    )
    factor_names = list(model.factor_names_)
    target_order = sorted(set(CELL_TYPE_MAPPING.values()))
    if set(factor_names) != set(target_order):
        raise ValueError(
            "Spatial model factors differ from official NCEM cell types: "
            f"{factor_names}"
        )

    reorder = [factor_names.index(cell_type) for cell_type in target_order]
    mu = expected["mu"]
    n_types = len(target_order)
    n_spots = adata_vis.n_obs
    cell_expression = sparse.vstack(
        [mu[index].astype(np.float32) for index in reorder],
        format="csr",
    )

    proportions = np.tile(
        adata_vis.obs[target_order].to_numpy(dtype=np.float32),
        (n_types, 1),
    )
    node_types = np.repeat(
        np.eye(n_types, dtype=np.float32),
        n_spots,
        axis=0,
    )
    spatial = np.tile(
        adata_vis.obsm["spatial"].astype(np.float32),
        (n_types, 1),
    )

    adata_ncem = AnnData(cell_expression)
    adata_ncem.obsm["proportions"] = proportions
    adata_ncem.obsm["node_types"] = node_types
    adata_ncem.obsm["spatial"] = spatial
    adata_ncem.uns["node_type_names"] = {
        cell_type: cell_type for cell_type in target_order
    }
    adata_ncem.var_names = adata_vis.var["SYMBOL"].astype(str).to_numpy()
    adata_ncem.var_names_make_unique()

    sc.pp.log1p(adata_ncem)
    sc.pp.highly_variable_genes(adata_ncem, n_top_genes=2000, subset=True)
    adata_ncem.obs["target_cell"] = np.repeat(target_order, n_spots)

    output_path = output_root / "cell2location_lymphnode.h5ad"
    adata_ncem.write_h5ad(output_path, compression="gzip")
    print(
        f"NCEM input written: {output_path} "
        f"({adata_ncem.n_obs} pseudo-cells x {adata_ncem.n_vars} HVGs)",
        flush=True,
    )
    return output_path


def write_metadata(args: argparse.Namespace, output_root: Path) -> None:
    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_device": torch.cuda.get_device_name(0)
        if torch.cuda.is_available()
        else None,
        "cell2location": cell2location.__version__,
        "scvi_tools": scvi.__version__,
        "scanpy": sc.__version__,
        "anndata": anndata.__version__,
        "reference_epochs": args.reference_epochs,
        "spatial_epochs": args.spatial_epochs,
        "reference_posterior_samples": args.reference_posterior_samples,
        "spatial_posterior_samples": args.spatial_posterior_samples,
        "spatial_posterior_mode": args.spatial_posterior_mode,
        "seed": args.seed,
    }
    (output_root / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA GPU is required for the official training schedule")
    torch.set_float32_matmul_precision("high")
    scvi.settings.seed = args.seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    args.output_root.mkdir(parents=True, exist_ok=True)
    write_metadata(args, args.output_root)
    sc_path, visium_root = require_inputs(args.raw_root)

    signatures_path = args.output_root / "reference_signatures.csv.gz"
    if args.force_reference and signatures_path.exists():
        signatures_path.unlink()
    adata_ref = load_reference(sc_path)
    signatures = train_reference(
        adata_ref,
        args.output_root,
        args.reference_epochs,
        args.reference_posterior_samples,
    )
    if args.stop_after == "reference":
        return

    adata_vis = load_visium(visium_root)
    spatial_path = args.output_root / "cell2location_spatial_q05.h5ad"
    if args.force_spatial and spatial_path.exists():
        spatial_path.unlink()
    adata_vis, model = train_spatial(
        adata_vis,
        signatures,
        args.output_root,
        args.spatial_epochs,
        args.spatial_posterior_samples,
        args.spatial_posterior_mode,
        args.force_spatial,
    )
    if args.stop_after == "spatial":
        return
    build_ncem_input(adata_vis, model, args.output_root)


if __name__ == "__main__":
    main()
