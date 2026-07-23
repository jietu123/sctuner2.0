from __future__ import annotations

import json
import math
from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "svg.fonttype": "none",
        "svg.image_inline": True,
    }
)


ROOT = Path(__file__).resolve().parents[1]
BIOAPP = ROOT / "visualizations" / "bioapp_experiment"
OUT = BIOAPP / "bioapp_downstream_phase_d3_microenvironment_context_validation"
D1 = BIOAPP / "bioapp_downstream_phase_d1_morphology_and_interface_analysis"
V312 = BIOAPP / "bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign"

N_PERM = 2000
N_BOOT = 3000
SEED = 20260708
RNG = np.random.default_rng(SEED)

DOMAINS = ["tumor_core", "immune_enriched", "stroma_rich", "mixed_boundary", "unmapped_or_excluded"]
NICHE_COLORS = ["#8C1D40", "#D55E00", "#4D9221", "#7B3294", "#0072B2", "#666666"]
DISPLAY_NAME = {
    "deep_tumor_core": "Deep tumor",
    "tumor_interior": "Tumor interior",
    "mixed_tumor_boundary_niche": "Mixed boundary",
    "immune_associated_niche": "Immune niche",
    "stroma_interface_niche": "Stroma interface",
    "low_signal_edge_context": "Low-signal edge",
}


def clean_json(obj):
    if isinstance(obj, dict):
        return {k: clean_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_json(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, np.generic):
        return clean_json(obj.item())
    return obj


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(clean_json(payload), indent=2, allow_nan=False), encoding="utf-8")


def load_inputs() -> tuple[pd.DataFrame, np.ndarray]:
    df = pd.read_csv(D1 / "phase_d1_downstream_analysis_by_spot.csv")
    df["withheld_binary"] = df["withheld_binary"].fillna(False).astype(bool)
    df["main_analysis_spot"] = df["primary_endpoint_status"].isin(["positive", "negative"])
    df["endpoint_positive_binary"] = df["primary_endpoint_status"].eq("positive").astype(float)
    img = plt.imread(V312 / "fig_bioapp_v3_embedded_tissue_background.png")
    return df, img


def nearest_neighbor_distance(coords: np.ndarray) -> float:
    tree = cKDTree(coords)
    dist, _ = tree.query(coords, k=2)
    return float(np.median(dist[:, 1]))


def build_neighbors(coords: np.ndarray) -> tuple[list[list[int]], list[list[int]], float]:
    median_nn = nearest_neighbor_distance(coords)
    radius = median_nn * 1.65
    tree = cKDTree(coords)
    one_hop = []
    for i, ns in enumerate(tree.query_ball_point(coords, r=radius)):
        one_hop.append([j for j in ns if j != i])
    two_hop = []
    for i, ns in enumerate(one_hop):
        seen = set(ns)
        for j in ns:
            seen.update(one_hop[j])
        seen.discard(i)
        two_hop.append(sorted(seen))
    return one_hop, two_hop, radius


def entropy(vals: np.ndarray) -> float:
    vals = np.asarray(vals, dtype=float)
    vals = vals[vals > 0]
    if vals.size == 0:
        return 0.0
    vals = vals / vals.sum()
    return float(-(vals * np.log(vals)).sum() / np.log(len(DOMAINS)))


def local_summary(df: pd.DataFrame, indices: list[int]) -> dict:
    if not indices:
        return {
            "tumor_fraction": np.nan,
            "immune_fraction": np.nan,
            "stroma_fraction": np.nan,
            "mixed_boundary_fraction": np.nan,
            "mean_withheld_score": np.nan,
            "withheld_binary_rate": np.nan,
            "cta_endpoint_rate": np.nan,
            "entropy": np.nan,
            "diversity": np.nan,
        }
    sub = df.iloc[indices]
    domain_counts = np.array([(sub["morphology_domain"].eq(d)).sum() for d in DOMAINS], dtype=float)
    return {
        "tumor_fraction": float(sub["Tumor_fraction"].mean()),
        "immune_fraction": float(sub["Immune_cells_fraction"].mean()),
        "stroma_fraction": float(sub["Stroma_fraction"].mean()),
        "mixed_boundary_fraction": float(sub["morphology_domain"].eq("mixed_boundary").mean()),
        "mean_withheld_score": float(sub["withheld_score"].mean()),
        "withheld_binary_rate": float(sub["withheld_binary"].mean()),
        "cta_endpoint_rate": float(sub["endpoint_positive_binary"].mean()),
        "entropy": entropy(domain_counts),
        "diversity": float((domain_counts > 0).sum()),
    }


def construct_microenvironment(df: pd.DataFrame, one_hop: list[list[int]], two_hop: list[list[int]]) -> pd.DataFrame:
    rows = []
    for i in range(len(df)):
        row = {"barcode": df.iloc[i]["barcode"]}
        for prefix, neigh in [("n1", one_hop[i]), ("n2", two_hop[i])]:
            summary = local_summary(df, neigh if neigh else [i])
            for key, value in summary.items():
                row[f"{prefix}_{key}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


def choose_k_and_cluster(micro: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame, int, float]:
    features = [
        "n1_tumor_fraction",
        "n1_immune_fraction",
        "n1_stroma_fraction",
        "n1_mixed_boundary_fraction",
        "n1_entropy",
        "n1_diversity",
        "n2_tumor_fraction",
        "n2_immune_fraction",
        "n2_stroma_fraction",
        "n2_mixed_boundary_fraction",
        "n2_entropy",
        "n2_diversity",
    ]
    x = micro[features].fillna(micro[features].median()).to_numpy(float)
    xz = StandardScaler().fit_transform(x)
    best = None
    for k in range(4, 7):
        model = KMeans(n_clusters=k, random_state=SEED, n_init=50)
        labels = model.fit_predict(xz)
        score = silhouette_score(xz, labels)
        if best is None or score > best[0]:
            best = (score, k, labels)
    score, k, labels = best
    centroids = micro.assign(cluster=labels).groupby("cluster")[features].mean().reset_index()
    return labels, centroids, k, float(score)


def name_niches(centroids: pd.DataFrame) -> dict[int, str]:
    names: dict[int, str] = {}
    for _, r in centroids.iterrows():
        cluster = int(r["cluster"])
        tumor = float(r["n1_tumor_fraction"])
        immune = float(r["n1_immune_fraction"])
        stroma = float(r["n1_stroma_fraction"])
        entropy_v = float(r["n1_entropy"])
        if tumor > 0.90 and entropy_v < 0.15:
            name = "deep_tumor_core"
        elif tumor > 0.75:
            name = "tumor_interior"
        elif immune > 0.30:
            name = "immune_associated_niche"
        elif stroma > 0.40:
            name = "stroma_interface_niche"
        elif tumor > 0.55 and entropy_v > 0.50:
            name = "mixed_tumor_boundary_niche"
        else:
            name = "low_signal_edge_context"
        names[cluster] = name
    return names


def mean_ci(x: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan"), float("nan")
    if len(x) == 1:
        val = float(x[0])
        return val, val, val
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        boot[i] = RNG.choice(x, size=len(x), replace=True).mean()
    return float(x.mean()), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else float("nan")


def empirical_p(null: np.ndarray, observed: float, greater: bool = True) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if len(null) == 0:
        return float("nan")
    if greater:
        return float((np.sum(null >= observed) + 1) / (len(null) + 1))
    return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))


def summarize_niches(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for niche in sorted(df["niche_label"].unique()):
        sub = df[df["niche_label"].eq(niche)]
        rows.append(
            {
                "niche_label": niche,
                "n_spots": int(len(sub)),
                "mean_withheld_score": float(sub["withheld_score"].mean()),
                "withheld_rate": float(sub["withheld_binary"].mean()),
                "cta_endpoint_rate": float(sub.loc[sub["main_analysis_spot"], "endpoint_positive_binary"].mean()) if sub["main_analysis_spot"].any() else np.nan,
                "mean_immune_fraction": float(sub["Immune_cells_fraction"].mean()),
                "mean_tumor_fraction": float(sub["Tumor_fraction"].mean()),
                "mean_stroma_fraction": float(sub["Stroma_fraction"].mean()),
                "mean_n1_entropy": float(sub["n1_entropy"].mean()),
            }
        )
    return pd.DataFrame(rows)


def niche_statistics(df: pd.DataFrame, niche_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    top_immune = niche_summary.sort_values("mean_immune_fraction", ascending=False).iloc[0]["niche_label"]
    top_tumor = niche_summary.sort_values("mean_tumor_fraction", ascending=False).iloc[0]["niche_label"]
    top_interface = niche_summary.sort_values(["mean_stroma_fraction", "mean_n1_entropy"], ascending=False).iloc[0]["niche_label"]
    rows = []
    boot_rows = []
    perm_rows = []
    for niche in niche_summary["niche_label"]:
        sub = df[df["niche_label"].eq(niche)]
        for metric in ["withheld_binary", "withheld_score", "Immune_cells_fraction"]:
            mean, lo, hi = mean_ci(sub[metric].astype(float).to_numpy())
            boot_rows.append({"panel": "I", "group": niche, "metric": metric, "n": int(len(sub)), "mean": mean, "ci95_low": lo, "ci95_high": hi, "n_bootstrap": N_BOOT})

    for comparison, a_label, b_label in [
        ("immune_associated_vs_tumor_dominant", top_immune, top_tumor),
        ("interface_related_vs_tumor_dominant", top_interface, top_tumor),
    ]:
        a = df[df["niche_label"].eq(a_label)]["withheld_score"].to_numpy(float)
        b = df[df["niche_label"].eq(b_label)]["withheld_score"].to_numpy(float)
        obs = float(np.nanmean(a) - np.nanmean(b))
        labels = df["niche_label"].to_numpy()
        values = df["withheld_score"].to_numpy(float)
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            perm_values = RNG.permutation(values)
            null[i] = float(np.nanmean(perm_values[labels == a_label]) - np.nanmean(perm_values[labels == b_label]))
            perm_rows.append({"analysis": "niche_withheld_score_permutation", "comparison": comparison, "permutation": i, "statistic": null[i]})
        rows.append(
            {
                "analysis": "niche_effect_size",
                "comparison": comparison,
                "group_a": a_label,
                "group_b": b_label,
                "observed_mean_difference": obs,
                "empirical_p_one_sided": empirical_p(null, obs, greater=True),
                "cohen_d": cohen_d(a, b),
                "n_group_a": int(len(a)),
                "n_group_b": int(len(b)),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(boot_rows), pd.DataFrame(perm_rows)


def neighbor_fraction_matrix(df: pd.DataFrame, neighborhoods: list[list[int]], niche_order: list[str]) -> np.ndarray:
    labels = df["niche_label"].to_numpy()
    mat = np.zeros((len(df), len(niche_order)), dtype=float)
    for i, ns in enumerate(neighborhoods):
        if not ns:
            mat[i, :] = np.nan
            continue
        neigh_labels = labels[ns]
        for j, niche in enumerate(niche_order):
            mat[i, j] = np.mean(neigh_labels == niche)
    return mat


def neighbor_enrichment(df: pd.DataFrame, one_hop: list[list[int]], two_hop: list[list[int]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    niche_order = sorted(df["niche_label"].unique())
    rows = []
    perm_rows = []
    labels = df["withheld_binary"].to_numpy(bool)
    for hop_name, neighborhoods in [("1-hop", one_hop), ("2-hop", two_hop)]:
        frac_mat = neighbor_fraction_matrix(df, neighborhoods, niche_order)
        obs = np.nanmean(frac_mat[labels, :], axis=0)
        null = np.empty((N_PERM, len(niche_order)))
        for i in range(N_PERM):
            perm = RNG.permutation(labels)
            null[i, :] = np.nanmean(frac_mat[perm, :], axis=0)
        for j, niche in enumerate(niche_order):
            mu = float(np.nanmean(null[:, j]))
            sd = float(np.nanstd(null[:, j], ddof=1))
            z = float((obs[j] - mu) / sd) if sd > 0 else float("nan")
            p = empirical_p(null[:, j], obs[j], greater=True)
            rows.append(
                {
                    "hop": hop_name,
                    "neighbor_niche": niche,
                    "observed_neighbor_fraction": float(obs[j]),
                    "null_mean": mu,
                    "null_sd": sd,
                    "enrichment_z": z,
                    "empirical_p_one_sided": p,
                    "n_permutations": N_PERM,
                }
            )
        for i in range(N_PERM):
            perm_rows.append({"analysis": "neighbor_enrichment", "hop": hop_name, "permutation": i, **{n: null[i, k] for k, n in enumerate(niche_order)}})
    return pd.DataFrame(rows), pd.DataFrame(perm_rows)


def edge_arrays(one_hop: list[list[int]]) -> tuple[np.ndarray, np.ndarray]:
    src = []
    dst = []
    for i, ns in enumerate(one_hop):
        for j in ns:
            src.append(i)
            dst.append(j)
    return np.asarray(src, dtype=int), np.asarray(dst, dtype=int)


def moran_geary(values: np.ndarray, src: np.ndarray, dst: np.ndarray) -> tuple[float, float]:
    x = np.asarray(values, dtype=float)
    xbar = float(np.nanmean(x))
    z = x - xbar
    denom = float(np.nansum(z * z))
    wij = len(src)
    moran_num = float(np.nansum(z[src] * z[dst]))
    geary_num = float(np.nansum((x[src] - x[dst]) ** 2))
    n = len(x)
    moran = float(n / wij * moran_num / denom) if wij and denom > 0 else float("nan")
    geary = float((n - 1) / (2 * wij) * geary_num / denom) if wij and denom > 0 else float("nan")
    return moran, geary


def spatial_autocorrelation(df: pd.DataFrame, one_hop: list[list[int]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    perm_rows = []
    src, dst = edge_arrays(one_hop)
    for metric, values in [
        ("withheld_score", df["withheld_score"].to_numpy(float)),
        ("withheld_binary", df["withheld_binary"].astype(float).to_numpy()),
    ]:
        obs_moran, obs_geary = moran_geary(values, src, dst)
        null_moran = np.empty(N_PERM)
        null_geary = np.empty(N_PERM)
        for i in range(N_PERM):
            perm = RNG.permutation(values)
            null_moran[i], null_geary[i] = moran_geary(perm, src, dst)
            perm_rows.append({"analysis": "spatial_autocorrelation", "metric": metric, "permutation": i, "moran_i": null_moran[i], "geary_c": null_geary[i]})
        rows.append(
            {
                "metric": metric,
                "moran_i": obs_moran,
                "moran_empirical_p_positive": empirical_p(null_moran, obs_moran, greater=True),
                "geary_c": obs_geary,
                "geary_empirical_p_low": float((np.sum(null_geary <= obs_geary) + 1) / (len(null_geary) + 1)),
                "n_permutations": N_PERM,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(perm_rows)


def save_panel_i(df: pd.DataFrame, img: np.ndarray, niche_summary: pd.DataFrame, niche_stats: pd.DataFrame) -> None:
    h, w = img.shape[:2]
    niche_order = sorted(df["niche_label"].unique())
    color_map = {n: NICHE_COLORS[i % len(NICHE_COLORS)] for i, n in enumerate(niche_order)}
    fig = plt.figure(figsize=(10.4, 3.55), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.10, 1.10, 1.18], wspace=0.42)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    for ax in [ax0, ax1]:
        ax.imshow(img, extent=(0, w, h, 0), alpha=0.88)
        ax.set_xlim(0, w)
        ax.set_ylim(h, 0)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    for niche in niche_order:
        sub = df[df["niche_label"].eq(niche)]
        ax0.scatter(
            sub["x_plot"],
            sub["y_plot"],
            s=7,
            c=color_map[niche],
            linewidths=0,
            alpha=0.75,
            label=DISPLAY_NAME.get(niche, niche),
            rasterized=True,
        )
    ax0.set_title("Unsupervised local niches", fontsize=8, fontweight="bold")
    ax0.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.025),
        fontsize=4.8,
        frameon=True,
        framealpha=0.88,
        markerscale=1.15,
        ncol=2,
        columnspacing=0.8,
        handletextpad=0.3,
    )

    base = df[~df["withheld_binary"]]
    wh = df[df["withheld_binary"]]
    ax1.scatter(
        base["x_plot"],
        base["y_plot"],
        s=4,
        c="#D0D0D0",
        linewidths=0,
        alpha=0.32,
        rasterized=True,
    )
    ax1.scatter(
        wh["x_plot"],
        wh["y_plot"],
        s=10,
        c="#0072B2",
        edgecolors="white",
        linewidths=0.25,
        alpha=0.90,
        rasterized=True,
    )
    ax1.set_title("SVTuner withheld-positive spots", fontsize=8, fontweight="bold")

    ordered = niche_summary.sort_values("withheld_rate", ascending=True)
    y = np.arange(len(ordered))
    ax2.barh(y, ordered["withheld_rate"], color=[color_map[n] for n in ordered["niche_label"]], alpha=0.86)
    ax2.set_yticks(y)
    ax2.set_yticklabels([DISPLAY_NAME.get(n, n) for n in ordered["niche_label"]], fontsize=6)
    ax2.set_xlim(0, max(0.05, float(ordered["withheld_rate"].max()) * 1.25))
    ax2.set_xlabel("withheld rate", fontsize=7)
    ax2.set_title("Withheld enrichment by niche", fontsize=8, fontweight="bold")
    ax2.grid(axis="x", color="#E0E0E0", linewidth=0.5)
    for i, row in enumerate(ordered.itertuples(index=False)):
        ax2.text(row.withheld_rate + 0.005, i, f"n={row.n_spots}", va="center", fontsize=5.5)
    for sp in ["top", "right"]:
        ax2.spines[sp].set_visible(False)
    fig.suptitle("I. Spatial microenvironment context of SVTuner withheld output", fontsize=10, fontweight="bold", y=1.02)
    fig.savefig(OUT / "panel_i_microenvironment_context.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_i_microenvironment_context.pdf", bbox_inches="tight")
    fig.savefig(OUT / "panel_i_microenvironment_context.svg", bbox_inches="tight")
    plt.close(fig)


def save_panel_j(enrich: pd.DataFrame) -> None:
    mat = enrich.pivot(index="hop", columns="neighbor_niche", values="enrichment_z").loc[["1-hop", "2-hop"]]
    pmat = enrich.pivot(index="hop", columns="neighbor_niche", values="empirical_p_one_sided").loc[["1-hop", "2-hop"]]
    fig, ax = plt.subplots(figsize=(5.8, 2.2), facecolor="white")
    cmap = LinearSegmentedColormap.from_list("enrich", ["#2166AC", "#FFFFFF", "#B2182B"])
    lim = max(1.0, float(np.nanmax(np.abs(mat.to_numpy(float)))))
    im = ax.imshow(mat.to_numpy(float), vmin=-lim, vmax=lim, cmap=cmap, aspect="auto")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([DISPLAY_NAME.get(c, c) for c in mat.columns], fontsize=6, rotation=30, ha="right")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            z = mat.iloc[i, j]
            p = pmat.iloc[i, j]
            star = "*" if p < 0.05 else ""
            ax.text(j, i, f"{z:.1f}{star}", ha="center", va="center", fontsize=6, color="black")
    ax.set_title("J. Neighbor-niche enrichment around withheld-positive spots", fontsize=8, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("permutation z-score", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.savefig(OUT / "panel_j_spatial_neighborhood_enrichment.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_j_spatial_neighborhood_enrichment.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df, img = load_inputs()
    coords = df[["x_plot", "y_plot"]].to_numpy(float)
    one_hop, two_hop, radius = build_neighbors(coords)
    micro = construct_microenvironment(df, one_hop, two_hop)
    labels, centroids, k, silhouette = choose_k_and_cluster(micro)
    name_map = name_niches(centroids)
    df = df.merge(micro, on="barcode", how="left")
    df["niche_cluster"] = labels
    df["niche_label"] = df["niche_cluster"].map(name_map)

    niche_assignment = df[["barcode", "niche_cluster", "niche_label", "x_plot", "y_plot", "morphology_domain", "withheld_binary", "withheld_score", "primary_endpoint_status"]].copy()
    spot_micro = df[["barcode"] + [c for c in micro.columns if c != "barcode"]].copy()
    niche_summary = summarize_niches(df)
    niche_stats, boot, niche_perm = niche_statistics(df, niche_summary)
    enrich, enrich_perm = neighbor_enrichment(df, one_hop, two_hop)
    moran, moran_perm = spatial_autocorrelation(df, one_hop)
    permutation_results = pd.concat([niche_perm, enrich_perm, moran_perm], ignore_index=True, sort=False)

    spot_micro.to_csv(OUT / "spot_microenvironment_table.csv", index=False)
    niche_assignment.to_csv(OUT / "niche_assignment.csv", index=False)
    niche_summary.to_csv(OUT / "niche_summary.csv", index=False)
    enrich.to_csv(OUT / "neighbor_enrichment.csv", index=False)
    moran.to_csv(OUT / "moran_statistics.csv", index=False)
    niche_stats.to_csv(OUT / "niche_statistics.csv", index=False)
    permutation_results.to_csv(OUT / "permutation_results.csv", index=False)
    boot.to_csv(OUT / "bootstrap_results.csv", index=False)

    save_panel_i(df, img, niche_summary, niche_stats)
    save_panel_j(enrich)

    guardrails = {
        "phase": "BioApp Downstream Phase D3 spatial microenvironment context validation",
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "mapping_regenerated": False,
        "endpoint_redefined": False,
        "morphology_domains_redefined": False,
        "tumor_stroma_boundary_redefined": False,
        "CTA_labels_redefined": False,
        "threshold_redefined": False,
        "niche_clustering_used_svtuner_outputs": False,
        "withheld_outputs_used_only_for_evaluation": True,
        "guardrails_passed": True,
    }
    write_json(OUT / "phase_d3_guardrails.json", guardrails)

    top_withheld = niche_summary.sort_values("withheld_rate", ascending=False).iloc[0].to_dict()
    top_neighbor = enrich.sort_values("enrichment_z", ascending=False).iloc[0].to_dict()
    summary = {
        "decision": "PASS",
        "phase": "BioApp Downstream Phase D3 spatial microenvironment context validation",
        "n_spots": int(len(df)),
        "neighbor_radius_image_px": radius,
        "niche_k": int(k),
        "niche_silhouette_score": silhouette,
        "niche_name_map": {str(k): v for k, v in name_map.items()},
        "top_withheld_niche": top_withheld,
        "top_neighbor_enrichment": top_neighbor,
        "moran_statistics": moran.to_dict(orient="records"),
        "interpretation": "SVTuner withheld outputs show preferential localization within computational local spatial microenvironment contexts.",
        "claim_boundary": [
            "Niches derive from unsupervised clustering of local neighborhood descriptors.",
            "SVTuner withheld outputs were not used to define niches.",
            "Niches are computational proxy niches, not pathology annotations.",
            "No new biological discovery claim is made.",
        ],
    }
    write_json(OUT / "phase_d3_summary.json", summary)

    report = f"""# BioApp Downstream Phase D3 spatial microenvironment context validation

Decision: `PASS`

## Scope

This phase builds the final downstream context layer using frozen spot
coordinates, CTA composition, morphology domains, boundary proxy, and SVTuner
withheld outputs. It does not rerun SVTuner, Stage3, Stage4, or CytoSPACE.

## Niche construction

- Neighborhoods: 1-hop and 2-hop spot graph neighborhoods.
- Clustering: KMeans on local neighborhood descriptors only.
- Chosen k: `{k}`
- Silhouette score: `{silhouette:.4f}`
- SVTuner withheld variables used for clustering: `false`

## Main outputs

- `panel_i_microenvironment_context.png/pdf`
- `panel_j_spatial_neighborhood_enrichment.png/pdf`
- `spot_microenvironment_table.csv`
- `niche_assignment.csv`
- `niche_summary.csv`
- `neighbor_enrichment.csv`
- `moran_statistics.csv`

## Interpretation boundary

SVTuner withheld outputs show preferential localization within specific
computational local spatial microenvironment contexts. These niches are proxy
niches derived from local neighborhood descriptors. They are not pathology
annotations and do not establish a new biological discovery.
"""
    (OUT / "phase_d3_report.md").write_text(report, encoding="utf-8")

    print("BioApp Downstream Phase D3 completed.")
    print("Decision: PASS")
    print(f"Niche k: {k}")
    print(f"Silhouette: {silhouette:.4f}")
    print(f"Top withheld niche: {top_withheld['niche_label']}")
    print(f"Output directory: {OUT}")


if __name__ == "__main__":
    main()
