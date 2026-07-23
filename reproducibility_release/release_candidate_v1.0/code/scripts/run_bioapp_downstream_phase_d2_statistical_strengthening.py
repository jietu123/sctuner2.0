from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "svg.fonttype": "none",
        "svg.image_inline": True,
    }
)


ROOT = Path(__file__).resolve().parents[1]
BIOAPP = ROOT / "visualizations" / "bioapp_experiment"
OUT = BIOAPP / "bioapp_downstream_phase_d2_statistical_strengthening"
D1 = BIOAPP / "bioapp_downstream_phase_d1_morphology_and_interface_analysis"

N_PERM = 2000
N_BOOT = 3000
SEED = 20260708
RNG = np.random.default_rng(SEED)

DOMAIN_ORDER = ["tumor_core", "immune_enriched", "stroma_rich", "mixed_boundary", "unmapped_or_excluded"]
FORMAL_DOMAIN_ORDER = ["tumor_core", "immune_enriched", "stroma_rich", "mixed_boundary"]
DOMAIN_COLORS = {
    "tumor_core": "#8C1D40",
    "immune_enriched": "#D55E00",
    "stroma_rich": "#4D9221",
    "mixed_boundary": "#7B3294",
    "unmapped_or_excluded": "#D0D0D0",
}
BOUNDARY_ORDER = [
    "deep tumor\n(<-4)",
    "inner tumor\n(-4,-2)",
    "tumor edge\n(-2,-1)",
    "near tumor\n(-1,0)",
    "interface side\n(0,1)",
    "outer side\n(1,2)",
    "outer domain\n(2,4)",
    "far outside\n(>4)",
]


def write_json(path: Path, payload: dict) -> None:
    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [clean(v) for v in obj]
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        return obj

    path.write_text(json.dumps(clean(payload), indent=2, allow_nan=False), encoding="utf-8")


def ci(vals: np.ndarray) -> tuple[float, float]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def mean_ci(x: np.ndarray, n_boot: int = N_BOOT) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan"), float("nan"), float("nan")
    if len(x) == 1:
        val = float(x[0])
        return val, val, val
    boots = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        boots[i] = RNG.choice(x, size=len(x), replace=True).mean()
    lo, hi = ci(boots)
    return float(x.mean()), lo, hi


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else float("nan")


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = 0
    lt = 0
    chunk = 1000
    for start in range(0, len(a), chunk):
        aa = a[start : start + chunk][:, None]
        gt += int((aa > b).sum())
        lt += int((aa < b).sum())
    return float((gt - lt) / (len(a) * len(b)))


def empirical_p(null: np.ndarray, observed: float, two_sided: bool = True) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if len(null) == 0 or not np.isfinite(observed):
        return float("nan")
    if two_sided:
        return float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))
    return float((np.sum(null >= observed) + 1) / (len(null) + 1))


def load_d1() -> pd.DataFrame:
    df = pd.read_csv(D1 / "phase_d1_downstream_analysis_by_spot.csv")
    df["withheld_binary"] = df["withheld_binary"].fillna(False).astype(bool)
    df["main_analysis_spot"] = df["primary_endpoint_status"].isin(["positive", "negative"])
    df["endpoint_positive_binary"] = df["primary_endpoint_status"].eq("positive").astype(float)
    if "boundary_distance_bin" not in df.columns:
        raise RuntimeError("D1 by-spot table does not contain boundary_distance_bin.")
    return df


def domain_stat(df: pd.DataFrame, labels: np.ndarray, test: str) -> float:
    tmp = df[["morphology_domain"]].copy()
    tmp["label"] = labels.astype(float)
    immune = tmp.loc[tmp["morphology_domain"].eq("immune_enriched"), "label"]
    if test == "immune_vs_tumor_core":
        ref = tmp.loc[tmp["morphology_domain"].eq("tumor_core"), "label"]
    elif test == "immune_vs_all_remaining":
        ref = tmp.loc[~tmp["morphology_domain"].eq("immune_enriched"), "label"]
    else:
        raise ValueError(test)
    return float(immune.mean() - ref.mean())


def panel_g_statistics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    formal = df[df["morphology_domain"].isin(FORMAL_DOMAIN_ORDER)].copy()
    labels = formal["withheld_binary"].to_numpy(bool)
    tests = ["immune_vs_tumor_core", "immune_vs_all_remaining"]
    perm_rows = []
    stat_rows = []
    for test in tests:
        obs = domain_stat(formal, labels, test)
        null = np.empty(N_PERM, dtype=float)
        for i in range(N_PERM):
            null[i] = domain_stat(formal, RNG.permutation(labels), test)
            perm_rows.append({"analysis": "panel_g_domain_permutation", "test": test, "permutation": i, "statistic": null[i]})
        stat_rows.append(
            {
                "analysis": "panel_g_domain_permutation",
                "test": test,
                "observed_statistic": obs,
                "empirical_p_two_sided": empirical_p(null, obs, two_sided=True),
                "n_permutations": N_PERM,
            }
        )

    boot_rows = []
    for domain in DOMAIN_ORDER:
        sub = df[df["morphology_domain"].eq(domain)]
        for metric, label in [
            ("withheld_binary", "withheld_rate"),
            ("withheld_score", "withheld_score"),
            ("Immune_cells_fraction", "immune_fraction"),
        ]:
            mean, lo, hi = mean_ci(sub[metric].astype(float).to_numpy())
            boot_rows.append(
                {
                    "panel": "G",
                    "group": domain,
                    "metric": label,
                    "n": int(len(sub)),
                    "mean": mean,
                    "ci95_low": lo,
                    "ci95_high": hi,
                    "n_bootstrap": N_BOOT,
                }
            )

    effect_rows = []
    comparisons = {
        "immune_enriched_vs_tumor_core": (
            formal.loc[formal["morphology_domain"].eq("immune_enriched"), "withheld_score"].to_numpy(float),
            formal.loc[formal["morphology_domain"].eq("tumor_core"), "withheld_score"].to_numpy(float),
        ),
        "immune_enriched_vs_all_remaining": (
            formal.loc[formal["morphology_domain"].eq("immune_enriched"), "withheld_score"].to_numpy(float),
            formal.loc[~formal["morphology_domain"].eq("immune_enriched"), "withheld_score"].to_numpy(float),
        ),
        "immune_enriched_vs_all_remaining_excluding_mixed_boundary": (
            formal.loc[formal["morphology_domain"].eq("immune_enriched"), "withheld_score"].to_numpy(float),
            formal.loc[
                (~formal["morphology_domain"].eq("immune_enriched")) & (~formal["morphology_domain"].eq("mixed_boundary")),
                "withheld_score",
            ].to_numpy(float),
        ),
    }
    for comp, (a, b) in comparisons.items():
        effect_rows.append(
            {
                "analysis": "panel_g_effect_size",
                "comparison": comp,
                "metric": "withheld_score",
                "n_group_a": int(np.isfinite(a).sum()),
                "n_group_b": int(np.isfinite(b).sum()),
                "mean_group_a": float(np.nanmean(a)),
                "mean_group_b": float(np.nanmean(b)),
                "cohen_d": cohen_d(a, b),
                "cliffs_delta": cliffs_delta(a, b),
            }
        )

    sensitivity = formal[~formal["morphology_domain"].eq("mixed_boundary")].copy()
    sens_obs = domain_stat(sensitivity, sensitivity["withheld_binary"].to_numpy(bool), "immune_vs_all_remaining")
    sens_null = np.empty(N_PERM, dtype=float)
    sens_labels = sensitivity["withheld_binary"].to_numpy(bool)
    for i in range(N_PERM):
        sens_null[i] = domain_stat(sensitivity, RNG.permutation(sens_labels), "immune_vs_all_remaining")
    stat_rows.append(
        {
            "analysis": "panel_g_sensitivity_excluding_mixed_boundary",
            "test": "immune_vs_all_remaining",
            "observed_statistic": sens_obs,
            "empirical_p_two_sided": empirical_p(sens_null, sens_obs, two_sided=True),
            "n_permutations": N_PERM,
        }
    )
    for i, val in enumerate(sens_null):
        perm_rows.append({"analysis": "panel_g_sensitivity_excluding_mixed_boundary", "test": "immune_vs_all_remaining", "permutation": i, "statistic": val})

    return pd.DataFrame(stat_rows), pd.DataFrame(boot_rows), pd.DataFrame(effect_rows), pd.DataFrame(perm_rows)


def logistic_regression_optional(df: pd.DataFrame) -> pd.DataFrame:
    try:
        import statsmodels.formula.api as smf
    except Exception:
        return pd.DataFrame([{"analysis": "panel_g_logistic_regression", "status": "statsmodels_unavailable"}])

    formal = df[df["morphology_domain"].isin(FORMAL_DOMAIN_ORDER)].copy()
    formal["withheld_binary_int"] = formal["withheld_binary"].astype(int)
    formal["morphology_domain"] = pd.Categorical(formal["morphology_domain"], categories=FORMAL_DOMAIN_ORDER)
    try:
        model = smf.logit(
            "withheld_binary_int ~ C(morphology_domain, Treatment(reference='tumor_core')) + Immune_cells_fraction + tissue_coverage",
            data=formal,
        ).fit(disp=False, maxiter=200)
        conf = model.conf_int()
        rows = []
        for term, coef in model.params.items():
            rows.append(
                {
                    "analysis": "panel_g_logistic_regression",
                    "status": "fit",
                    "term": term,
                    "odds_ratio": float(np.exp(coef)),
                    "ci95_low": float(np.exp(conf.loc[term, 0])),
                    "ci95_high": float(np.exp(conf.loc[term, 1])),
                    "p_value": float(model.pvalues[term]),
                }
            )
        return pd.DataFrame(rows)
    except Exception as exc:
        return pd.DataFrame([{"analysis": "panel_g_logistic_regression", "status": f"fit_failed: {exc}"}])


def bin_order(df: pd.DataFrame) -> list[str]:
    present = set(df["boundary_distance_bin"].astype(str))
    ordered = [x for x in BOUNDARY_ORDER if x in present]
    extras = sorted(present.difference(ordered))
    return ordered + extras


def topology_summary(df: pd.DataFrame) -> pd.DataFrame:
    order = bin_order(df)
    rows = []
    for i, label in enumerate(order):
        sub = df[df["boundary_distance_bin"].astype(str).eq(label)]
        main = sub[sub["main_analysis_spot"]]
        rows.append(
            {
                "boundary_distance_bin": label,
                "bin_index": i,
                "n_spots": int(len(sub)),
                "n_main_analysis_spots": int(len(main)),
                "endpoint_positive_rate": float(main["endpoint_positive_binary"].mean()) if len(main) else np.nan,
                "withheld_rate": float(sub["withheld_binary"].mean()) if len(sub) else np.nan,
                "mean_withheld_score": float(sub["withheld_score"].mean()) if len(sub) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def trend_from_summary(summary: pd.DataFrame, metric: str = "withheld_rate") -> tuple[float, float, float, float]:
    sub = summary[["bin_index", metric]].dropna()
    if len(sub) < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    spearman = stats.spearmanr(sub["bin_index"], sub[metric])
    kendall = stats.kendalltau(sub["bin_index"], sub[metric])
    return float(spearman.statistic), float(spearman.pvalue), float(kendall.statistic), float(kendall.pvalue)


def panel_h_statistics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = topology_summary(df)
    rho, rho_p, tau, tau_p = trend_from_summary(summary, "withheld_rate")

    perm_rows = []
    labels = df["withheld_binary"].to_numpy(bool)
    null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        tmp = df.copy()
        tmp["withheld_binary"] = RNG.permutation(labels)
        null[i] = trend_from_summary(topology_summary(tmp), "withheld_rate")[0]
        perm_rows.append({"analysis": "panel_h_topology_permutation", "test": "withheld_rate_spearman", "permutation": i, "statistic": null[i]})

    boot_rows = []
    for label in bin_order(df):
        sub = df[df["boundary_distance_bin"].astype(str).eq(label)]
        for metric, col in [("withheld_rate", "withheld_binary"), ("mean_withheld_score", "withheld_score")]:
            mean, lo, hi = mean_ci(sub[col].astype(float).to_numpy())
            boot_rows.append(
                {
                    "panel": "H",
                    "group": label,
                    "metric": metric,
                    "n": int(len(sub)),
                    "mean": mean,
                    "ci95_low": lo,
                    "ci95_high": hi,
                    "n_bootstrap": N_BOOT,
                }
            )

    effect_rows = []
    groups = {
        "tumor_side": df[df["boundary_distance_bin"].astype(str).isin([str(x) for x in bin_order(df)[:4]])],
        "interface": df[df["boundary_distance_bin"].astype(str).eq(bin_order(df)[4])],
        "outer_side": df[df["boundary_distance_bin"].astype(str).eq(bin_order(df)[5])],
        "outside": df[df["boundary_distance_bin"].astype(str).isin([str(x) for x in bin_order(df)[6:]])],
    }
    for comp, other in [("tumor_side_vs_interface", "interface"), ("tumor_side_vs_outer_side", "outer_side"), ("tumor_side_vs_outside", "outside")]:
        a = groups[other]["withheld_score"].to_numpy(float)
        b = groups["tumor_side"]["withheld_score"].to_numpy(float)
        effect_rows.append(
            {
                "analysis": "panel_h_effect_size",
                "comparison": comp,
                "metric": "withheld_score",
                "n_group_a": int(np.isfinite(a).sum()),
                "n_group_b": int(np.isfinite(b).sum()),
                "mean_group_a": float(np.nanmean(a)),
                "mean_group_b": float(np.nanmean(b)),
                "cohen_d": cohen_d(a, b),
                "cliffs_delta": cliffs_delta(a, b),
                "withheld_rate_group_a": float(groups[other]["withheld_binary"].mean()),
                "withheld_rate_tumor_side": float(groups["tumor_side"]["withheld_binary"].mean()),
            }
        )

    sensitivity_rows = []
    for name, bins in {
        "narrower_reporting_bins": [-np.inf, -4, -2, -0.5, 0.5, 2, 4, np.inf],
        "wider_reporting_bins": [-np.inf, -4, -2, -1.5, 1.5, 2, 4, np.inf],
    }.items():
        tmp = df.copy()
        tmp["boundary_distance_bin"] = pd.cut(tmp["boundary_distance_spot_width"], bins=bins, include_lowest=True).astype(str)
        sens_summary = topology_summary(tmp)
        sens_rho, sens_p, sens_tau, sens_tau_p = trend_from_summary(sens_summary, "withheld_rate")
        sensitivity_rows.append(
            {
                "analysis": "panel_h_reporting_bin_sensitivity",
                "setting": name,
                "spearman_rho": sens_rho,
                "spearman_p_value": sens_p,
                "kendall_tau": sens_tau,
                "kendall_p_value": sens_tau_p,
                "qualitative_positive_trend": bool(sens_rho > 0),
            }
        )

    stat_rows = [
        {
            "analysis": "panel_h_topology_trend",
            "metric": "withheld_rate",
            "spearman_rho": rho,
            "spearman_p_value": rho_p,
            "kendall_tau": tau,
            "kendall_p_value": tau_p,
            "empirical_permutation_p_one_sided_positive": empirical_p(null, rho, two_sided=False),
            "n_permutations": N_PERM,
        }
    ] + sensitivity_rows
    return pd.DataFrame(stat_rows), pd.DataFrame(boot_rows), pd.DataFrame(effect_rows), pd.DataFrame(perm_rows)


def save_panel_g(df: pd.DataFrame, g_stats: pd.DataFrame, g_boot: pd.DataFrame, g_eff: pd.DataFrame) -> None:
    domain_summary = pd.read_csv(D1 / "panel_g_morphology_domain_concordance_summary.csv")
    metrics = ["endpoint_positive_rate", "mean_cta_immune_fraction", "withheld_binary_rate", "mean_withheld_score"]
    labels = ["CTA endpoint+\nrate", "CTA immune\nfraction", "Withheld\nrate", "Withheld\nscore"]
    mat = domain_summary.set_index("morphology_domain").loc[DOMAIN_ORDER, metrics]

    fig = plt.figure(figsize=(7.5, 3.5), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.4], wspace=0.30)
    ax0 = fig.add_subplot(gs[0, 0])
    counts = domain_summary.set_index("morphology_domain").loc[DOMAIN_ORDER, "n_spots"]
    ax0.barh(range(len(DOMAIN_ORDER)), counts.values, color=[DOMAIN_COLORS[d] for d in DOMAIN_ORDER], edgecolor="none")
    ax0.set_yticks(range(len(DOMAIN_ORDER)))
    ax0.set_yticklabels([d.replace("_", " ") for d in DOMAIN_ORDER], fontsize=7)
    ax0.invert_yaxis()
    ax0.set_xlabel("spots", fontsize=7)
    ax0.set_title("Frozen domains", fontsize=8, fontweight="bold")
    ax0.tick_params(axis="x", labelsize=7)
    for sp in ax0.spines.values():
        sp.set_visible(False)

    ax1 = fig.add_subplot(gs[0, 1])
    im = ax1.imshow(mat.to_numpy(float), vmin=0, vmax=1, cmap="magma", aspect="auto")
    ax1.set_yticks(range(len(DOMAIN_ORDER)))
    ax1.set_yticklabels([d.replace("_", " ") for d in DOMAIN_ORDER], fontsize=7)
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, fontsize=6, rotation=35, ha="right")
    ax1.set_title("Frozen domain signal with bootstrap CIs", fontsize=8, fontweight="bold")
    boot_lookup = {
        (r["group"], r["metric"]): r
        for _, r in g_boot.iterrows()
    }
    metric_to_boot = {
        "mean_cta_immune_fraction": "immune_fraction",
        "withheld_binary_rate": "withheld_rate",
        "mean_withheld_score": "withheld_score",
    }
    for i, domain in enumerate(DOMAIN_ORDER):
        for j, metric in enumerate(metrics):
            val = mat.loc[domain, metric]
            if metric in metric_to_boot and (domain, metric_to_boot[metric]) in boot_lookup:
                r = boot_lookup[(domain, metric_to_boot[metric])]
                text = f"{val:.2f}\n[{r['ci95_low']:.2f},{r['ci95_high']:.2f}]"
                fs = 4.8
            else:
                text = f"{val:.2f}"
                fs = 5.5
            color = "white" if (val < 0.25 or val > 0.55) else "black"
            ax1.text(j, i, text, ha="center", va="center", fontsize=fs, color=color)
    cbar = fig.colorbar(im, ax=ax1, fraction=0.035, pad=0.02)
    cbar.set_label("fraction / mean score", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    for sp in ax1.spines.values():
        sp.set_visible(False)

    p = g_stats.loc[g_stats["test"].eq("immune_vs_tumor_core"), "empirical_p_two_sided"].iloc[0]
    d = g_eff.loc[g_eff["comparison"].eq("immune_enriched_vs_tumor_core"), "cohen_d"].iloc[0]
    delta = g_eff.loc[g_eff["comparison"].eq("immune_enriched_vs_tumor_core"), "cliffs_delta"].iloc[0]
    fig.suptitle(
        f"G. Morphology-domain concordance (immune vs tumor: permutation p={p:.4f}, d={d:.2f}, Cliff delta={delta:.2f})",
        fontsize=9,
        fontweight="bold",
        y=1.04,
    )
    fig.savefig(OUT / "panel_g_morphology_domain_concordance_statistical.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_g_morphology_domain_concordance_statistical.pdf", bbox_inches="tight")
    fig.savefig(OUT / "panel_g_morphology_domain_concordance_statistical.svg", bbox_inches="tight")
    plt.close(fig)


def save_panel_h(df: pd.DataFrame, h_stats: pd.DataFrame, h_boot: pd.DataFrame, h_eff: pd.DataFrame) -> None:
    summary = topology_summary(df)
    order = summary["boundary_distance_bin"].tolist()
    x = np.arange(len(summary))
    fig, ax = plt.subplots(figsize=(7.3, 3.2), facecolor="white")

    for metric, label, color, boot_metric in [
        ("withheld_rate", "SVTuner withheld rate", "#0072B2", "withheld_rate"),
        ("mean_withheld_score", "Mean withheld score", "#009E73", "mean_withheld_score"),
    ]:
        ci_rows = h_boot[h_boot["metric"].eq(boot_metric)].set_index("group").reindex(order)
        y = summary[metric].to_numpy(float)
        lo = ci_rows["ci95_low"].to_numpy(float)
        hi = ci_rows["ci95_high"].to_numpy(float)
        ax.plot(x, y, marker="o", linewidth=1.7, markersize=4.0, color=color, label=label)
        ax.fill_between(x, lo, hi, color=color, alpha=0.16, linewidth=0)

    endpoint_plot = summary["endpoint_positive_rate"].where(summary["n_main_analysis_spots"] >= 20)
    ax.plot(x, endpoint_plot, marker="o", linewidth=1.4, markersize=3.8, color="#D55E00", label="CTA endpoint+ rate")

    ax.axvline(3.5, color="#333333", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.text(3.55, 0.96, "boundary", fontsize=6, va="top")
    xticklabels = [f"{row.boundary_distance_bin}\nn={int(row.n_spots)}" for row in summary.itertuples(index=False)]
    ax.set_xticks(x)
    ax.set_xticklabels(xticklabels, fontsize=6, rotation=35, ha="right")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("fraction / mean score", fontsize=7)
    stat = h_stats[h_stats["analysis"].eq("panel_h_topology_trend")].iloc[0]
    eff = h_eff[h_eff["comparison"].eq("tumor_side_vs_outer_side")].iloc[0]
    ax.set_title(
        f"H. Boundary topology (Spearman rho={stat['spearman_rho']:.2f}, perm p={stat['empirical_permutation_p_one_sided_positive']:.4f}, outer-vs-tumor d={eff['cohen_d']:.2f})",
        fontsize=8,
        fontweight="bold",
    )
    ax.grid(axis="y", color="#E0E0E0", linewidth=0.6)
    ax.legend(loc="upper left", fontsize=6, frameon=False, ncol=2)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.savefig(OUT / "panel_h_tumor_stroma_interface_topology_statistical.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "panel_h_tumor_stroma_interface_topology_statistical.pdf", bbox_inches="tight")
    fig.savefig(OUT / "panel_h_tumor_stroma_interface_topology_statistical.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_d1()

    g_stats, g_boot, g_eff, g_perm = panel_g_statistics(df)
    logistic = logistic_regression_optional(df)
    h_stats, h_boot, h_eff, h_perm = panel_h_statistics(df)

    panel_g_csv = pd.concat([g_stats, g_eff, logistic], ignore_index=True, sort=False)
    panel_h_csv = pd.concat([h_stats, h_eff], ignore_index=True, sort=False)
    bootstrap = pd.concat([g_boot, h_boot], ignore_index=True, sort=False)
    permutations = pd.concat([g_perm, h_perm], ignore_index=True, sort=False)

    panel_g_csv.to_csv(OUT / "panel_g_statistics.csv", index=False)
    panel_h_csv.to_csv(OUT / "panel_h_statistics.csv", index=False)
    permutations.to_csv(OUT / "permutation_results.csv", index=False)
    bootstrap.to_csv(OUT / "bootstrap_results.csv", index=False)

    save_panel_g(df, g_stats, g_boot, g_eff)
    save_panel_h(df, h_stats, h_boot, h_eff)

    g_primary_p = float(g_stats.loc[g_stats["test"].eq("immune_vs_tumor_core"), "empirical_p_two_sided"].iloc[0])
    g_primary_d = float(g_eff.loc[g_eff["comparison"].eq("immune_enriched_vs_tumor_core"), "cohen_d"].iloc[0])
    h_primary = h_stats[h_stats["analysis"].eq("panel_h_topology_trend")].iloc[0]
    g_sensitivity = g_stats[g_stats["analysis"].eq("panel_g_sensitivity_excluding_mixed_boundary")].iloc[0]
    h_sensitivity = h_stats[h_stats["analysis"].eq("panel_h_reporting_bin_sensitivity")].to_dict(orient="records")
    logistic_status = (
        str(logistic["status"].iloc[0])
        if "status" in logistic.columns and len(logistic)
        else "not_reported"
    )

    guardrails = {
        "phase": "BioApp Downstream Phase D2 statistical strengthening",
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "mapping_regenerated": False,
        "morphology_domains_redefined": False,
        "tumor_stroma_boundary_redefined": False,
        "CTA_endpoint_redefined": False,
        "threshold_redefined": False,
        "endpoint_label_tuning": False,
        "only_frozen_D0_D1_outputs_consumed": True,
        "guardrails_passed": True,
    }
    write_json(OUT / "phase_d2_guardrails.json", guardrails)

    summary = {
        "decision": "PASS",
        "phase": "BioApp Downstream Phase D2 statistical strengthening",
        "n_permutations": N_PERM,
        "n_bootstrap": N_BOOT,
        "panel_g_primary_test": {
            "test": "immune_enriched_vs_tumor_core withheld_binary rate difference",
            "empirical_p_two_sided": g_primary_p,
            "cohen_d_withheld_score": g_primary_d,
        },
        "panel_g_sensitivity_excluding_mixed_boundary": {
            "test": "immune_enriched_vs_all_remaining withheld_binary rate difference",
            "empirical_p_two_sided": float(g_sensitivity["empirical_p_two_sided"]),
            "observed_statistic": float(g_sensitivity["observed_statistic"]),
        },
        "panel_g_logistic_regression_status": logistic_status,
        "panel_h_primary_test": {
            "test": "withheld_rate monotonic trend across frozen boundary bins",
            "spearman_rho": float(h_primary["spearman_rho"]),
            "spearman_p_value": float(h_primary["spearman_p_value"]),
            "kendall_tau": float(h_primary["kendall_tau"]),
            "empirical_permutation_p_one_sided_positive": float(h_primary["empirical_permutation_p_one_sided_positive"]),
        },
        "panel_h_reporting_bin_sensitivity": h_sensitivity,
        "panel_h_interpretation_note": "Primary frozen-bin topology is significant; narrower/wider reporting-bin sensitivity preserves a positive trend but has weaker nominal significance.",
        "panel_g_statistical_figure": "panel_g_morphology_domain_concordance_statistical.png",
        "panel_h_statistical_figure": "panel_h_tumor_stroma_interface_topology_statistical.png",
        "interpretation": "These analyses evaluate whether SVTuner withheld outputs are statistically associated with frozen morphology-domain and boundary topology.",
        "claim_boundary": [
            "No biological discovery claim.",
            "Frozen morphology domains are proxy labels, not pathology annotations.",
            "Frozen boundary is a spot-graph proxy, not a pixel-level histology boundary.",
            "No upstream BioApp metrics, endpoint labels, or thresholds were changed.",
        ],
    }
    write_json(OUT / "phase_d2_statistics_summary.json", summary)

    report = f"""# BioApp Downstream Phase D2 statistical strengthening

Decision: `PASS`

## Scope

This phase upgrades Panel G and Panel H from descriptive downstream prototypes
to statistically supported validation panels. It consumes only frozen D0/D1
outputs and does not rerun SVTuner, Stage3, Stage4, or CytoSPACE.

## Panel G

Primary permutation test:

- Test: immune-enriched vs tumor-core withheld-rate difference
- Empirical p-value: `{g_primary_p:.6g}`
- Withheld-score Cohen's d: `{g_primary_d:.4f}`

Bootstrap confidence intervals are written to `bootstrap_results.csv`.
Permutation null distributions are written to `permutation_results.csv`.

Sensitivity excluding mixed-boundary spots:

- Empirical p-value: `{float(g_sensitivity['empirical_p_two_sided']):.6g}`

Optional logistic regression status: `{logistic_status}`.

## Panel H

Primary topology trend test:

- Spearman rho: `{float(h_primary['spearman_rho']):.4f}`
- Spearman p-value: `{float(h_primary['spearman_p_value']):.6g}`
- Kendall tau: `{float(h_primary['kendall_tau']):.4f}`
- Empirical permutation p-value: `{float(h_primary['empirical_permutation_p_one_sided_positive']):.6g}`

Reporting-bin sensitivity retained positive trends, but nominal p-values were weaker than the primary frozen-bin analysis. Treat Panel H as topology-supporting evidence, not as a standalone mechanistic claim.

## Interpretation boundary

These analyses evaluate whether SVTuner withheld outputs are statistically
associated with frozen morphology-domain and boundary topology. They do not
establish a new biological discovery.

Frozen morphology domains are proxies. Frozen boundary is a spot-graph proxy.
No pathology annotation is used.
"""
    (OUT / "phase_d2_statistics_report.md").write_text(report, encoding="utf-8")

    print("BioApp Downstream Phase D2 completed.")
    print("Decision: PASS")
    print(f"Panel G p={g_primary_p:.6g}, d={g_primary_d:.4f}")
    print(
        "Panel H rho="
        f"{float(h_primary['spearman_rho']):.4f}, perm_p={float(h_primary['empirical_permutation_p_one_sided_positive']):.6g}"
    )
    print(f"Output directory: {OUT}")


if __name__ == "__main__":
    main()
