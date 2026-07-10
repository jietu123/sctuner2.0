# BioApp Main Figure V3.11 Top-Row A/B/C Scientific Visualization Redesign Report

## Purpose

V3.11 is a layout-only redesign of the top row. Panel A was redrawn as a four-step experimental design schematic, Panel B was changed from metric cards into an endpoint-concordance lollipop plot, and Panel C was harmonized as a publication-style burden-reduction bar chart. D/E/F biological content is preserved from V3.10.

## Inputs Used

- Base figure version: `V3.10`.
- V3.10 D/E/F tissue-overlay panels.
- Phase 8 spot-level SVTuner-vs-endpoint outputs.
- Phase 9 final biological-application audit.
- Coordinate registration recovery audit.
- Embedded image exported from `ST_data.RData::bcsa@images[['BCSA2TumB1']]@image`.

## Tissue Background

V3.11 continues to use the verified Seurat image-slot embedded tissue image. No raw external H&E, internet image, screenshot, or unverified Space Ranger bundle was used.

## Coordinate Mapping

```text
x = imagecol * lowres
y = imagerow * lowres
```

Lowres factor used: `0.0416`

Embedded image dimensions: `[600, 594, 3]`

Bounds check passed: `True`

## Panel A: Four-Step Experimental Design Schematic

Panel A was redrawn to remove text and arrow overlap and to make the experimental logic readable:

1. Frozen external CTA Immune endpoint, evaluation-only.
2. Immune-all-dropout reference with ST spots unchanged.
3. Baseline forced call versus SVTuner withheld call.
4. Frozen endpoint evaluation, Phase 9 PASS, biological evidence allowed.

Guardrail labels were added in compact form: endpoint locked, no endpoint tuning, no rerun, and no threshold change.

## Panel B: Endpoint-Concordance Visualization

The previous KPI-style metric cards were removed. Panel B is now a two-group withheld-call-rate comparison:

- CTA endpoint-negative withheld rate: `2.17%`
- CTA endpoint-positive withheld rate: `44.78%`

Supporting metrics are shown only as a small annotation:

- AUROC: `0.8305`
- AUPRC: `0.2837`
- enrichment delta: `0.1453`

This panel should be interpreted as endpoint-concordant withholding, not proof that every withheld spot is an immune cell.

## Panel C: Burden Reduction Style Harmonization

Panel C remains a burden-reduction bar chart using the frozen Phase 8 values:

- Baseline forced = `1.00`
- Prevented binary = `44.8%`
- Prevented continuous = `0.52`

The bars were narrowed and typography was aligned with Panel B.

## D2 Footprint-Only Redesign

D2 is intentionally shown as a footprint rather than a burden gradient, because the baseline forced burden is near-saturated / near-constant across the analyzed region, making a continuous or binned intensity map visually misleading.

- display mapping: footprint-only single-color display
- footprint fill: `#66BB6A`
- footprint alpha: `0.38`
- colorbar removed: `true`
- gradient removed: `true`
- quantile / rank bins removed: `true`
- exact constant values: `True`
- near-constant values: `False`
- unique value count: `3`

中文说明：由于 baseline forced burden 在分析区域内近似饱和、缺乏可解释的空间强弱梯度，因此 D2 以 footprint 形式展示其覆盖范围，而不是以连续热图或分箱强度图展示。

The D2 panel is `fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_11.png`.

The visual emphasis is the spatial footprint of baseline forced non-immune assignment, not nonexistent spatial intensity variation.

## D/E/F Preservation

D/E/F were preserved from V3.10. In the main figure, D uses the D2/D3 no-external-marker version so that:

- D1 shows the CTA endpoint.
- D2 shows the baseline forced-call footprint.
- D3 shows the SVTuner withheld score.
- D4 shows endpoint-only / withheld-only / endpoint+withheld overlap.

The D4 three-class logic is unchanged:
- endpoint only
- withheld only
- endpoint + withheld

D4 overlap semantics were preserved.

## Preserved Semantics

- D1 still shows the frozen CTA Immune endpoint.
- D2 now shows the baseline forced-call footprint under immune-all-dropout.
- D3 still shows SVTuner withheld score and keeps the V3.5/V3.2 blue logic.
- D4 still shows endpoint-only, withheld-only, and endpoint + withheld overlap categories.
- E/F local FOV panels keep the same interpretation and remain visualization-only.

D2 remains a baseline forced-call footprint, not a biological identity map and not a quantitative burden-gradient map. D3 remains a withheld-score display, not proof of immune-cell identity. D4 overlap remains a visual overlap category, not a newly computed biological metric.

## Overlap Preservation Check

- Whole section endpoint-only: `74`
- Whole section withheld-only: `89`
- Whole section endpoint + withheld: `60`
- Count consistency passed: `True`

These are figure-level display checks only, not new biological metrics.

## Unchanged Elements

No data, metrics, endpoint, thresholds, biological labels, or conclusions were changed. The Seurat embedded tissue image and coordinate mapping remain:

```text
x = imagecol * lowres
y = imagerow * lowres
```

## Recommended Visual Interpretation

The top row provides the experimental design, endpoint-concordance summary, and burden-reduction summary; the lower panels provide tissue-background spatial evidence. V3.11 changes only the A/B/C design and preserves D/E/F data and interpretation.

## Allowed Interpretation

SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune assignment burden relative to the baseline mapping.

## Disallowed Interpretations

- SVTuner discovered a new immune niche.
- SVTuner proves the biological identity of all withheld spots.
- SVTuner definitively improves all biological interpretations.
- SVTuner replaces spatial transcriptomics annotation.
- The withheld regions are confirmed immune cells without external validation.
- Any claim based on newly tuned thresholds or redefined endpoints.
