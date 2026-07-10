# BioApp Main Figure V3.10 D2 Light-Green Footprint Correction Report

## Purpose

V3.10 is a layout-only correction of V3.9. The intended color change was for D2, not D3. V3.10 therefore changes D2 to a light-green footprint, restores D3 to the previous blue withheld-score colormap, and preserves the D1/D4 marker-size reductions plus the additional D-panel export where D2 and D3 do not show external endpoint reference markers.

## Inputs Used

- Base figure version: `V3.9`.
- V3.9 D1/D4 marker refinement and extra D-panel export.
- Phase 8 spot-level SVTuner-vs-endpoint outputs.
- Phase 9 final biological-application audit.
- Coordinate registration recovery audit.
- Embedded image exported from `ST_data.RData::bcsa@images[['BCSA2TumB1']]@image`.

## Tissue Background

V3.10 continues to use the verified Seurat image-slot embedded tissue image. No raw external H&E, internet image, screenshot, or unverified Space Ranger bundle was used.

## Coordinate Mapping

```text
x = imagecol * lowres
y = imagerow * lowres
```

Lowres factor used: `0.0416`

Embedded image dimensions: `[600, 594, 3]`

Bounds check passed: `True`

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

The D2 panel is `fig_bioapp_panel_D2_baseline_forced_call_footprint_v3_10.png`.

The visual emphasis is the spatial footprint of baseline forced non-immune assignment, not nonexistent spatial intensity variation.

## D1/D2/D3/D4 Refinement

- D1 external endpoint markers were reduced to limit tissue-background occlusion.
- D2 footprint now uses a light-green footprint color: `#66BB6A`.
- D3 withheld score is restored to the original blue colormap: `#F7F7F7` to `#0072B2`.
- D4 markers were reduced again to limit overlap clutter.
- An extra D-panel export removes external endpoint reference markers from D2 and D3 only.

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

The D panel should be read as a tissue-background spatial comparison between the frozen external CTA Immune endpoint, baseline forced-call footprint, SVTuner withheld score, and endpoint-withheld overlap. V3.10 changes only D2 footprint color, D3 colormap restoration, marker size, and reference-marker visibility, not data or conclusions.

## Allowed Interpretation

SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune assignment burden relative to the baseline mapping.

## Disallowed Interpretations

- SVTuner discovered a new immune niche.
- SVTuner proves the biological identity of all withheld spots.
- SVTuner definitively improves all biological interpretations.
- SVTuner replaces spatial transcriptomics annotation.
- The withheld regions are confirmed immune cells without external validation.
- Any claim based on newly tuned thresholds or redefined endpoints.
