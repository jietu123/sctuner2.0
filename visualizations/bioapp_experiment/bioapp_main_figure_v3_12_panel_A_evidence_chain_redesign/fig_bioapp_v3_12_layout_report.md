# BioApp Main Figure V3.12 Panel A Evidence-Chain Schematic Redesign Report

## Purpose

V3.12 is a layout-only redesign of Panel A. V3.11 solved text overlap, but the A panel became too simplified and compressed the actual BioApp experiment into a generic four-step schematic. V3.12 expands Panel A into a compact biological-application evidence chain while preserving B/C/D/E/F.

## Inputs Used

- Base figure version: `V3.11`.
- Phase 8 frozen SVTuner-vs-endpoint evaluation outputs.
- Phase 9 final biological-application audit outputs.
- V3.11 B/C/D/E/F panel logic and accepted D2 footprint-only display.
- Verified Seurat image-slot tissue background: `ST_data.RData::bcsa@images[['BCSA2TumB1']]@image`.

## Panel A Redesign

Panel A now uses a five-block left-to-right evidence chain:

1. Frozen CTA endpoint as the evaluation anchor.
2. Immune-all-dropout reference perturbation with ST spots unchanged.
3. Two same-input outputs: baseline forced non-immune call and SVTuner withheld / unsupported call.
4. Phase 8 frozen endpoint comparison.
5. Phase 9 Golden Rules audit and evidence-allowed decision.

The middle block deliberately splits into two branches so that baseline and SVTuner are shown as outputs from the same immune-all-dropout input, rather than as unrelated steps.

## Endpoint Leakage Prevention

The CTA endpoint is shown only as a frozen evaluation anchor. The schematic does not imply that CTA endpoint labels were used to run SVTuner, tune thresholds, rerun mapping, or redefine the endpoint.

Guardrail labels shown in Panel A:

- endpoint locked
- no endpoint tuning
- no rerun
- no threshold change

## Preserved Panels

- Panel B is preserved as the V3.11 endpoint-concordant withholding plot.
- Panel C is preserved as the V3.11 forced-burden reduction plot.
- Panel D is preserved with D1 endpoint, D2 footprint-only baseline forced-call display, D3 SVTuner withheld score, and D4 endpoint/withheld overlap.
- Panels E/F are preserved as local FOV tissue-background views.

## D2 Preservation

D2 remains footprint-only. No gradient, colorbar, quantile binning, or rank display was reintroduced.

## Tissue Background and Coordinate Mapping

V3.12 continues to use the verified Seurat image-slot embedded tissue image. Coordinate mapping remains:

```text
x = imagecol * lowres
y = imagerow * lowres
```

Lowres factor used: `0.0416`

Bounds check passed: `True`

## Interpretation

Panel A summarizes the frozen endpoint-anchored biological application design; Panels B and C summarize endpoint concordance and burden reduction; Panels D/E/F provide tissue-background spatial evidence.

Allowed interpretation: Panel A explains the frozen biological-application evaluation design in which an immune-all-dropout reference perturbation creates a stress test, baseline forced calls and SVTuner withheld calls are compared against an independent frozen CTA Immune endpoint, and final Golden Rules audit authorizes bounded biological-application evidence.

Disallowed interpretations:

- SVTuner used the CTA endpoint to decide where to withhold.
- SVTuner discovered all immune spots.
- The endpoint was used for threshold tuning.
- Phase 8 or Phase 9 reran mapping.
- Any biological conclusion changed in V3.12.

## Guardrail Confirmation

- Data changed: `false`
- Metrics changed: `false`
- Endpoint redefined: `false`
- Threshold changed: `false`
- Conclusions changed: `false`
- Panel A evidence chain expanded: `true`
- Panel A text overlap removed: `true`
- B/C/D/E/F preserved: `true`
- D2 footprint-only preserved: `true`
