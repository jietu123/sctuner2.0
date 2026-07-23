# Task 4D Final Decision

## 1. Decision

`CONDITIONAL PASS`

The formal Stage3B implementation, execution paths, solver function, two formal
solver environments and mathematical behavior are confirmed. The decision is
conditional because the mainline Stage1 logs do not contain a run-local
`sessionInfo()` or `packageVersion()` snapshot.

## 2. Seurat

### Mainline Stage1

- Formal use: confirmed by 10 successful Stage1 logs.
- Environment pin: R `4.5.2`, Seurat `5.3.1`, SeuratObject `5.2.0`.
- Exact formal-run version: **not independently recoverable**.

### CTA biological application

- R: `4.5.1` (`x86_64-w64-mingw32`).
- Seurat: `5.3.1`.
- SeuratObject: `5.2.0`.
- Evidence status: direct formal package/runtime records.

## 3. Stage3B NNLS

- Implementation: `src/stages/stage3b_st_unsupported.py`.
- Function: `fit_nonnegative_mixtures`.
- Caller: `run_stage3b`.
- Solver: `scipy.optimize.nnls`.
- Fitting: per spot with design `B.T` (`n_genes x n_types`).
- Intercept: none.
- Regularisation: none.
- Explicit tolerance/maxiter: none; SciPy defaults apply.
- Fallback: none.
- Coefficients: sum-normalised only when the sum is positive.
- Zero solution: retained as zeros; reconstruction is zero.
- Failure: exception propagates; no silent replacement.

## 4. Formal solver versions

- 5%/10% composite simulation routes: Python `3.10.19`, NumPy `1.26.4`, SciPy `1.11.4`.
- High-resolution/BioApp base routes: Python `3.13.5`, NumPy `2.1.0`, SciPy `1.15.3`.

## 5. Mathematical consistency

`PASS_WITH_DOCUMENTED_DEFAULTS_AND_PREPROCESSING_CONTEXT`

The objective orientation, non-negative constraint, conditional coefficient
normalisation, reconstruction and residual match the supplied Methods
definition. The manual Methods wording should define `x_i` and `B` as the
linearised, row-composition-normalised profiles and should not claim a custom
tolerance or fallback.

## 6. Remaining unresolved items

- Mainline Stage1 exact run-local R/Seurat/SeuratObject versions.
- Whether the optional Seurat enrichment backend generated the formal Fig. 2 panel.

## 7. Guardrails

```text
Formal manuscript accessed: false
Formal manuscript modified: false
Formal bibliography accessed: false
Formal bibliography modified: false
Experimental code modified: false
Experimental outputs modified: false
Figures modified: false
Formal stages rerun: false
GitHub repository access required: false
```
