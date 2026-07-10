# BioApp Downstream Phase D2 statistical strengthening

Decision: `PASS`

## Scope

This phase upgrades Panel G and Panel H from descriptive downstream prototypes
to statistically supported validation panels. It consumes only frozen D0/D1
outputs and does not rerun SVTuner, Stage3, Stage4, or CytoSPACE.

## Panel G

Primary permutation test:

- Test: immune-enriched vs tumor-core withheld-rate difference
- Empirical p-value: `0.00049975`
- Withheld-score Cohen's d: `1.3153`

Bootstrap confidence intervals are written to `bootstrap_results.csv`.
Permutation null distributions are written to `permutation_results.csv`.

Sensitivity excluding mixed-boundary spots:

- Empirical p-value: `0.00049975`

Optional logistic regression status: `fit_failed: Singular matrix`.

## Panel H

Primary topology trend test:

- Spearman rho: `0.8503`
- Spearman p-value: `0.00747141`
- Kendall tau: `0.7638`
- Empirical permutation p-value: `0.0069965`

Reporting-bin sensitivity retained positive trends, but nominal p-values were weaker than the primary frozen-bin analysis. Treat Panel H as topology-supporting evidence, not as a standalone mechanistic claim.

## Interpretation boundary

These analyses evaluate whether SVTuner withheld outputs are statistically
associated with frozen morphology-domain and boundary topology. They do not
establish a new biological discovery.

Frozen morphology domains are proxies. Frozen boundary is a spot-graph proxy.
No pathology annotation is used.
