# NCEM Figure 2A-F reproduction

## Scope

This directory reproduces Figure 2A-F from:

Fischer et al., *Modeling intercellular communication in tissues using spatial graphs of cells*.

The analysis uses the downloaded human lymph-node scRNA-seq and Visium data, the
official cell2location preprocessing workflow, and the NCEM 0.1.5 statistical
definitions. No cell-type whitelist or manually supplied unsupported-cell-type
list is used.

## Panels

- Figure 2A: B-cell, FDC, and Mast pseudo-cell expression subclusters and mean
  cell-type abundance in their spatial neighborhoods.
- Figure 2B: directed cell-type coupling network. Edges require at least 200
  significant genes; width represents the L1 norm of significant effects.
- Figure 2C: sender effects for the B-cell receiver.
- Figure 2D: receiver effects for the FDC sender, retaining the tensor indexing
  behavior of the published NCEM 0.1.5 implementation.
- Figure 2E: B-cell/FDC sender-receiver volcano plot.
- Figure 2F: correlation of sender-effect profiles for the B-cell receiver.

## Reproducibility

Environment:

```text
conda environment: ncem_fig2_tf39
Python: 3.9
NCEM: 0.1.5
TensorFlow: 2.11.1
Scanpy: 1.9.3
```

Run:

```powershell
C:\Users\Aimyon\.conda\envs\ncem_fig2_tf39\python.exe `
  sctuner2.0\scripts\reproduce_ncem_fig2_human_lymph_node.py
```

Force recomputation:

```powershell
C:\Users\Aimyon\.conda\envs\ncem_fig2_tf39\python.exe `
  sctuner2.0\scripts\reproduce_ncem_fig2_human_lymph_node.py `
  --force-panel-a --force-statistics
```

The optimized OLS/Wald implementation was checked against the official NCEM
functions on a balanced subset. Maximum absolute differences were:

```text
design matrix: 0
OLS coefficient: 1.76e-14
p-value: 2.46e-13
q-value: 2.01e-13
```

## Data-version note

The current downloaded Visium input contains 4,039 spots, producing 64,624
cell-type-specific pseudo-cells. The historical output embedded in the official
notebook reports 4,035 spots and 64,560 pseudo-cells. Therefore, the reproduced
panel structure and statistical conclusions follow the official workflow, but
individual genes and exact values are not expected to be pixel-identical to the
archived paper figure.
