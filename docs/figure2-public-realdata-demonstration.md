# Figure 2 Public Real-Data Demonstration

## Figure-Level Claim

Figure 2 demonstrates that FibroDynMix can run on a balanced public raw-count
fibroblast subset and produce coherent model outputs:

```text
raw UMI count matrices
-> NB objective optimization
-> latent state composition
-> cell-level state mixtures
-> FPI
-> cross-sectional transition flow
```

## Claim Boundary

```text
This is a public real-data demonstration, not a disease-mechanism validation.
```

The data are pooled mouse breast fibroblast/CAF count matrices from DOI
`10.6071/M3238R`. The figure should not be interpreted as a replicated disease
biology result because donor-level replication, full QC, cell curation,
sensitivity analysis, and external validation are still required.

## Panels

- **A. Balanced public raw-count subset**
  - Condition-stratified subset used for the demonstration.

- **B. NB objective improvement**
  - Raw-count objective trace from the FibroDynMix NB optimizer.

- **C. Disease-minus-normal shift**
  - Difference in fitted state composition between disease and normal labels.

- **D. State composition**
  - Condition-level mean latent state weights.

- **E. FPI by condition**
  - Cell-level Fibroblast Plasticity Index.

- **F. Inferred normal-to-disease flow**
  - Cross-sectional optimal-transport state flow.

- **G. Cell-level latent mixtures**
  - Heatmap of fitted simplex weights ordered by condition and dominant state.

## Role in the Manuscript

Figure 2 bridges the architecture figure and simulation benchmarks by showing
that the current implementation runs on real public raw UMI count matrices.
The strong methodological claims should still come from simulation benchmarks,
calibration, uncertainty, and transition-flow validation.
