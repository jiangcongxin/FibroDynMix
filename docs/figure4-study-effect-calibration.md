# Figure 4 Study-Effect Calibration

## Figure-Level Claim

Figure 4 is a simulation-calibration figure for the first hierarchical
extension of FibroDynMix. Its claim boundary is:

```text
Ridge-penalized study effects can be calibrated under batch-confounded
simulation to balance state recovery, raw-count likelihood, study-effect
magnitude, and optimizer stability.
```

This is not a real-cohort biological result and does not establish disease
mechanisms.

## Panels

- **A. State recovery under batch confounding**
  - State-weight RMSE across study-effect penalties and marker penalties.
  - Dashed line shows the no-study-effect NB baseline.

- **B. Raw-count NB objective**
  - Best NB objective across the same penalty grid.
  - Shows the likelihood tradeoff when stronger shrinkage is imposed.

- **C. Study-effect magnitude**
  - L2 norm of the fitted study-by-gene effect.
  - Confirms larger `study_l2` shrinks the study effect.

- **D. Penalty tradeoff score**
  - Weighted selection score combining RMSE, NB objective, study-effect
    magnitude, and rollback rate.
  - Circled point marks the selected penalty pair.

## Current Selected Default

The current deterministic sensitivity run selects:

```text
study_l2 = 5
marker_l2 = 0.05
```

This selection is simulation-calibrated and should be revisited when the
benchmark size, confounding strength, or real-data target changes.
