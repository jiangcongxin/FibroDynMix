# Supplementary Bootstrap Uncertainty Figure

## Command

```bash
Rscript scripts/make_supplementary_uncertainty_figure.R
```

## Outputs

- `figures/supplementary_uncertainty/exports/supplementary_uncertainty.png`
- `figures/supplementary_uncertainty/exports/supplementary_uncertainty.pdf`
- `figures/supplementary_uncertainty/exports/supplementary_uncertainty.svg`
- `figures/supplementary_uncertainty/exports/supplementary_uncertainty.tiff`
- `figures/supplementary_uncertainty/source_data/fig5_sample_composition_uncertainty.tsv`
- `figures/supplementary_uncertainty/source_data/fig5_cell_state_uncertainty.tsv`
- `figures/supplementary_uncertainty/source_data/fig5_entropy_uncertainty.tsv`
- `figures/supplementary_uncertainty/source_data/fig5_uncertainty_vs_entropy.tsv`
- `figures/supplementary_uncertainty/source_data/fig5_marker_program_stability.tsv`

## Figure-Level Claim

The supplementary uncertainty figure summarizes the first uncertainty layer for
FibroDynMix:

```text
Cell-bootstrap replicates provide uncertainty summaries for sample-level state
composition, cell entropy, state weights, and marker-program stability.
```

The claim boundary is explicit:

```text
These intervals are bootstrap uncertainty summaries from simulated data, not
full Bayesian posterior credible intervals.
```

## Panels

- **A. Sample-level state composition intervals**
  - Bootstrap mean and 95% interval for selected fibroblast states.

- **B. Entropy uncertainty by cell**
  - Cell-level entropy interval width against mean entropy.

- **C. Marker-program bootstrap stability**
  - Top genes per state ranked by mean absolute `beta`, with bootstrap
    intervals.

- **D. State-weight uncertainty across mixed states**
  - Mean state-weight interval width evaluated against mean entropy.

## Role in the Manuscript

This supplementary figure answers the reviewer question:

```text
Does FibroDynMix expose uncertainty rather than only point-estimate state labels?
```

It is a pragmatic uncertainty layer. The eventual full Bayesian/VI version
should replace or complement this with posterior credible intervals.
