# Engineering Quality Gates

The project should pass three engineering gates before manuscript-facing
handoff:

```bash
Rscript -e 'testthat::test_local()'
Rscript scripts/check_project_integrity.R
R CMD build . && R CMD check --no-manual FibroDynMix_0.0.0.9000.tar.gz
```

## What The Integrity Gate Checks

- Required figure packages exist for Figure 1-6 and the public real-data smoke
  supplement.
- PDF/SVG/PNG/TIFF/contact-sheet exports exist.
- PNG/contact-sheet dimensions are large enough for manuscript review.
- `figure_manifest.tsv` contains non-empty `primary_claim` and
  `claim_boundary`.
- `panel_source_data_manifest.tsv` points to existing source-data files.
- Main figure legends are present.
- `.DS_Store` files are absent.
- Large public data and generated analysis/figure folders are excluded from R
  package builds through `.Rbuildignore`.
- Package maturity files and submission-readiness audit outputs are present.

## Current Claim Discipline

The code and figures should distinguish implemented functionality from planned
functionality:

- implemented: raw-count NB likelihood, latent simplex optimization,
  ridge-penalized study and donor effects, simulation benchmark, bootstrap
  uncertainty, logistic-normal VI posterior skeleton and simulation calibration,
  transition flow, public raw-count demonstration, leave-study-out transfer
  benchmark, and public real-data transfer smoke;
- planned: full amortized or fully hierarchical variational posterior over all
  parameters and completed human disease-atlas or cross-atlas generalization.
