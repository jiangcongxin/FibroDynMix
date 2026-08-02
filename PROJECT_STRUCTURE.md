# Project Structure Map

A quick map of what each top-level item is, and what is deliverable vs scratch.
Last updated: 2026-08-03.

## R package (the shippable core)

| Path | What it is |
|---|---|
| `DESCRIPTION`, `NAMESPACE`, `LICENSE` | R package metadata |
| `R/` | Package source (21 functions: NB model, VI, transfer, transition flow, plots) |
| `man/` | Generated `.Rd` documentation (do not edit by hand) |
| `tests/` | testthat suite (359 tests) |
| `inst/` | Installed extras (incl. `CITATION`) |
| `NEWS.md`, `README.md` | Release notes and package readme |
| `_pkgdown.yml` | pkgdown site config |
| `cran-comments.md` | CRAN-style check notes |

The released source snapshot is available from the GitHub
[`v0.1.0` release](https://github.com/jiangcongxin/FibroDynMix/releases/tag/v0.1.0)
(commit `d945dd1`) and Zenodo version DOI
[10.5281/zenodo.21761171](https://doi.org/10.5281/zenodo.21761171).

## Research evidence (build-ignored, but load-bearing for the paper)

| Path | What it is |
|---|---|
| `analysis/` | Benchmark and validation outputs; many are required by the integrity gate |
| `figures/` | Working figure packages (source data, manifests, legends, QC) |
| `data/` | Public input datasets (Dryad breast, GSE246215, GSE167339). Required by scripts |
| `projects/` | Per-dataset subprojects with fitted Seurat objects (~7.5 GB). Referenced by transfer scripts |
| `scripts/` | Run scripts for every analysis (see `docs/reproducibility-runbook.md`) |
| `docs/` | Method/validation write-ups and runbook |
| `MT3_CAFs_raw.txt`, `Normal_mammary_fibroblasts_raw.txt` | Root-path public count inputs required by smoke/runtime scripts |

## Submission deliverable

| Path | What it is |
|---|---|
| `submission/iscience_submission_package_final_v4/` | Current synchronized iScience submission-freeze candidate: manuscript, figures, front matter, cover letter, source-data manifests, and QC records. |
| `submission/iscience_review_package/` | Historical review package; do not mix its files into the current submission candidate. |

## Scratch / regenerable (moved aside)

| Path | What it is |
|---|---|
| `archive/` | Items moved out of the top level during cleanup. See `archive/ARCHIVE_MANIFEST.md`. Nothing deleted; everything restorable. Build-ignored. |
| `.venv-scvi` | (now in `archive/regenerable_env/`) regenerable Python env for scVI baseline |

## Where to look first

- Want the current submission candidate? → `submission/iscience_submission_package_final_v4/`
- Want to re-run analyses? → `scripts/` + `docs/reproducibility-runbook.md`
- Want the validation gate? → `scripts/check_project_integrity.R`
- Submission status? → `submission/iscience_submission_package_final_v4/08_quality_control/final_freeze_status.md`

## Note on heavy data

`data/` (~580 MB) and `projects/` (~7.5 GB) are large but cannot be relocated
without breaking script paths and the integrity gate (they are referenced by
relative paths). They are already excluded from the package build via
`.Rbuildignore`. To move them out cleanly later, refactor the referencing
scripts to accept a configurable base path first.
