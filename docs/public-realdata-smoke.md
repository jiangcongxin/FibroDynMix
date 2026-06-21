# Public Real-Data Smoke Analysis

## Purpose

The public smoke analysis checks that FibroDynMix can run on real raw UMI
counts rather than only simulated matrices:

```text
public raw counts
-> prepared FibroDynMix input
-> NB optimizer
-> state composition
-> FPI and optional transition flow
-> fit diagnostics and marker-coverage audit
```

This is a smoke test, not a biological claim figure.

## Public Dataset

The first public source is:

```text
Hum et al. / Sebastian et al. breast fibroblast scRNA-seq raw counts
DOI: 10.6071/M3238R
Zenodo record: https://zenodo.org/records/3977255
Dryad record: https://datadryad.org/dataset/doi:10.6071/M3238R
```

The script uses two raw-count files from the public record:

```text
MT3_CAFs_raw.txt
Normal_mammary_fibroblasts_raw.txt
```

These are mouse Cell Ranger UMI count matrices. The public runner labels them as
`disease` and `normal` for technical smoke testing.

## Commands

Automatic public runner:

```bash
Rscript scripts/run_public_realdata_smoke.R \
  --data-dir=data/public_dryad_breast_fibroblast \
  --out=analysis/public_realdata_smoke \
  --max-cells=240 \
  --max-genes=700
```

If the public portal blocks command-line content downloads, download the two
files manually from the Zenodo or Dryad record and place them in:

```text
data/public_dryad_breast_fibroblast/
```

Then rerun the same command.

Generic local runner:

```bash
Rscript scripts/run_realdata_smoke.R \
  --counts=data/public_dryad_breast_fibroblast/dryad_breast_fibroblast_counts.rds \
  --metadata=data/public_dryad_breast_fibroblast/dryad_breast_fibroblast_metadata.tsv \
  --out=analysis/public_realdata_smoke \
  --condition-col=condition \
  --study-col=study_id \
  --donor-col=donor_id \
  --cell-id-col=cell_id
```

For h5ad inputs, first export to Matrix Market and TSV files:

```bash
python3 scripts/export_h5ad_to_fibrodynmix_inputs.py \
  --h5ad=input.h5ad \
  --out=data/h5ad_export \
  --layer=X \
  --cell-type-col=cell_type \
  --cell-type-regex=fibro
```

Then run:

```bash
Rscript scripts/run_realdata_smoke.R \
  --counts=data/h5ad_export/counts.mtx \
  --genes=data/h5ad_export/genes.tsv \
  --cells=data/h5ad_export/cells.tsv \
  --metadata=data/h5ad_export/metadata.tsv \
  --out=analysis/h5ad_realdata_smoke
```

## Expected Outputs

- `cell_state_weights.tsv`
- `state_composition.tsv`
- `cell_fpi.tsv`
- `transition_flow.tsv`, when normal and disease labels are available
- `transition_summary.tsv`
- `marker_coverage.tsv`
- `filter_summary.tsv`
- `nb_fit_diagnostics.tsv`
- `run_manifest.tsv`

## Claim Boundary

The public smoke analysis only demonstrates that the raw-count FibroDynMix
pipeline executes on a public real scRNA-seq count matrix. It should not be used
as a disease-mechanism result without full QC, sample-level replication,
cell-type curation, and sensitivity analysis.
