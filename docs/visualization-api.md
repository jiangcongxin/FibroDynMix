# Visualization API

FibroDynMix now exposes package-level plotting helpers for model outputs and
validation tables. These functions return `ggplot` objects; they do not write
files or impose a manuscript layout.

## Core Plots

- `plot_state_composition()` visualizes sample, donor, condition, or dataset
  state fractions from a long composition table.
- `plot_cell_state_heatmap()` visualizes cell-level simplex weights as a
  cell-by-state heatmap.
- `plot_transfer_diagnostics()` visualizes leave-dataset-out or leave-donor-out
  transfer convergence and related diagnostics.
- `plot_transition_flow()` visualizes optimal-transport state flow as a
  source-state by target-state heatmap.
- `plot_fpi_distribution()` visualizes Fibroblast Plasticity Index or entropy
  distributions across conditions, datasets, or donors.
- `plot_benchmark_rankings()` visualizes method performance across simulation
  stress scenarios.
- `plot_purity_qc()` visualizes fibroblast purity margins from the
  bioinformatics validation layer.
- `plot_pathway_enrichment()` visualizes state-program pathway enrichment.
- `plot_fibroblast_annotation()` visualizes FibroDynMix or user-supplied
  fibroblast annotations on Seurat embeddings.
- `plot_fibroblast_marker_dot()` visualizes marker expression by FibroDynMix
  state or another fibroblast annotation column.

## Example

```r
library(FibroDynMix)

composition <- read.delim(
  "analysis/independent_geo_gse167339_validation/multi_public_dataset_state_composition.tsv",
  check.names = FALSE
)

p <- plot_state_composition(
  composition,
  group_col = "dataset_id",
  state_col = "state",
  value_col = "composition"
)

p
```

## Seurat Annotation Plots

```r
plot_fibroblast_annotation(
  seurat_object,
  reduction = "umap",
  annotation_col = "fibrodynmix_dominant_state"
)

plot_fibroblast_marker_dot(
  seurat_object,
  features = c("COL1A1", "DCN", "ACTA2", "CXCL12", "HLA-DRA"),
  group_col = "fibrodynmix_dominant_state"
)
```

The plotting API is intentionally thin. Statistical quantities such as state
weights, uncertainty, transition flow, FPI, purity margins, and enrichment
statistics are computed upstream by model or validation functions; plotting
functions only check columns, reshape simple wide state-weight tables when
needed, and render reviewable figures.
