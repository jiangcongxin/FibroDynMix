# GSE243716 Scar Fibroblast Annotation

This package case study documents a compact public-data workflow for applying
FibroDynMix to fibroblasts from pathological scar tissue.

The dataset is GEO `GSE243716`, a human scRNA-seq study comparing one keloid
sample and one hypertrophic scar sample from the same patient. It is smaller
than multi-donor keloid atlases and is therefore a better package tutorial
dataset. The commands below are written as a reproducible analysis template and
are not run during package checks because they download public raw matrices.

The complete runnable project script is:

```sh
Rscript projects/scar_fibro_gse243716/scripts/run_scar_fibro_analysis.R
```

## Dataset

- GEO accession: `GSE243716`
- Title: Single-cell analysis reveals transcriptomic difference between human
  hypertrophic scar and keloid
- Samples:
  - `GSM7794710`: keloid scRNA-seq
  - `GSM7794711`: hypertrophic scar scRNA-seq
- Supplementary file: `GSE243716_RAW.tar`, provided as MTX/TSV files

## Download Raw Data

```r
project <- file.path("projects", "scar_fibro_gse243716")
dir.create(file.path(project, "data"), recursive = TRUE, showWarnings = FALSE)

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243716&format=file",
  file.path(project, "data", "GSE243716_RAW.tar"),
  mode = "wb"
)
```

The archive contains two 10x-style matrices:

- `GSM7794710_K-barcodes.tsv.gz`
- `GSM7794710_K-features.tsv.gz`
- `GSM7794710_K-matrix.mtx.gz`
- `GSM7794711_H-barcodes.tsv.gz`
- `GSM7794711_H-features.tsv.gz`
- `GSM7794711_H-matrix.mtx.gz`

Extract them with:

```r
untar(file.path(project, "data", "GSE243716_RAW.tar"), list = TRUE)
untar(
  file.path(project, "data", "GSE243716_RAW.tar"),
  exdir = file.path(project, "data", "raw_matrices")
)
```

## Conventional Single-cell Workflow

Read the two samples, create a Seurat object, and run standard fibroblast
analysis steps. If the public archive includes author metadata with cell types,
use that metadata to retain fibroblasts. Otherwise, annotate fibroblasts by
canonical markers such as `DCN`, `LUM`, `COL1A1`, `COL1A2`, `PDGFRA`, and
`COL14A1`.

```r
library(Matrix)
library(Seurat)
library(FibroDynMix)

read_prefixed_10x <- function(prefix, condition, sample_id) {
  raw_dir <- file.path(project, "data", "raw_matrices")
  counts <- Matrix::readMM(gzfile(file.path(raw_dir, paste0(prefix, "-matrix.mtx.gz"))))
  features <- read.delim(gzfile(file.path(raw_dir, paste0(prefix, "-features.tsv.gz"))), header = FALSE)
  barcodes <- read.delim(gzfile(file.path(raw_dir, paste0(prefix, "-barcodes.tsv.gz"))), header = FALSE)
  rownames(counts) <- make.unique(features[[2]])
  colnames(counts) <- paste(sample_id, barcodes[[1]], sep = "_")
  metadata <- data.frame(condition = condition, sample_id = sample_id, row.names = colnames(counts))
  list(counts = counts, metadata = metadata)
}

keloid <- read_prefixed_10x("GSM7794710_K", "keloid", "GSM7794710")
hs <- read_prefixed_10x("GSM7794711_H", "hypertrophic_scar", "GSM7794711")
counts <- do.call(cbind, list(keloid$counts, hs$counts))
metadata <- rbind(keloid$metadata, hs$metadata)

object <- CreateSeuratObject(counts = counts, meta.data = metadata)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
object <- subset(object, subset = nFeature_RNA >= 500 & nFeature_RNA <= 6000 & percent.mt <= 20)
object <- NormalizeData(object)
object <- FindVariableFeatures(object, nfeatures = 2500)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 30)
object <- FindNeighbors(object, dims = 1:20)
object <- FindClusters(object, resolution = 0.4)
object <- RunUMAP(object, dims = 1:20)
```

If author cell-type labels are unavailable, one lightweight approach is to
retain clusters with broad fibroblast marker expression and low epithelial,
endothelial, immune, and melanocyte marker expression. For formal analysis,
manual inspection and marker plots should be retained in the project outputs.

## FibroDynMix Marker Priors

Use a compact human fibroblast state marker prior:

```r
fibro_markers <- list(
  resident = c("DCN", "LUM", "COL14A1", "PDGFRA", "PI16"),
  inflammatory = c("IL6", "CXCL12", "CXCL14", "CCL2", "CXCL2"),
  myofibroblast = c("ACTA2", "TAGLN", "MYL9", "TPM2", "CNN1"),
  `ECM-remodeling` = c("COL1A1", "COL1A2", "FN1", "POSTN", "MMP2"),
  `antigen-presenting` = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"),
  `IFN-stress` = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1")
)
```

## Fit and Transfer FibroDynMix

For a tutorial-sized dataset, fit directly on fibroblasts or on a balanced
subset, then transfer the fitted state program back to all fibroblasts. The
transfer call uses the package's warm-start, chunking, and optional parallel
support.

```r
fibro_object <- object
fibro_counts <- GetAssayData(fibro_object, assay = "RNA", layer = "counts")
fibro_metadata <- fibro_object[[]]
fibro_metadata$cell_id <- rownames(fibro_metadata)

prepared <- prepare_fibrodynmix_data(
  counts = fibro_counts,
  cell_metadata = fibro_metadata,
  marker_index = fibro_markers,
  cell_id_col = "cell_id",
  study_col = "condition",
  donor_col = "sample_id",
  min_cells_per_gene = 0,
  min_counts_per_gene = 0
)

fit <- fit_fibrodynmix_prepared(
  prepared,
  fit_study_effect = TRUE,
  fit_donor_effect = FALSE,
  n_outer = 1,
  initializer_args = list(n_iter = 2),
  maxit_beta = 20,
  maxit_z = 35
)

warm_start <- score_marker_baseline(
  counts = prepared$counts,
  marker_index = prepared$marker_index,
  library_size = prepared$library_size
)$z_pred

transfer <- fit_fibrodynmix_transfer(
  counts = prepared$counts,
  fit = fit,
  library_size = prepared$library_size,
  z_init = warm_start[, rownames(fit$beta_hat), drop = FALSE],
  chunk_size = 500,
  parallel = TRUE,
  n_workers = 4,
  maxit_z = 75
)
```

Attach results back to Seurat:

```r
fibro_object <- add_fibrodynmix_to_seurat(
  fibro_object,
  fit = transfer,
  cells = rownames(transfer$z_hat),
  prefix = "fibrodynmix_"
)
```

## Annotation Confidence

Evaluate annotation confidence before interpreting dominant states:

```r
weights <- data.frame(
  cell_id = rownames(transfer$z_hat),
  transfer$z_hat,
  check.names = FALSE
)

evaluation <- evaluate_fibrodynmix_annotation(
  weights,
  state_cols = names(fibro_markers),
  metadata = fibro_metadata,
  metadata_cell_col = "cell_id",
  cluster_col = "seurat_clusters",
  condition_col = "condition",
  expression = GetAssayData(fibro_object, assay = "RNA", layer = "data"),
  marker_index = fibro_markers
)

evaluation$state_summary
evaluation$marker_support
plot_marker_support(evaluation$marker_support)
plot_cluster_state_agreement(evaluation$cluster_agreement)
```

The `state_summary` table also reports support labels:

- `support_label`: `supported`, `exploratory`, or `unsupported`
- `formal_ready`: whether cell count, confidence, and marker support pass the
  configured thresholds
- `support_reason`: concise reason when a state remains exploratory

## Recommended Figures

```r
DimPlot(fibro_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(fibro_object, reduction = "umap", group.by = "fibrodynmix_dominant_state")
FeaturePlot(fibro_object, reduction = "umap", features = "fibrodynmix_max_weight")
FeaturePlot(fibro_object, reduction = "umap", features = "fibrodynmix_normalized_entropy")
```

The project script writes the main outputs to:

- `projects/scar_fibro_gse243716/results/fibro_state_support_summary.tsv`
- `projects/scar_fibro_gse243716/results/fibro_marker_support.tsv`
- `projects/scar_fibro_gse243716/results/fibro_transfer_cell_diagnostics.tsv`
- `projects/scar_fibro_gse243716/figures/fibro_clusters_and_fibrodynmix_umap.png`

## Interpretation Boundary

This compact dataset is useful for a package tutorial and API demonstration, not
for definitive population-level scar biology. Because it contains only one
keloid and one hypertrophic scar sample, condition differences should be treated
as example workflow outputs. Dominant FibroDynMix states should be interpreted
with marker support, cluster agreement, state confidence, and transfer
diagnostics.
