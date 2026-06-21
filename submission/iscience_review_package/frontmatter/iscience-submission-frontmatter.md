# iScience Submission Front Matter Package

## Purpose

This file collects the graphical abstract plan, eTOC blurb, and cover-letter fit
paragraph for the iScience-facing FibroDynMix manuscript. The wording is aligned
with the current manuscript title, Summary, Results headings, STAR Methods, and
figure claim boundaries.

## Official Constraints Used

- iScience is a broad Cell Press journal for original research across the life,
  physical, earth, social, and health sciences.
- iScience article guidance advises concise articles and places detailed methods
  in STAR Methods.
- Cell Press graphical abstracts should be a single-panel image that gives an
  immediate view of the paper's main finding. The working technical target is a
  1200 x 1200 px square, Arial 12-16 pt, and TIFF, PDF, or JPG for final upload.
- The graphical abstract should be distinct from main-text figures, use sparse
  labels, read left-to-right or top-to-bottom, emphasize the new findings, and
  avoid embedding the eTOC blurb or Highlights in the image.

## Graphical Abstract

### Working Visual Claim

FibroDynMix turns public fibroblast raw UMI count matrices into auditable
simplex state mixtures, with simulation truth and public datasets used to bound
state recovery, transfer, and scar-context interpretation.

### Figure Strategy

- Study style: bioinformatics/omics method figure.
- Format: one square panel, 1200 x 1200 px, left-to-right flow.
- Evidence boundary: method-first computational study; no clinical diagnostic,
  causal disease-mechanism, or lineage-tracing claim.
- Source artwork draft:
  `figures/graphical_abstract/graphical_abstract.svg`.

### Panel Flow

1. **Input layer**: raw UMI count matrices from public fibroblast datasets.
2. **Model layer**: negative-binomial count model with weak marker orientation,
   simplex state weights, and study/donor terms.
3. **Output layer**: state mixtures, uncertainty, transfer, and
   cross-sectional transition-flow/FPI.
4. **Evidence layer**: simulation truth, public transfer, marker-gradient
   recovery, and directional scar-module shifts.

### Text Allowed Inside Image

Keep text sparse. Recommended labels:

- Raw UMI counts
- FibroDynMix
- Simplex state mixtures
- Uncertainty + transfer
- Cross-sectional state flow
- Simulation truth
- Public fibroblast datasets

### Text to Avoid Inside Image

- Exact RMSE, AUC, q values, or other data items.
- eTOC blurb or Highlights.
- "Validated", "confirmed", "diagnostic", "therapeutic", or causal disease
  language.
- Dense method details such as optimizer settings.

### Visual Style

- Layout: three main vertical blocks plus a thin evidence band at the bottom.
- Color: neutral background; teal for raw counts, blue for model fitting,
  green for outputs, and gray for claim boundaries.
- Annotation: direct arrows from input to model to outputs; no causal arrows for
  biological disease mechanisms.
- Typography: Arial; keep image labels large enough for thumbnail reading.

### Alt Text

Raw fibroblast UMI count matrices enter FibroDynMix, which fits a
negative-binomial model with weak marker orientation, simplex state weights, and
study/donor terms. The fitted model returns state mixtures, uncertainty,
transfer, and cross-sectional transition-flow summaries, with simulation truth
and public fibroblast datasets defining the claim boundary.

## eTOC Blurb

### Primary Version

FibroDynMix estimates fibroblast state mixtures directly from raw UMI counts,
combining weak marker orientation, study/donor effects, uncertainty, transfer,
and cross-sectional transition flow to benchmark state recovery and test bounded
biological utility in public fibroblast datasets.

Word count: 35.

### Shorter Backup Version

FibroDynMix models raw single-cell fibroblast counts as simplex state mixtures,
separating likelihood fit from state recovery in simulations and supporting
bounded transfer, marker-gradient recovery, and directional scar-module signals
across public fibroblast datasets.

Word count: 32.

## Cover-Letter Fit Paragraph

We are submitting "FibroDynMix models fibroblast state plasticity as raw-count
latent mixtures" to iScience because the study combines a computational model,
single-cell fibroblast biology, public data reuse, and reproducible software.
FibroDynMix addresses a common limitation in fibroblast state analysis: marker
scores are useful for annotation but do not estimate mixed state weights under a
raw count-generating model or separate state structure from study and donor
effects. The manuscript reports a negative-binomial mixture model, simulation
benchmarks with known latent truth, validation-aware model selection, public
raw-count execution, transfer across public fibroblast datasets, and explicit
source-data and claim-boundary tracking. We believe this method-first,
source-backed study will interest iScience readers who use single-cell data to
quantify cell-state plasticity across disease contexts.

## One-Sentence Fit Statement

FibroDynMix fits iScience as a reproducible single-cell computational method
that links raw-count modeling, benchmarked state recovery, and bounded
public-data validation for fibroblast state plasticity.

