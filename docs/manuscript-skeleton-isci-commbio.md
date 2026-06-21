# FibroDynMix Manuscript Skeleton

## Working Title

FibroDynMix: a raw-count generative mixture framework for fibroblast state
plasticity, uncertainty, and cross-cohort transfer

## Target Journal Framing

The strongest current framing is a method-first manuscript for iScience or
Communications Biology. The biological framing should remain conservative:
FibroDynMix is presented as a probabilistic framework for fibroblast state
plasticity, not as a completed disease atlas, cross-atlas generalization study,
or lineage-tracing study.

## Abstract Draft

Single-cell studies often represent fibroblast phenotypes using fixed marker
scores, which can obscure mixed or transitional cell states and confound
biological variation with cohort effects. We introduce FibroDynMix, a
raw-count negative-binomial generative framework that models fibroblast states
as latent compositional programs on a simplex. The current implementation
combines latent state optimization, weak marker-prior orientation, ridge
study/donor effects, posterior uncertainty summaries, transition-flow
estimation, and cross-cohort transfer diagnostics. In simulations with known
latent truth, FibroDynMix quantifies recovery of state mixtures, robustness to
batch confounding, rare transition detection, marker-program recovery, and
posterior interval calibration against marker-scoring and NMF/topic baselines.
On public mouse breast fibroblast raw UMI counts, the workflow executes
end-to-end and supports auditable state composition, plasticity,
transition-flow, transfer smoke analyses, and registry-driven multi-public
leave-dataset-out transfer diagnostics. An independent GSE246215 human
fibroblast count-matrix validation extends this real-data check to four
cancer-type fibroblast subsets, while GSE167339, GSE156326, and GSE181316
provide bounded external biological-utility checks. The
package includes source-data-backed figure generation, release metadata, and
R CMD check validation. FibroDynMix provides a reproducible method-development
foundation for modeling fibroblast plasticity beyond direct gene-set scoring.

## Results Outline

### Result 1. FibroDynMix Defines a Raw-Count Generative Architecture

Primary figure: Figure 1

Core message:

- FibroDynMix models raw UMI counts under an NB likelihood.
- Fibroblast states are latent simplex variables, not direct marker scores.
- Study/donor effects, weak marker priors, transition flow, uncertainty, and
  transfer are explicit model layers.

Claim boundary:

- Figure 1 is an architecture and evidence map, not standalone biological
  validation.

### Result 2. Public Raw-Count Smoke Execution Demonstrates End-to-End Feasibility

Primary figures: Figure 2 and public real-data smoke figure

Core message:

- FibroDynMix runs on public raw count matrices.
- The workflow produces NB diagnostics, state composition, FPI, and
  cross-sectional transition-flow outputs.

Claim boundary:

- Public mouse pooled-count analyses are technical smoke tests and should not be
  presented as disease-mechanism discovery.

### Result 3. Simulations Quantify State-Mixture Recovery and VI Behavior

Primary figure: Figure 3
Primary tables: `analysis/extended_method_benchmark/`

Core message:

- Benchmark scenarios include continuous mixtures, discrete states, batch
  confounding, and rare transitions.
- FibroDynMix is evaluated against known latent state truth and compared with
  marker scoring plus a mature-package NMF/topic baseline.
- The VI posterior mean and calibrated interval diagnostics are included in the
  simulation benchmark.
- Marker scoring remains competitive in marker-aligned simulations; FibroDynMix
  should be positioned by its generative likelihood, hierarchy, uncertainty,
  transfer, and transition outputs rather than uniform RMSE dominance.

Claim boundary:

- Simulation results validate method behavior under known generative truth; they
  do not establish real-tissue transitions.

### Result 4. Study-Effect Calibration Addresses Cohort Confounding

Primary figure: Figure 4

Core message:

- Ridge-penalized study effects reduce confounding sensitivity while preserving
  state-mixture recovery.
- The selected penalty is justified by RMSE/objective/effect-size tradeoffs.

Claim boundary:

- The penalty selection is simulation-calibrated and should be re-evaluated for
  new real cohorts.

### Result 5. Independent Human Fibroblast Datasets Validate Biological Utility

Primary figure: Figure 5

Core message:

- GSE167339 supports donor-level transfer robustness in an independent
  perturbation/scar-system fibroblast dataset.
- GSE156326 and GSE181316 recover ECM-rich scar-context composition shifts and
  directionally positive but non-FDR-significant disease-control module-score
  deltas after transfer.
- Bootstrap and VI uncertainty summaries are moved to uncertainty/supplementary
  material.

Claim boundary:

- Figure 5 supports external reproducibility and biological utility, not
  diagnostic performance, causal mechanism, lineage tracing, or therapeutic
  efficacy.

### Result 6. Transition Flow and FPI Summarize Plasticity

Primary figure: Figure 6

Core message:

- State-flow estimation uses optimal transport over fitted state compositions.
- FPI combines mixed-state entropy and transition potential.

Claim boundary:

- Transition flow is cross-sectional and should not be interpreted as lineage
  tracing or observed temporal conversion.

### Result 7. Cross-Cohort Transfer Provides a Migration Diagnostic

Primary outputs:

- `analysis/cross_cohort_transfer/`
- `analysis/public_realdata_transfer/`
- `analysis/multi_public_realdata_validation/`
- `analysis/independent_geo_gse246215_validation/`

Core message:

- A fitted NB state program can be frozen and transferred to held-out cells.
- Simulation transfer reports recovery metrics and z convergence.
- Public real-data transfer smoke runs bidirectionally between the two public
  conditions.
- Multi-public validation runs a pooled public raw-count fit and
  leave-dataset-out transfer across registry-defined public count matrices.
- Independent GSE246215 validation runs the same workflow on human fibroblast
  cancer-type subsets from GC, HCC, NSCLC, and TNBC.
- GSE246215 sensitivity repeats validation across random downsampling seeds and
  q95 library-size trimming, with transfer convergence and composition
  variability summaries.

Claim boundary:

- This supports transfer mechanics, not full real-human multi-cohort atlas
  generalization.
- The current default multi-public registry uses two condition-specific matrices
  from one Dryad record; a completed cross-study human disease-atlas claim
  remains future work.
- GSE246215 adds an independent human source but remains a processed
  count-matrix and within-study cancer-type validation.
- HCC remains a library-size outlier in the sampled public matrix; HCC-specific
  state composition should be treated as QC-sensitive and hypothesis-generating.

## Methods Outline

### Data and Preprocessing

- Raw gene-by-cell UMI matrices are required.
- Metadata alignment and marker-prior retention are handled by
  `prepare_fibrodynmix_data()`.
- Public smoke analyses use DOI 10.6071/M3238R with explicit claim boundaries.

### Generative Model

Describe the NB mean model:

```text
x_ig ~ NB(mu_ig, phi_g)
log(mu_ig) = log(l_i) + alpha_g + z_i beta_g + b_sg + u_jg
```

Describe the latent simplex:

```text
eta_i ~ logistic-normal
z_i = softmax(eta_i)
```

### Optimization

- Initialize state weights and programs with `fit_fibrodynmix_initializer()`.
- Optimize NB model with alternating gene-level and cell-level updates.
- Fit optional study/donor effects with ridge penalties.
- Track objective traces, rollback, early stopping, and effect norms.

### Posterior and Uncertainty

- Bootstrap uncertainty uses cell resampling.
- VI posterior uses diagonal logistic-normal state-logit approximation around
  the fitted NB mode.
- Simulation interval calibration uses known latent truth and is not available
  for real-data claims.

### Transition Flow and FPI

- Compute state-program transcriptional cost from fitted beta.
- Estimate entropy-regularized optimal transport between condition-level state
  compositions.
- Compute FPI from entropy plus transition potential.

### Cross-Cohort Transfer

- Freeze alpha, beta, and phi from training fit.
- Optimize held-out z under the raw-count NB likelihood.
- Report held-out objective, convergence rate, and composition summaries.

### Software and Reproducibility

- R package version: 0.0.0.9000.
- Quality gates:
  - `testthat::test_local()`
  - `scripts/check_project_integrity.R`
  - `R CMD check --no-manual`
- Project maturity audit:
  - `scripts/run_project_maturity_audit.R`

## Data and Code Availability Draft

All code required to reproduce the analyses is organized in the FibroDynMix R
package. Public raw-count smoke analyses use the Hum/Sebastian breast fibroblast
record, DOI 10.6071/M3238R. Generated analysis outputs are stored under
`analysis/`, and manuscript figure packages with source data, manifests, legends,
and export QC are stored under `figures/`. Large public raw-count files and
generated analysis/figure directories are excluded from the R source package
tarball but remain part of the local reproducibility workspace.

## Limitations Draft

FibroDynMix currently implements a pragmatic raw-count NB optimizer and a
lightweight logistic-normal VI posterior over cell state logits. It does not yet
implement full amortized inference or a joint variational posterior over all
model parameters, study effects, donor effects, and dispersions. Public
real-data analyses are bounded count-matrix validations and case studies rather
than diagnostic, causal, therapeutic, or completed human disease-atlas evidence.
Transition flow is inferred from cross-sectional state composition and
transcriptional cost, not from lineage tracing or observed temporal
trajectories.
