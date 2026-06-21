# FibroDynMix Manuscript Draft

## Title Options

1. FibroDynMix models fibroblast state plasticity as raw-count latent mixtures
2. FibroDynMix: a negative-binomial generative framework for fibroblast state mixtures, uncertainty, and transfer
3. Raw-count generative modeling of fibroblast state plasticity with FibroDynMix
4. A hierarchical mixture model for fibroblast state plasticity from single-cell counts

Recommended title: **FibroDynMix: a negative-binomial generative framework for
fibroblast state mixtures, uncertainty, and transfer**.

## Abstract

Fibroblast phenotypes are often analyzed with fixed marker scores, but such
scores can blur continuous state mixtures, transitional phenotypes, and
cohort-level expression shifts. We introduce FibroDynMix, a raw-count
negative-binomial generative framework that represents fibroblast states as
latent compositional programs on a simplex. The current implementation combines
weak marker-prior orientation with count-likelihood optimization, ridge
study/donor effects, bootstrap and logistic-normal posterior state uncertainty,
transition-flow summaries, and cross-cohort transfer diagnostics. In simulations
with known latent truth, FibroDynMix quantifies state-weight recovery,
dominant-state assignment, mixed-state entropy preservation, rare-transition
detection, learned marker-program recovery, and posterior interval calibration.
The package also executes end-to-end on public mouse breast fibroblast raw UMI
counts, producing auditable negative-binomial diagnostics, state composition,
plasticity, transition-flow, bidirectional transfer smoke-test outputs, and
registry-driven multi-public leave-dataset-out transfer diagnostics. An
independent public human fibroblast count-matrix validation using GSE246215
extends this test to four cancer-type fibroblast subsets, while GSE167339,
GSE156326, and GSE181316 provide bounded external biological-utility checks.
FibroDynMix is distributed as a reproducible R package with documented functions,
source-data-backed figures, release metadata, and project integrity gates. These
results support FibroDynMix as a method-development framework for modeling
fibroblast plasticity beyond direct gene-set scoring, while leaving full
all-parameter Bayesian inference and completed human disease-atlas or
cross-atlas generalization as future extensions.

## Introduction

Fibroblasts occupy a broad spectrum of tissue-resident, inflammatory,
matrix-remodeling, myofibroblast, antigen-presenting, and stress-associated
phenotypes. Single-cell RNA-seq has made these phenotypes observable at scale,
but the computational representation of fibroblast state plasticity remains a
limiting step. In many analyses, cells are assigned to discrete labels or scored
against fixed marker sets. These approaches are easy to interpret, yet they are
not generative models of raw counts and do not naturally represent cells that
lie between states.

This limitation matters because fibroblast plasticity is often expressed as a
continuum rather than a set of mutually exclusive categories. A cell can carry
resident and inflammatory features simultaneously, or occupy a transitional
state between inflammatory, myofibroblast, and extracellular-matrix programs.
Marker scores can detect enrichment of predefined genes, but they do not model
the count-generating process, quantify uncertainty in the state mixture, or
separate biological state variation from study- and donor-level expression
shifts.

Several probabilistic and latent-factor strategies address parts of this
problem, including topic models, matrix factorization, and deep latent-variable
models. However, a fibroblast-focused method-development framework also needs to
make biological state programs identifiable, retain raw-count likelihoods,
account for cross-cohort effects, estimate transition-like structure, report
uncertainty, and support transfer to held-out cohorts. These requirements are
especially important for methods intended to support disease and tissue
comparison across public single-cell studies.

Here we introduce FibroDynMix, a hierarchical raw-count mixture framework for
fibroblast state plasticity. FibroDynMix models each cell as a latent mixture of
fibroblast state programs under a negative-binomial count likelihood. Marker
sets are used as weak orientation priors rather than as direct scores. The
implemented R package includes simulation, negative-binomial optimization,
study/donor effect modeling, bootstrap uncertainty, a lightweight
logistic-normal variational posterior over cell state logits, transition-flow
estimation, cross-cohort transfer, and public raw-count smoke analyses. We
evaluate the method in simulations with known latent truth and organize all
figures, source data, package metadata, and claim boundaries as a reproducible
method-development workflow.

## Results

### FibroDynMix Defines a Raw-Count Generative Framework for Fibroblast State Mixtures

FibroDynMix was designed to distinguish latent state inference from direct
gene-set scoring. The model starts from a raw gene-by-cell UMI count matrix and
uses a negative-binomial likelihood with a library-size offset. Each cell is
represented by a state-mixture vector on the simplex, and state-specific gene
programs contribute to the log mean of the count model. Study- and donor-level
gene effects can be added on the same linear predictor scale. Figure 1
summarizes this architecture and maps each implemented layer to the current
evidence base.

The current R package exposes the likelihood, deviance, and optimization
objective through `fibrodynmix_nb_loglik()`, `fibrodynmix_nb_deviance()`, and
`fibrodynmix_nb_objective()`. The fitting workflow initializes latent state
weights and state-gene programs, then alternates negative-binomial updates for
gene-level coefficients and cell-level simplex weights. This implementation
supports raw-count optimization and explicit hierarchy effects, but it remains a
pragmatic optimizer rather than a full joint Bayesian inference engine.

### Public Raw-Count Smoke Analyses Demonstrate End-to-End Execution

We next evaluated whether the workflow could run on public raw UMI count data
rather than only on simulated matrices. Using public mouse breast fibroblast and
CAF count files from DOI 10.6071/M3238R, the public smoke workflow prepared raw
counts, retained weak-prior markers, fit the negative-binomial model, estimated
state composition, computed Fibroblast Plasticity Index (FPI), and inferred a
cross-sectional normal-to-disease state-flow summary. Figure 2 and the public
real-data smoke figure provide source-data-backed outputs for this technical
demonstration.

These analyses support real raw-count execution and auditability. They do not
support disease-mechanism claims because the public input is a pooled mouse
smoke-test dataset rather than a replicated human disease-atlas validation.
This boundary is recorded in the figure manifests and public-data
documentation.

### Simulations Quantify State Recovery, Entropy Preservation, and VI Behavior

To evaluate method behavior under known truth, we built four simulation
scenarios: continuous state mixtures, discrete states, batch confounding, and
rare transition-like cells. The benchmark compares marker scoring, an
NMF/topic baseline, the FibroDynMix initializer, the negative-binomial
optimizer, study-effect NB, and the FibroDynMix VI posterior mean. It reports
state-weight RMSE, dominant-state assignment, entropy preservation,
rare-transition detection, recovery of true state-marker programs, and
donor-grouped downstream prediction of the simulated disease label from each
method's inferred feature representation.

The simulation results provide the main method-validation evidence because the
true latent state mixture is known. They also expose the limits of the current
posterior layer. The extended method benchmark shows that ordinary NMF/topic
factorization is weaker than the FibroDynMix-family methods in the current
state-recovery simulations. Marker scoring remains competitive in several
marker-aligned settings, so the manuscript should not claim uniform superiority
over scoring. The updated downstream task benchmark shows that FibroDynMix-family
features also improve mean simulated disease-label balanced accuracy and
macro-F1 over the marker/NMF baseline average, although marker scoring remains
strong in clean-prior scenarios. A new NB convergence benchmark shows that the
current bounded `n_outer = 2` setting should be treated as a smoke/runtime
configuration rather than a convergence-adequate optimizer setting: objective
values and inferred `z` continue to change at larger outer-iteration budgets.
The stronger claim is therefore limited to the implemented framework and
bounded benchmark behavior until primary NB results are rerun with explicit
convergence or model-selection criteria. FibroDynMix estimates latent mixtures
under a raw-count likelihood and adds hierarchy, uncertainty, transfer, and
transition diagnostics that scoring and ordinary NMF/topic baselines do not
provide. A higher-outer core rerun with `n_outer = 10` and `20` confirms the
issue: lower NB training objectives consistently worsen `z` RMSE in the current
bounded benchmark, so objective minimization alone is not an adequate
state-recovery selection rule. The lightweight VI posterior reports raw interval coverage and
simulation-calibrated interval coverage; in the current benchmark, calibration
raises mean interval coverage from 0.276 to 0.924 at the cost of wider intervals.
This supports uncertainty calibration as a measured diagnostic rather than a
qualitative claim, while making clear that the calibration uses simulated truth
and is not directly available in real data.

### Study-Effect Calibration Addresses Cohort Confounding

Cross-cohort fibroblast analyses must distinguish state variation from study or
batch effects. FibroDynMix therefore includes ridge-penalized study and donor
gene effects. Figure 4 evaluates the study-effect penalty under simulated batch
confounding, measuring state-weight RMSE, negative-binomial objective, fitted
study-effect magnitude, and a weighted tradeoff score. The selected default
penalty balances state recovery with shrinkage of cohort-level expression
effects.

This calibration supports the use of hierarchy-aware effects in the current
optimizer. It does not imply that the selected penalty is universally optimal.
For new cohorts, the same penalty-sensitivity workflow should be rerun, and the
current ridge estimates should be interpreted as regularized effects rather than
full posterior random effects.

### External Human Fibroblast Datasets Validate Biological Utility

To move beyond a single cancer-atlas case study, we promoted three independent
human fibroblast datasets into the main biological validation layer. Figure 5
combines donor-level perturbation evidence from GSE167339 with external
scar-context transfer and module validation from GSE156326 and GSE181316.
GSE167339 supports donor-level transfer robustness, while the scar cohorts test
whether transferred FibroDynMix states recover ECM-rich disease contexts and
directionally positive scar-module shifts.

This evidence supports external reproducibility and disease-context biological
utility. The scar-module contrasts are descriptive small-sample pseudobulk
effects and do not pass FDR significance, so they should not be interpreted as
diagnostic performance, causal mechanism, lineage tracing, or therapeutic
efficacy. The bootstrap uncertainty
and lightweight VI posterior summaries remain available as supplementary
uncertainty evidence: they report state-weight, entropy, sample-composition, and
marker-program intervals, but they do not constitute a full joint posterior over
all model parameters.

### Transition Flow and FPI Summarize Plasticity

To summarize plasticity from fitted state mixtures, FibroDynMix computes a
state-program cost matrix and estimates an entropy-regularized optimal-transport
flow between condition-level state compositions. Figure 6 uses this flow to
derive a cell-level Fibroblast Plasticity Index that combines mixed-state
entropy with transition potential.

The transition-flow layer is intended as a cross-sectional summary of
composition shifts and transcriptional cost. It should not be interpreted as
lineage tracing, temporal observation, or proof of a physical state transition.
This distinction is important for maintaining a defensible biological claim.

### Cross-Cohort Transfer Provides a Diagnostic for Program Migration

FibroDynMix can freeze a fitted state program and optimize held-out cell state
weights under the raw-count negative-binomial likelihood. In leave-study-out
simulation benchmarks, the transfer workflow reports held-out state recovery,
negative-binomial objective, entropy, and convergence diagnostics. In public
real-data smoke analyses, the fitted program was transferred bidirectionally
between the public normal and disease count files, with convergence rates of
1.000 and 0.99375 in the two directions.

These results support the transfer mechanics and diagnostics of the current
implementation. They do not yet establish completed human cross-atlas
generalization, which would require harmonized cell-type curation, independent
cohort replication, and stronger biological validation.

## Discussion

FibroDynMix reframes fibroblast state analysis as raw-count latent mixture
modeling rather than fixed marker scoring. This distinction is central to the
method: state weights are inferred as latent variables under a count likelihood,
while marker sets orient state programs without directly defining state scores.
The current package combines this modeling strategy with hierarchy-aware effects,
uncertainty summaries, transition-flow analysis, and transfer diagnostics.

The simulation benchmark is the strongest current evidence because it evaluates
state mixtures against known latent truth. These simulations show how the method
behaves under continuous mixtures, discrete states, batch confounding, and rare
transition-like cells. The VI calibration results are especially useful because
they show both the raw interval behavior and the width cost required to reach
target coverage under simulated truth. This helps avoid an overconfident
uncertainty claim. The NB convergence benchmark also shows that current
`n_outer = 2` NB fits are not sufficient for convergence-level claims; existing
NB, VI, transfer, and study-effect results should be presented as bounded
prototype evidence unless rerun with a convergence grid or a justified
early-stopping/model-selection rule.

The public raw-count analyses provide a complementary engineering check. They
show that the workflow can ingest real raw count matrices and produce auditable
outputs. The multi-public validation layer further fits a pooled model across
condition-specific public count matrices with study effects and then performs
leave-dataset-out transfer diagnostics. We also prepared GSE246215, an
independent public human fibroblast atlas count matrix, into four cancer-type
datasets and ran the same pooled-fit and leave-dataset-out transfer workflow.
We further stress-tested GSE246215 with random downsampling and q95
library-size trimming. Across four sensitivity runs, the minimum transfer
convergence rate was 0.9833 and the maximum state-composition standard deviation
across seed repeats was 0.0420. However, HCC remained a major library-size
outlier after simple trimming, so HCC-specific composition should be treated as
QC-sensitive. We also added a GSE246215 downstream benchmark that predicts
cancer type from patient- or sample-level fibroblast representations. In the
current bounded run, NMF topics are the strongest cancer-type classifier,
FibroDynMix `z` is competitive with marker-derived representations but does not
dominate them, and canonical myofibroblast, ECM, and inflammatory marker
gradients correlate positively with their expected FibroDynMix state weights.
More broadly, the GSE246215 layer uses processed public count matrices and
sampled cancer-type subsets from one GEO study. Without replicated independent
human studies, FASTQ-level reprocessing, and full biological validation, these
outputs should not be interpreted as disease mechanisms or atlas-level findings.

The external validation layer addresses this limitation in a bounded way by
adding GSE167339, GSE156326, and GSE181316 as independent human fibroblast
datasets. These analyses support the method-first manuscript with external
biological-utility evidence by testing conserved and context-specific
fibroblast state plasticity in non-cancer scar and perturbation settings. The
claim remains deliberately bounded because the source data are public processed
matrices and the GSE156326/GSE181316 module statistics are small-sample
pseudobulk validations with non-significant adjusted q values.

Several limitations remain. The negative-binomial optimizer is not yet a full
joint Bayesian inference procedure. The VI posterior is restricted to cell state
logits around the fitted mode and does not yet cover all parameters. The
transition-flow layer is cross-sectional rather than temporal. The transfer
analysis demonstrates mechanics and smoke-test feasibility rather than definitive
cross-atlas generalization. These limitations define the next development
steps: full amortized or fully hierarchical VI, stronger harmonized human
multi-cohort validation, and broader comparisons against deep generative
single-cell baselines.

Overall, FibroDynMix provides a reproducible method-development framework for
probabilistic fibroblast state analysis. Its current contribution is not a final
fibroblast disease atlas, but a mature R package and evidence-backed modeling
workflow for moving beyond direct gene-set scoring toward raw-count state
mixture inference.

## Methods

### Software Implementation

FibroDynMix is implemented as an R package. The package exposes functions for
simulation, raw-count likelihood evaluation, negative-binomial optimization,
real-data preparation, bootstrap uncertainty, VI posterior summaries,
transition-flow estimation, cross-cohort transfer, and benchmarking. The package
includes documentation for all exported functions, release notes, citation
metadata, a function index, and project integrity checks.

### Simulation Model

The simulator generates gene-by-cell raw counts under a negative-binomial model.
For cell `i` and gene `g`, the mean is modeled as a library-size offset plus a
gene baseline, a latent state-mixture contribution, and optional study/donor
effects:

```text
x_ig ~ NB(mu_ig, phi_g)
log(mu_ig) = log(l_i) + alpha_g + z_i beta_g + b_sg + u_jg
```

The latent state vector `z_i` lies on the simplex and is generated by a
logistic-normal layer. Marker genes are simulated as weakly oriented
state-specific programs, allowing marker recovery to be evaluated against known
truth.

### Negative-Binomial Optimization

The optimizer initializes state mixtures and state-gene programs, then performs
alternating updates of gene-level coefficients and cell-level state logits under
the negative-binomial likelihood. Optional study and donor effects are fit with
ridge penalties. The optimizer records initial, final, and best objectives,
best iteration, rollback flags, stopping reason, and effect norms.

NB optimizer convergence is evaluated separately by running outer-iteration
grids and comparing objective traces plus inferred state weights. Current
bounded results show that `n_outer = 2` is not convergence adequate for primary
NB claims.

### Posterior and Uncertainty

Bootstrap uncertainty is estimated by resampling cells and refitting the selected
model. The VI posterior layer fixes the fitted negative-binomial program and
places a diagonal logistic-normal approximation over cell state logits. Posterior
draws are transformed through softmax to state weights, from which state-weight,
entropy, and sample-composition intervals are computed. In simulation,
posterior interval coverage is evaluated against true latent state weights, and
post-hoc interval-scale calibration is reported as a diagnostic.

### Transition Flow and FPI

State-to-state transcriptional cost is computed from fitted state programs.
Condition-level state compositions are connected by entropy-regularized optimal
transport. FPI is computed from normalized state-mixture entropy plus
transition-potential contribution from the inferred flow.

### Cross-Cohort Transfer

Transfer fixes the fitted training program parameters and optimizes held-out
cell state weights under the raw-count likelihood. The workflow reports held-out
negative-binomial objective, log-likelihood, shared genes, convergence status,
cell-level convergence rate, state composition, and entropy.

### Public Data and Multi-Public Validation

Public raw-count analyses use the Hum/Sebastian breast fibroblast count record
with DOI 10.6071/M3238R. The first analyses are run as technical smoke tests
for raw-count execution and bidirectional transfer. A second registry-driven
analysis fits a pooled model across multiple public count matrices, reports
dataset-level composition, estimates a condition-level transition flow, and
performs leave-dataset-out transfer. We additionally prepare the GSE246215
public human fibroblast atlas count matrix into GC, HCC, NSCLC, and TNBC
fibroblast subsets and run the same leave-dataset-out validation workflow.
GSE246215 sensitivity analyses repeat this workflow under random downsampling
and q95 library-size trimming, reporting transfer convergence and
state-composition variability. A downstream GSE246215 benchmark compares
FibroDynMix state composition, marker-scoring features, and NMF topics for
patient- or sample-level cancer-type classification, and records canonical
marker-gradient correlations against expected FibroDynMix states.
These outputs are not treated as mechanistic disease-discovery analyses or
completed human multi-study atlas validation.

### Reproducibility and Quality Gates

The project includes source-data-backed figures, panel manifests, figure
legends, export QC, release metadata, and a project maturity audit. The standard
validation gate is:

```bash
Rscript -e 'testthat::test_local()'
Rscript scripts/check_project_integrity.R
R CMD build .
R CMD check --no-manual FibroDynMix_0.0.0.9000.tar.gz
```

## Data and Code Availability

The FibroDynMix R package contains all code required to reproduce the simulation
benchmarks, model fits, posterior diagnostics, transition-flow analyses,
cross-cohort transfer analyses, and figure packages. Public raw-count smoke
analyses use the Hum/Sebastian breast fibroblast record, DOI 10.6071/M3238R.
Generated analysis outputs are stored under `analysis/`, and figure source data,
manifests, legends, and export QC are stored under `figures/`. Large public
raw-count files and generated analysis directories are excluded from the R source
package tarball through `.Rbuildignore`.

## Limitations

FibroDynMix is currently a method-development package. It implements a
raw-count negative-binomial optimizer, bootstrap uncertainty, and a lightweight
VI posterior over state logits, but it is not yet a full all-parameter Bayesian
inference engine. Public real-data analyses are bounded count-matrix
validations and case studies rather than diagnostic, causal, therapeutic, or
completed human disease-atlas evidence. Transition flow is cross-sectional
optimal transport and should not be interpreted as lineage tracing. These
limitations should remain explicit in the title, abstract, results, discussion,
and figure legends.
