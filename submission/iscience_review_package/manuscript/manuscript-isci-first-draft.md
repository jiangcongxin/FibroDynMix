# FibroDynMix models fibroblast state plasticity as raw-count latent mixtures

## Summary

Single-cell fibroblast studies often reduce continuous, plastic states to fixed
marker scores that cannot jointly model count noise, cohort effects, and
uncertainty. We present FibroDynMix, an R package that infers simplex fibroblast
state weights from raw UMI counts under a negative-binomial likelihood, coupling
weak marker orientation, study/donor terms, uncertainty, transfer, and
transition flow in one auditable workflow. Using simulations with known truth, we
find that a lower training count objective can yield worse state recovery,
motivating validation-aware selection that scores held-out fit, state stability,
marker-gradient preservation, and downstream utility rather than training
objective alone. In a benchmark against marker and topic baselines,
FibroDynMix-family features improved state-weight recovery (RMSE 0.2296 versus
0.2568) and downstream classification, and validation-aware defaults further
improved recovery and transfer. It generalized across public mouse and human
fibroblast matrices through transfer, recovered marker gradients, and
directional scar-module signals. FibroDynMix makes fibroblast state inference a
reproducible, uncertainty-aware count model.

## Highlights

- FibroDynMix infers fibroblast states from raw UMI counts, not marker scores
- A single workflow couples count likelihood, uncertainty, and transfer
- A lower training objective can worsen latent-state recovery
- Validation-aware selection and transfer generalize to public fibroblast data

## Graphical Abstract

The graphical abstract draft is provided as a separate 1200 x 1200 px
single-panel artwork file. Raw UMI count matrices enter FibroDynMix, which fits
a negative-binomial mixture model with weak marker orientation, simplex state
weights, and study/donor terms. Fitted state mixtures then feed uncertainty,
transfer, cross-sectional transition-flow, and FPI analyses, with simulation
truth and public fibroblast datasets defining the claim boundary.

## Introduction

Fibroblasts occupy resident, inflammatory, matrix-remodeling,
myofibroblast-like, antigen-presenting, and stress-associated states across
healthy and diseased tissues. Cross-tissue and disease-focused single-cell
atlases have identified fibroblast activation states in cancer, chronic
inflammation, fibrosis, and wound repair [1-7]. Many analyses still represent
fibroblast identity with discrete labels or fixed marker scores [8,9]. These
summaries are useful for annotation, but they compress continuous mixtures and
do not model the raw count-generating process.

Mixed fibroblast states are central to plasticity. A cell can carry resident and
inflammatory features at the same time, or sit between inflammatory,
ECM-remodeling, and myofibroblast programs. Matrix factorization studies have
shown that single-cell profiles can contain mixtures of identity and activity
gene programs [10-12]. Marker scores can detect enrichment of predefined genes,
but they do not estimate latent state proportions, account for count noise, or
separate state variation from study- and donor-level expression shifts. Deep
latent models and negative-binomial count models address parts of this problem,
yet many single-cell fibroblast studies still need a model that keeps state
programs biologically oriented while preserving raw-count likelihoods and
transparent diagnostics [13-16].

FibroDynMix addresses this gap by modeling each cell as a composition of
fibroblast state programs under a negative-binomial count model. Marker genes
orient the state programs but do not directly define the final state weights.
Single-cell integration methods commonly correct batch effects by aligning cells
in expression or embedding spaces [17,18]; FibroDynMix instead fits optional
study and donor terms on the same log-mean scale as the state programs. The
package includes simulation, model fitting, study and donor effects, bootstrap
uncertainty, a lightweight logistic-normal posterior over state logits,
transition-flow estimation, and transfer to held-out cells or datasets. These
components give fibroblast studies a raw-count model whose outputs can be
audited against marker scoring, topic models, simulation truth, and public count
matrices.

We evaluated FibroDynMix in four simulation scenarios with known latent state
truth, public mouse fibroblast raw UMI counts, the GSE246215 human fibroblast
count matrix, the GSE167339 perturbation/scar-system fibroblast dataset, and two
external scar/fibrosis-context datasets, GSE156326 and GSE181316 [3,19-23].
Because count-likelihood improvement did not always track latent-state recovery,
the final workflow uses validation-aware NB selection rather than objective-only
stopping. Each main figure is linked to source data, claim boundaries, and
quality gates. This evidence sequence evaluates FibroDynMix as a method-first
package with bounded external biological validation.

## Results

### FibroDynMix estimates simplex fibroblast states from raw UMI counts

FibroDynMix starts from a gene-by-cell UMI count matrix and represents each cell
with a state-weight vector on the simplex. The count mean includes a library-size
offset, a gene baseline, state-program contributions, and optional study and
donor effects. Figure 1 summarizes the model and the evidence linked to each
implemented component.

The package exposes the negative-binomial likelihood and objective through
documented R functions and fits the model by alternating updates of gene-level
coefficients and cell-level state logits. Marker genes enter as weak orientation
terms. This design makes the fitted state weights different from direct
signature scores: the state mixture is inferred under a raw-count likelihood and
can be evaluated against likelihood diagnostics, marker-program recovery, and
held-out transfer behavior.

The same fit returns optimizer diagnostics for the biologically oriented mixture
model. It includes study and donor effects, bootstrap intervals, a
logistic-normal posterior approximation around the fitted state logits, and
transition-flow summaries. Cell-bootstrap uncertainty summaries are shown in
Supplementary Figure S1. These components make FibroDynMix an auditable R
package rather than a black-box annotation workflow.

### Public raw-count matrices test package execution and bidirectional transfer

We first tested FibroDynMix on public raw UMI count matrices. The public mouse
breast fibroblast and CAF count files from DOI 10.6071/M3238R were processed as
a technical real-data test [19,20]. The Figure 2 workflow retained a balanced
subset of 600 cells and 1,217 genes, fit the negative-binomial model, estimated
state composition, computed FPI, and inferred a cross-sectional
normal-to-disease state-flow summary. The NB objective improved by 0.0900 and
the transition-flow calculation converged.

Supplementary Figure S2 repeated the workflow on 240 cells and 720 genes and
produced an NB objective improvement of 0.1227. A bidirectional transfer
analysis then froze the fitted state program in one condition and optimized
held-out state weights in the other condition. Transfer convergence was high in
both directions in the current run. These analyses show that the package can
ingest public raw counts, fit the model, and export source-backed diagnostics
across independent count matrices.

### Simulation truth separates state recovery from marker and topic baselines

Simulation provides the primary test of latent state recovery because the true
state mixture is known. Figure 3 visualizes the core simulation readouts for
state-weight recovery, entropy preservation, rare transition-like cells, and
marker-program recovery. The extended benchmark added NMF/topic modeling,
study-effect NB, and the FibroDynMix VI posterior mean [10-13].

Across 48 rows in the bounded extended benchmark, FibroDynMix-family methods
had lower mean state-weight RMSE than the marker/NMF baseline average (0.2296
versus 0.2568). Mean dominant-state accuracy was also higher for the
FibroDynMix-family methods (0.3815 versus 0.3086). In the simulated downstream
disease-label task, FibroDynMix-family features had higher mean balanced
accuracy (0.7383 versus 0.7090) and macro-F1 (0.7371 versus 0.7049).

Together, the simulations show that FibroDynMix can recover latent-state
structure while adding count likelihoods, study/donor effects, uncertainty,
transfer, and transition summaries. Marker scoring remained strong in
marker-aligned simulations, which made the benchmark useful as a calibration
step rather than a simple method-ranking exercise. This motivated the
validation-aware selection layer tested next: the final FibroDynMix workflow
requires fitted count models to preserve held-out fit, state stability,
marker-gradient structure, and downstream utility instead of relying on the NB
training objective alone.

### Validation-aware selection separates NB fit from state recovery

Initial NB benchmarks used fast smoke settings, so we reran core methods under
higher outer-iteration budgets. The core-converged benchmark tested
`n_outer = 10` and `20` across four scenarios. In the 80-row run, marker scoring
and the initializer remained stronger than full NB in RMSE for most scenarios.
Mean NB-family RMSE was 0.2622 compared with 0.2514 for the marker/NMF baseline
average, and mean downstream balanced accuracy was 0.6338 compared with 0.6836.

The core-converged runs exposed a model-selection problem (Supplementary Figure S3).
All 16 paired `n_outer = 10` to `20` contrasts had lower NB objective and
worse latent-state RMSE. In the five-replicate scVI-inclusive rerun,
Supplementary Figure S4 showed the same pattern in 40 of 40 objective/RMSE
tradeoff pairs. A lower training count objective therefore did not guarantee
better recovery of the simulated latent state.

FibroDynMix therefore uses a validation-aware selector that combines held-out NB
objective, state-weight stability, marker-gradient preservation, and downstream
validation. Selected defaults improved the final VI and transfer summaries:
mean selected-NB VI RMSE was 0.202, mean selected-NB VI downstream balanced
accuracy was 0.766, mean selected-NB transfer RMSE was 0.215, and transfer
state-weight convergence was 1.000. These results support explicit validation
criteria for choosing NB fits instead of objective-only stopping.

### Study-effect terms reduce simulated cohort-shift sensitivity

Study and donor effects are fitted on the same log-mean scale as the state
programs. Figure 4 evaluates ridge penalties under simulated batch confounding
and tracks state-weight RMSE, NB objective, study-effect magnitude, and marker
orientation. The selected penalty pair was study L2 = 5 and marker L2 = 0.05.

This calibration shows that FibroDynMix can absorb cohort-level expression
shifts while retaining state-mixture recovery in simulation. The study-effect
term tests whether cohort terms alter inferred fibroblast states. New datasets
require the same penalty grid rather than an assumed penalty.

### Human fibroblast datasets support transfer and bounded biological utility

The GSE246215 human fibroblast count matrix was split into GC, HCC, NSCLC, and
TNBC fibroblast subsets and analyzed with the same leave-dataset-out transfer
workflow [3]. Sensitivity runs repeated the workflow under random downsampling and
q95 library-size trimming. Across four sensitivity runs, the minimum transfer
convergence rate was 0.9833 and the maximum state-composition standard
deviation across seed repeats was 0.0420. HCC remained a library-size outlier
after simple trimming, so HCC-specific composition is treated as
QC-sensitive.

GSE246215 also tested whether FibroDynMix state weights could serve as
downstream features. Patient-level cancer-type classification was strongest for
NMF topics, with balanced accuracy of 0.8100. FibroDynMix state weights reached
balanced accuracy of 0.6314, macro-F1 of 0.7585, and macro-AUROC of 0.8141,
exceeding raw marker-score features on macro-F1 and macro-AUROC. Marker-gradient
checks showed that expected fibroblast programs aligned with the corresponding
state weights: ECM-remodeling rho = 0.6714, myofibroblast rho = 0.4236, and
inflammatory rho = 0.3475.

Figure 5 extends the real-data layer to independent public human fibroblast
datasets. In GSE167339, leave-donor-out transfer convergence remained high
(Human1 = 1.000, Human2 = 1.000, Human3 = 0.978). Human3 was used only as
hash-unknown donor robustness evidence because the public files do not map
HumanHashTag IDs to treatment labels. External scar/fibrosis-context transfer
also converged in GSE156326 and GSE181316 (0.954 and 0.960) [21-23].
Scar-module deltas were positive in both cohorts, with GSE156326 delta = 0.028,
AUC = 0.78, q = 0.38 and GSE181316 delta = 0.085, AUC = 1.00, q = 0.061. These
external scar analyses were sample-limited and did not pass FDR significance, so
the module deltas are interpreted as directional effect-size evidence rather
than confirmed differential-state signals. These public datasets support
reproducibility and bounded biological utility across independent fibroblast
count matrices.

### Cross-sectional transition flow links state mixtures to FPI

FibroDynMix uses fitted state programs to compute a state-to-state
transcriptional cost matrix. Condition-level state compositions are then linked
by entropy-regularized optimal transport, with mathematical support from
optimal-transport algorithms and prior single-cell use in developmental
time-course analysis [24-26]. Figure 6 uses this state flow to derive a
cell-level Fibroblast Plasticity Index that combines state-mixture entropy with
transition potential.

The Figure 6 analysis converged and produced a flow entropy of 2.1602 with an
expected cost of 15.3812. The resulting FPI summarizes which cells carry mixed
state structure and lie near inferred condition-associated flow. We treat FPI
as a cross-sectional state-composition summary linked to the fitted model's
transcriptional cost structure.

## Discussion

FibroDynMix estimates fibroblast state mixtures directly from raw UMI counts
rather than from normalized marker scores alone. Across simulations and public
count matrices, the model returned simplex state weights, study/donor-adjusted
count fits, uncertainty summaries, transfer diagnostics, and transition-flow
outputs from one R workflow. This extends marker scoring and matrix
factorization approaches [8-12] by keeping marker orientation weak and auditable
while fitting a count likelihood.

Simulation truth showed where the method gains and where simpler baselines
remain strong. FibroDynMix-family features improved mean state-weight RMSE and
downstream classification in the extended benchmark, but marker scoring remained
competitive in clean marker-aligned simulations. The higher-outer NB reruns also
showed that lower training objectives could worsen latent-state RMSE. These
results support validation-aware selection because the best count fit is not
necessarily the best recovery of the biological state quantity being inferred.

Public fibroblast datasets tested execution, transfer, and bounded biological
utility rather than broad atlas-level discovery. The mouse breast fibroblast
matrices supported end-to-end raw-count fitting and bidirectional transfer.
GSE246215 recovered expected fibroblast marker gradients, although NMF topics
were the strongest cancer-type classifier and HCC composition remained
QC-sensitive. GSE167339 supported donor-level transfer robustness, and
GSE156326/GSE181316 showed directional scar-module shifts that were consistent
with extracellular-matrix remodeling in keloid scars [27] but not
FDR-significant. The human-data claim is therefore limited to transfer,
marker-gradient recovery, and directional scar-context evidence.

The package-level audit trail is part of the method claim. FibroDynMix includes
documented functions, source-data-backed figures, figure manifests, release
metadata, runtime records, project maturity audits, and package checks. This
structure makes each figure traceable to scripts and source tables, which is
essential for a computational method whose claims depend on reproducible model
fits rather than a single static annotation.

## Limitations of the study

FibroDynMix currently fits a pragmatic negative-binomial optimizer, not a full
joint Bayesian model over every parameter. The VI layer approximates posterior
uncertainty over cell state logits around the fitted mode and does not include a
full posterior over gene programs, study effects, donor effects, or dispersion.
The study-effect and donor-effect terms are penalized fixed-effect corrections
rather than full random-effect posteriors. The NB objective can improve while
latent-state RMSE worsens, so final model selection requires validation-aware
criteria. Public real-data analyses use processed or public count matrices and
bounded sampling. They do not establish diagnostic performance, causal disease
mechanisms, therapeutic relevance, or a completed human disease atlas.
Transition flow and FPI summarize cross-sectional state composition and
transcriptional cost; they do not measure observed temporal conversion or
lineage relationships.

## STAR Methods

### Key resources table

| Reagent or resource | Source | Identifier |
|---|---|---|
| Mouse breast fibroblast and CAF raw count matrices | Sebastian et al. [19]; Hum et al. [20] | DOI: 10.6071/M3238R |
| GSE246215 human fibroblast count matrix | Gao et al. [3] | GEO: GSE246215 |
| GSE167339 collagen scar-system fibroblast dataset | Chen et al. [21] | GEO: GSE167339 |
| GSE156326 scar-formation dataset | Vorstandlechner et al. [22] | GEO: GSE156326 |
| GSE181316 keloid/scar fibroblast dataset | Direder et al. [23] | GEO: GSE181316 |
| FibroDynMix R package source | This paper | Local R package source; version 0.0.0.9000 |
| R statistical environment | R Core Team [30] | R 4.6.0; https://www.r-project.org/ |
| Single-cell data conventions | Amezquita et al. [28] | Bioconductor-compatible single-cell object conventions |
| Analysis outputs and figure source data | This paper | `analysis/`; `figures/` |

### Resource availability

#### Lead contact

Requests for resources and further information will be directed to the
corresponding author listed in the final submitted manuscript.

#### Materials availability

This study did not generate new biological materials.

#### Data and code availability

This paper analyzes existing, publicly available data. These accession numbers
and DOIs are listed in the key resources table. Public mouse breast fibroblast
count data are from DOI 10.6071/M3238R [19,20]. Public human analyses used
GSE246215, GSE167339, GSE156326, and GSE181316-derived project inputs as
documented in the analysis manifests [3,21-23].

All original code has been deposited as an R package and is publicly available as
of the date of publication. The package source, analysis scripts, generated
analysis summaries, figure source data, manifests, legends, and export QC files
are available at https://github.com/jiangcongxin/FibroDynMix and permanently
archived at https://doi.org/10.5281/zenodo.20787527. Generated analyses are
stored under `analysis/`, and figure source data, manifests, legends, and export
QC are stored under `figures/`. Large public count files are excluded from the
source package tarball through `.Rbuildignore` and are available from the public
repositories and accessions listed above.
Analyses followed R/Bioconductor-compatible single-cell data conventions and were
run in R 4.6.0 with project-local package and runtime records [28,30,38].

Any additional information required to reanalyze the data reported in this paper
is available from the lead contact upon request.

### Experimental model and study participant details

This study reanalyzed public, previously generated single-cell fibroblast count
matrices. No new human participants, animals, cell lines, or biological
specimens were enrolled, generated, or collected for this work. Public datasets
were used at the count-matrix or processed-count level documented in the
corresponding source records and local analysis manifests [3,19-23].
GSE167339 Human3 was retained only for donor robustness because the public files
do not map HumanHashTag IDs to treatment labels.

### Method details

#### FibroDynMix count model

FibroDynMix models gene-by-cell UMI counts directly rather than treating
normalized expression as the primary response [35]. The current implementation
uses a negative-binomial likelihood. For cell `i` and gene `g`, the log mean is
modeled as a library-size offset plus a gene baseline, latent state-program
contribution, and optional study and donor effects [14-16]:

```text
x_ig ~ NB(mu_ig, phi_g)
log(mu_ig) = log(l_i) + alpha_g + z_i beta_g + b_sg + u_jg
```

The state vector `z_i` lies on the simplex, so state weights are interpreted as
compositional quantities [31,32]. Marker sets orient the state programs through
weak penalties and alignment terms but do not directly define the final state
weights.

#### Model fitting

The fitting workflow initializes state weights and state-gene programs, then
alternates gene-level coefficient updates and cell-level state-logit updates
under the NB objective. The optimizer records objective values, executed
iterations, rollback status, stopping reason, and study or donor effect norms.
Optional study and donor effects are fitted with ridge penalties.

#### Simulation benchmark

The simulator generated four scenarios: continuous state mixtures, discrete
states, batch confounding, and rare transition-like cells. Benchmarks compared
marker scoring, NMF/topic modeling, the FibroDynMix initializer, FibroDynMix NB,
study-effect NB, and FibroDynMix VI [8-13]. The scVI latent baseline was
included in the local scVI run and projected to state-aligned simplex features
for comparison. Metrics included state-weight RMSE, dominant-state accuracy,
entropy preservation, marker-program recovery, NB optimizer
diagnostics, and donor-grouped downstream disease-label classification.

#### Validation-aware NB selection

Higher-outer NB runs evaluated `n_outer = 10` and `20`; the validation-aware
selector evaluated `n_outer = 2, 5, 10, 20`. Candidate fits were scored with
held-out NB objective, state-weight stability, marker-gradient preservation, and
downstream validation. Training objective was reported but not used alone for
final selection. Simulated truth metrics were retained as audit columns.

#### Posterior and uncertainty summaries

Bootstrap uncertainty was estimated by resampling cells and refitting the
selected model [33]. The VI posterior fixed the fitted NB program and used a
diagonal logistic-normal approximation over cell state logits [31,34].
Posterior draws were transformed through softmax to compute state-weight,
entropy, and sample-composition intervals. In simulations, interval coverage was
evaluated against known state weights, with post-hoc interval-scale calibration
reported as a diagnostic. Bootstrap outputs summarize state composition, cell
entropy, marker-program stability, and entropy-associated state-weight interval
width.

#### Public count-matrix analyses

Public mouse breast fibroblast and CAF count matrices from DOI 10.6071/M3238R
were used for raw-count execution, state composition, FPI, transition flow, and
bidirectional transfer analyses [19,20]. GSE246215 was split into GC, HCC, NSCLC,
and TNBC fibroblast subsets and evaluated with pooled fitting, leave-dataset-out
transfer, downsampling sensitivity, q95 library-size trimming, and downstream
cancer-type classification from patient- or sample-level features [3]. GSE167339,
GSE156326, and GSE181316 were used for bounded external fibroblast validation
and Figure 5 [21-23].

#### Transition flow and FPI

State-program transcriptional cost was computed from fitted state programs.
Condition-level state compositions were connected with entropy-regularized
optimal transport [24-26]. FPI combined normalized state-mixture entropy with
transition-potential contribution from the inferred flow.

### Quantification and statistical analysis

Simulation benchmarks used known latent state weights as truth for RMSE and
dominant-state accuracy. Downstream classification used grouped
cross-validation at donor, patient, or sample level depending on the analysis,
preserving sample-level grouping where available. GSE246215 biological-gradient
checks used Spearman correlations between marker gradients and expected
FibroDynMix state weights. Scar-module validation used small-sample pseudobulk
effect summaries, following sample-level single-cell differential-state
principles [36,37]. Figure 5 module contrasts were treated as directional
effect-size evidence because the adjusted q values did not pass FDR significance
[29].

### Quality gates

The project validation gate consists of:

```bash
Rscript scripts/run_project_maturity_audit.R
Rscript scripts/check_project_integrity.R
Rscript -e 'testthat::test_local()'
R CMD build .
R CMD check --no-manual --no-build-vignettes FibroDynMix_0.0.0.9000.tar.gz
```

The latest local validation reported complete figure source data, a passing
project integrity check, 359 passing unit tests, and `R CMD check` status OK.

## Acknowledgments

Funding, institutional support, and data-provider acknowledgments will be
added by the authors before submission.

## Author contributions

Author contribution statements will be added after the author list is fixed.

## Declaration of interests

The declaration of interests will be completed by the authors before
submission.

## References

1. Buechler, M.B., Pradhan, R.N., Krishnamurty, A.T., Cox, C.,
   Calviello, A.K., Wang, A.W., et al. (2021). Cross-tissue organization of the
   fibroblast lineage. Nature 593, 575-579.
   https://doi.org/10.1038/s41586-021-03549-5.
2. Lynch, M.D., and Watt, F.M. (2018). Fibroblast heterogeneity: implications
   for human disease. J. Clin. Invest. 128, 26-35.
   https://doi.org/10.1172/JCI93555.
3. Gao, Y., Li, J., Cheng, W., Diao, T., Liu, H., Bo, Y., et al. (2024).
   Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with
   distinct roles in immune modulation. Cancer Cell 42, 1764-1783.e10.
   https://doi.org/10.1016/j.ccell.2024.08.020.
4. Korsunsky, I., Wei, K., Pohin, M., Kim, E.Y., Barone, F., Major, T.,
   et al. (2022). Cross-tissue, single-cell stromal atlas identifies shared
   pathological fibroblast phenotypes in four chronic inflammatory diseases.
   Med 3, 481-518.e14. https://doi.org/10.1016/j.medj.2022.05.002.
5. Kuppe, C., Ibrahim, M.M., Kranz, J., Zhang, X., Ziegler, S.,
   Perales-Paton, J., et al. (2021). Decoding myofibroblast origins in human
   kidney fibrosis. Nature 589, 281-286.
   https://doi.org/10.1038/s41586-020-2941-1.
6. Guerrero-Juarez, C.F., Dedhia, P.H., Jin, S., Ruiz-Vega, R., Ma, D.,
   Liu, Y., et al. (2019). Single-cell analysis reveals fibroblast
   heterogeneity and myeloid-derived adipocyte progenitors in murine skin
   wounds. Nat. Commun. 10, 650.
   https://doi.org/10.1038/s41467-018-08247-x.
7. Valenzi, E., Bulik, M., Tabib, T., Morse, C., Sembrat, J.,
   Trejo Bittar, H., et al. (2019). Single-cell analysis reveals fibroblast
   heterogeneity and myofibroblasts in systemic sclerosis-associated
   interstitial lung disease. Ann. Rheum. Dis. 78, 1379-1387.
   https://doi.org/10.1136/annrheumdis-2018-214865.
8. Satija, R., Farrell, J.A., Gennert, D., Schier, A.F., and Regev, A.
   (2015). Spatial reconstruction of single-cell gene expression data.
   Nat. Biotechnol. 33, 495-502. https://doi.org/10.1038/nbt.3192.
9. Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E.,
   Mauck, W.M., et al. (2019). Comprehensive integration of single-cell data.
   Cell 177, 1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031.
10. Lee, D.D., and Seung, H.S. (1999). Learning the parts of objects by
    non-negative matrix factorization. Nature 401, 788-791.
    https://doi.org/10.1038/44565.
11. Gaujoux, R., and Seoighe, C. (2010). A flexible R package for nonnegative
    matrix factorization. BMC Bioinformatics 11, 367.
    https://doi.org/10.1186/1471-2105-11-367.
12. Kotliar, D., Veres, A., Nagy, M.A., Tabrizi, S., Hodis, E.,
    Melton, D.A., and Sabeti, P.C. (2019). Identifying gene expression programs
    of cell-type identity and cellular activity with single-cell RNA-Seq.
    eLife 8, e43803. https://doi.org/10.7554/eLife.43803.
13. Lopez, R., Regier, J., Cole, M.B., Jordan, M.I., and Yosef, N. (2018).
    Deep generative modeling for single-cell transcriptomics. Nat. Methods 15,
    1053-1058. https://doi.org/10.1038/s41592-018-0229-2.
14. Hafemeister, C., and Satija, R. (2019). Normalization and variance
    stabilization of single-cell RNA-seq data using regularized negative
    binomial regression. Genome Biol. 20, 296.
    https://doi.org/10.1186/s13059-019-1874-1.
15. Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold
    change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550.
    https://doi.org/10.1186/s13059-014-0550-8.
16. Ahlmann-Eltze, C., and Huber, W. (2021). glmGamPoi: fitting Gamma-Poisson
    generalized linear models on single cell count data. Bioinformatics 36,
    5701-5702. https://doi.org/10.1093/bioinformatics/btaa1009.
17. Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F.,
    Wei, K., et al. (2019). Fast, sensitive and accurate integration of
    single-cell data with Harmony. Nat. Methods 16, 1289-1296.
    https://doi.org/10.1038/s41592-019-0619-0.
18. Haghverdi, L., Lun, A.T.L., Morgan, M.D., and Marioni, J.C. (2018).
    Batch effects in single-cell RNA-sequencing data are corrected by matching
    mutual nearest neighbors. Nat. Biotechnol. 36, 421-427.
    https://doi.org/10.1038/nbt.4091.
19. Sebastian, A., Hum, N.R., Martin, K.A., Gilmore, S.F., Peran, I.,
    Byers, S.W., et al. (2020). Single-cell transcriptomic analysis of
    tumor-derived fibroblasts and normal tissue-resident fibroblasts reveals
    fibroblast heterogeneity in breast cancer. Cancers 12, 1307.
    https://doi.org/10.3390/cancers12051307.
20. Hum, N., Sebastian, A., Martin, K., Gilmore, S., Byers, S., Wheeler, E.,
    et al. (2020). Data from: Single-cell transcriptomic analysis of
    tumor-derived fibroblasts and normal tissue-resident fibroblasts reveals
    fibroblast heterogeneity in breast cancer. Dryad Dataset.
    https://doi.org/10.6071/M3238R.
21. Chen, K., Kwon, S.H., Henn, D., Kuehlmann, B.A., Tevlin, R.,
    Bonham, C.A., et al. (2021). Disrupting biological sensors of force promotes
    tissue regeneration in large organisms. Nat. Commun. 12, 5256.
    https://doi.org/10.1038/s41467-021-25410-z.
22. Vorstandlechner, V., Laggner, M., Copic, D., Klas, K., Direder, M.,
    Chen, Y., et al. (2021). The serine proteases dipeptidyl-peptidase 4 and
    urokinase are key molecules in human and mouse scar formation. Nat. Commun.
    12, 6242. https://doi.org/10.1038/s41467-021-26495-2.
23. Direder, M., Weiss, T., Copic, D., Vorstandlechner, V., Laggner, M.,
    Pfisterer, K., et al. (2022). Schwann cells contribute to keloid formation.
    Matrix Biol. 108, 55-76. https://doi.org/10.1016/j.matbio.2022.03.001.
24. Schiebinger, G., Shu, J., Tabaka, M., Cleary, B., Subramanian, V.,
    Solomon, A., et al. (2019). Optimal-transport analysis of single-cell gene
    expression identifies developmental trajectories in reprogramming. Cell
    176, 928-943.e22. https://doi.org/10.1016/j.cell.2019.01.006.
25. Cuturi, M. (2013). Sinkhorn distances: Lightspeed computation of optimal
    transport. Adv. Neural Inf. Process. Syst. 26, 2292-2300.
26. Peyre, G., and Cuturi, M. (2019). Computational optimal transport with
    applications to data sciences. Found. Trends Mach. Learn. 11, 355-607.
    https://doi.org/10.1561/2200000073.
27. Barallobre-Barreiro, J., Woods, E., Bell, R.E., Easton, J.A., Hobbs, C.,
    Eager, M., et al. (2019). Cartilage-like composition of keloid scar
    extracellular matrix suggests fibroblast mis-differentiation in disease.
    Matrix Biol. Plus 4, 100016.
    https://doi.org/10.1016/j.mbplus.2019.100016.
28. Amezquita, R.A., Lun, A.T.L., Becht, E., Carey, V.J., Carpp, L.N.,
    Geistlinger, L., et al. (2020). Orchestrating single-cell analysis with
    Bioconductor. Nat. Methods 17, 137-145.
    https://doi.org/10.1038/s41592-019-0654-x.
29. Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
    rate: a practical and powerful approach to multiple testing. J. R. Stat.
    Soc. Series B Stat. Methodol. 57, 289-300.
    https://doi.org/10.1111/j.2517-6161.1995.tb02031.x.
30. R Core Team (2026). R: A Language and Environment for Statistical
    Computing. R Foundation for Statistical Computing, Vienna, Austria.
31. Aitchison, J., and Shen, S.M. (1980). Logistic-normal distributions: Some
    properties and uses. Biometrika 67, 261-272.
    https://doi.org/10.1093/biomet/67.2.261.
32. Aitchison, J. (1982). The statistical analysis of compositional data. J. R.
    Stat. Soc. Series B Stat. Methodol. 44, 139-160.
    https://doi.org/10.1111/j.2517-6161.1982.tb01195.x.
33. Efron, B. (1979). Bootstrap methods: another look at the jackknife. Ann.
    Stat. 7, 1-26. https://doi.org/10.1214/aos/1176344552.
34. Blei, D.M., Kucukelbir, A., and McAuliffe, J.D. (2017). Variational
    inference: a review for statisticians. J. Am. Stat. Assoc. 112, 859-877.
    https://doi.org/10.1080/01621459.2017.1285773.
35. Townes, F.W., Hicks, S.C., Aryee, M.J., and Irizarry, R.A. (2019). Feature
    selection and dimension reduction for single-cell RNA-Seq based on a
    multinomial model. Genome Biol. 20, 295.
    https://doi.org/10.1186/s13059-019-1861-6.
36. Soneson, C., and Robinson, M.D. (2018). Bias, robustness and scalability in
    single-cell differential expression analysis. Nat. Methods 15, 255-261.
    https://doi.org/10.1038/nmeth.4612.
37. Crowell, H.L., Soneson, C., Germain, P.-L., Calini, D., Collin, L.,
    Raposo, C., et al. (2020). muscat detects subpopulation-specific state
    transitions from multi-sample multi-condition single-cell transcriptomics
    data. Nat. Commun. 11, 6077.
    https://doi.org/10.1038/s41467-020-19894-4.
38. Luecken, M.D., and Theis, F.J. (2019). Current best practices in
    single-cell RNA-seq analysis: a tutorial. Mol. Syst. Biol. 15, e8746.
    https://doi.org/10.15252/msb.20188746.
