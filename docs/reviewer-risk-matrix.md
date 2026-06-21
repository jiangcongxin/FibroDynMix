# Reviewer Risk Matrix

## Purpose

This matrix lists likely reviewer concerns for an iScience / Communications
Biology submission and points each concern to current evidence or remaining work.

| Reviewer concern | Current evidence | Residual risk | Response strategy |
|---|---|---|---|
| The method is only gene-set scoring. | Raw-count NB likelihood, latent simplex optimization, weak marker priors, Figure 1, Figure 3. | Marker priors still orient states. | Emphasize priors guide identifiability; z is inferred under count likelihood. |
| Marker scoring performs competitively in state-weight RMSE. | Extended method benchmark, marker-stress benchmark, and core converged benchmark report marker scoring, NMF/topic baselines, and FibroDynMix-family metrics with source tables. | Marker-aligned simulations naturally favor marker scoring for some RMSE summaries, and high-outer NB optimization does not currently rescue RMSE. | Do not claim uniform superiority over scoring; claim broader model outputs: count likelihood, hierarchy, uncertainty, transfer, transition diagnostics, and bounded downstream utility. |
| Marker scoring is enough; FibroDynMix z has no downstream value. | Simulation downstream benchmark shows higher mean FibroDynMix-family balanced accuracy and macro-F1 than the marker/NMF baseline average. GSE246215 downstream benchmark shows FibroDynMix z improves over raw marker-score features on macro-F1/AUROC and recovers expected myofibroblast, ECM-remodeling, and inflammatory gradients. | GSE246215 cancer-type classification is not a universal FibroDynMix win: NMF topics are strongest, and marker-score composition has higher balanced accuracy than FibroDynMix z in that task. | Use downstream evidence as a bounded rebuttal: z is useful and biologically interpretable, but not uniformly superior for every classifier or endpoint. |
| Ordinary topic models explain the results. | `fit_topic_nmf_baseline()` uses the mature `NMF` package when available and shows topic/NMF baseline metrics in `analysis/extended_method_benchmark`. | NMF can be sensitive to rank, initialization, and topic-state alignment. | Present NMF as a comparator and keep the claim focused on generative hierarchy and uncertainty. |
| The method ignores batch or cohort effects. | Study/donor effects in `fit_fibrodynmix_nb()`, Figure 4, study-effect sensitivity outputs. | Ridge effects are not full random-effect posterior. | State current hierarchy layer and future full posterior extension. |
| The uncertainty is not Bayesian enough. | Bootstrap uncertainty, `fit_fibrodynmix_vi()`, VI benchmark, calibrated coverage. | VI is over z logits around NB mode, not all parameters. | Present as lightweight posterior skeleton and explicitly bound claims. |
| Transition flow implies lineage tracing. | Figure 6 and docs state cross-sectional optimal transport boundary. | Readers may overinterpret arrows. | Use language: state-flow summary, not observed transition. |
| Public data results are underpowered or pooled. | Figure 2, public smoke figure, Dryad multi-public validation, and independent GSE246215 human validation disclose DOI/accession, cell/gene counts, state composition, and claim boundary. | GSE246215 validation uses processed public count matrices and sampled cancer-type subsets from one GEO study. | Frame as public count-matrix validation and pressure testing, not broad disease discovery. |
| HCC state composition may be a library-size artifact. | `analysis/gse246215_sensitivity` tests random downsampling and q95 library-size trimming; HCC remains a library-size outlier. | Simple q95 trimming does not fully normalize the HCC library-size distribution. | Treat HCC-specific composition as QC-sensitive and hypothesis-generating; avoid cancer-type biology claims. |
| Cross-cohort transfer is only simulated. | Simulation transfer, bidirectional public transfer smoke, Dryad leave-dataset-out transfer, GSE246215 leave-cancer-type-out transfer diagnostics, and GSE167339 donor transfer evidence. | Public transfer is not completed human cross-atlas generalization. | Keep transfer claim as mechanics and diagnostics; require harmonized independent cohorts before stronger atlas-level claims. |
| The R package is immature. | NEWS, CITATION, cran-comments, pkgdown config, function index, testthat suite passed, and R CMD check OK. | Development version 0.0.0.9000. | Present as reproducible research package; release version before submission. |
| Figures are not reproducible. | Figure source data, panel manifests, legends, export QC, project integrity gate. | Full figure regeneration depends on local analysis files. | Use `scripts/check_project_integrity.R` and audit outputs as reproducibility evidence. |
| Calibration is post-hoc. | VI benchmark reports raw and simulation-calibrated coverage/width. | Calibration uses known truth and cannot be applied directly to real data. | Present calibration only as simulation diagnostic. |
| Biological novelty is limited. | Fibroblast plasticity framing, FPI, transition flow, public smoke, and bounded external validation in GSE167339/GSE156326/GSE181316. | Method-first evidence remains stronger than disease-atlas evidence. | Target method-first framing and present external datasets as bounded biological-utility evidence. |

## Highest-Priority Remaining Risks

1. Validation-aware NB stopping/model selection is not yet finalized.
2. Atlas-level human multi-cohort generalization is not yet established.
3. Full amortized or all-parameter hierarchical VI is not yet implemented.
4. Public real-data analyses are smoke tests and case studies, not disease mechanism results.

These are not blockers for a method-first manuscript, but they constrain title,
abstract, and discussion language.
