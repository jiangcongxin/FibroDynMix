# Reference Article Structure Matrix for FibroDynMix

## Purpose

This matrix compares target-journal article structures before rewriting the
FibroDynMix manuscript. The goal is structural alignment, not textual imitation.
The selected papers emphasize computational methods, single-cell analysis,
Bayesian modeling, benchmarking, or resource-style presentation in iScience and
Communications Biology.

## Journal-Level Structure Notes

| Journal | Observed structure | Implication for FibroDynMix |
|---|---|---|
| iScience | Summary, graphical abstract, highlights, introduction, result-led sections, discussion, limitations of the study, STAR Methods, resource availability, data/code availability. | Use Summary and Highlights, then make each Results heading a direct figure-level claim. Put implementation and reproducibility details in STAR Methods rather than Results. |
| Communications Biology | Abstract, Introduction, Results, Discussion, Methods, Reporting Summary, Data availability, Code availability. | Useful backup style for Methods and availability statements, but the main draft should follow iScience because the current evidence is method-first with bounded biological validation. |

## Reference Matrix

| Paper | Journal | Article type for our purposes | Results-section pattern | Discussion/limitation pattern | Structural lesson for FibroDynMix |
|---|---|---|---|---|---|
| SCING: Inference of robust, interpretable gene regulatory networks from single cell and spatial transcriptomics | iScience, 2023, DOI: 10.1016/j.isci.2023.107124 | Computational method with benchmarking and applications | Starts with method overview, then speed/capacity, perturb-seq validation, robustness, network properties, disease subnetworks, and three application cases. | Discussion follows method performance and applications; limitations are explicit and separated. | Start FibroDynMix Results with the model schematic, then move through benchmarks, robustness, public-data applications, and bounded biological use. |
| Single-cell Bayesian deconvolution | iScience, 2023, DOI: 10.1016/j.isci.2023.107941 | Bayesian method with synthetic and experimental validation | Defines the theoretical problem first, then tests synthetic data, external-noise data, and multidimensional distributions. | Limitations are separated from Discussion. | Use a clean sequence: model problem, simulation truth, public count matrices, then external validation. Keep Bayesian boundaries in Limitations, not repeatedly in each paragraph. |
| SC-MO-GRN-DB: A comprehensive repository for single-cell multiomic gene regulatory networks | iScience, 2026, DOI: 10.1016/j.isci.2026.115323 | Resource/database article with method curation | Uses compact Results headings: reference networks and single-cell datasets. | Discussion connects resource scope to current gaps; limitations are explicit. | Keep package/resource maturity concise. Do not over-explain every audit file in Results; move package reproducibility to STAR Methods and Resource availability. |
| Single-cell aging trajectories reveal a dynamic coupling between nuclear size and proteasome concentration | iScience, 2026, DOI: 10.1016/j.isci.2026.114736 | Result-led single-cell analysis article | Each Results heading names the experiment and the finding: tracking dynamics, identifying predictors, modeling transport, and testing perturbation. | Discussion starts from findings rather than a broad review. | Write FibroDynMix Results headings as claims, not generic labels. Use dataset or analysis as the subject. |
| scPML: pathway-based multi-view learning for cell type annotation from single-cell RNA-seq data | Communications Biology, 2023, DOI: 10.1038/s42003-023-05634-z | Single-cell machine-learning method | Overview first, then cross-platform, cross-species, multi-view benefit, multiple training data, and unknown cell-type detection. | Methods are detailed and separated; statistics and reproducibility are explicit. | Keep comparator logic clear: marker scoring, NMF, scVI, initializer, NB, study-effect NB, and VI. Separate benchmark settings from claims. |
| VBASS enables integration of single cell gene expression data in Bayesian association analysis of rare variants | Communications Biology, 2023, DOI: 10.1038/s42003-023-05155-9 | Bayesian model with simulations and real disease examples | Model overview, simulation with bulk expression, simulation with single-cell expression, then disease applications. | Discussion explains why the model gains power and where data assumptions remain. | For FibroDynMix, pair each model claim with simulation evidence before public-data examples. |
| Advancing spatial cellular communication inference with ligand diffusion and transport model | Communications Biology, 2026, DOI: 10.1038/s42003-025-09413-w | Spatial single-cell computational method | Overview, benchmark superiority, global optimization interpretation, marker identification, subdomain analysis, and external validation. | Methods contain mathematical formulation and evaluation metrics. | FibroDynMix can use the same order: overview, benchmark, optimization caveat, marker-gradient sanity check, transfer, external validation. |
| SpaDC enables sequence-based integrative analysis and regulatory inference of spatial chromatin accessibility data | Communications Biology, 2026, DOI: 10.1038/s42003-026-10462-y | Recent computational model article | Abstract pattern: problem, omitted data modality, model, benchmark tasks, biological application. | The available early-access page emphasizes abstract-level claim density. | FibroDynMix Summary should be compressed to problem, model, simulations, public data, external validation, and reproducibility. |

## Common Structural Moves to Adopt

1. Start with a concrete unmet analytical problem.
2. State the method in one sentence before naming every component.
3. Put method overview before benchmark results.
4. Use simulation truth as the primary validation for latent variables.
5. Report where baselines remain competitive.
6. Use public datasets as application and stress-test evidence, not as overbroad biological proof.
7. Keep limitations in a dedicated section and make them direct.
8. Move implementation details, commands, and reproducibility checks into STAR Methods.

## Target Outline for the Rewritten FibroDynMix Manuscript

1. Title
2. Summary
3. Highlights
4. Graphical abstract placeholder
5. Introduction
6. Results
   - FibroDynMix estimates simplex fibroblast states from raw UMI counts.
   - Public raw-count matrices test package execution and bidirectional transfer.
   - Simulation truth separates state recovery from marker and topic baselines.
   - Validation-aware selection separates NB fit from state recovery.
   - Study-effect terms reduce simulated cohort-shift sensitivity.
   - Human fibroblast datasets support transfer and bounded biological utility.
   - Cross-sectional transition flow links state mixtures to FPI.
7. Discussion
8. Limitations of the study
9. STAR Methods
10. Resource availability
11. Data and code availability

## Sources

- iScience author information: https://www.cell.com/iscience/information-for-authors
- iScience article types: https://www.cell.com/iscience/information-for-authors/article-types
- Communications Biology submission guidelines: https://www.nature.com/commsbio/submit/submission-guidelines
- SCING: https://pmc.ncbi.nlm.nih.gov/articles/PMC10331489/
- Single-cell Bayesian deconvolution: https://pmc.ncbi.nlm.nih.gov/articles/PMC10579429/
- SC-MO-GRN-DB: https://pmc.ncbi.nlm.nih.gov/articles/PMC13018865/
- Single-cell aging trajectories: https://pmc.ncbi.nlm.nih.gov/articles/PMC12955095/
- scPML: https://pmc.ncbi.nlm.nih.gov/articles/PMC10721875/
- VBASS: https://pmc.ncbi.nlm.nih.gov/articles/PMC10368729/
- SCILD: https://pmc.ncbi.nlm.nih.gov/articles/PMC12855878/
- SpaDC: https://www.nature.com/articles/s42003-026-10462-y
