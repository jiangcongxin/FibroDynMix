# Cover Letter Draft

Dear iScience Editors,

Please consider our manuscript, "FibroDynMix models fibroblast state plasticity
as raw-count latent mixtures," for publication in iScience. The manuscript
presents FibroDynMix, an R package for estimating fibroblast state mixtures
directly from raw single-cell UMI counts using a negative-binomial model with
weak marker orientation, simplex state weights, study/donor terms, uncertainty
summaries, transfer, and cross-sectional transition-flow analysis.

FibroDynMix addresses a common limitation in fibroblast state analysis. Marker
scores are useful for annotation, but they do not estimate mixed state weights
under a raw count-generating model or separate state structure from study- and
donor-level expression shifts. In simulations with known latent truth,
FibroDynMix-family features reduced mean state-weight RMSE versus marker/NMF
baselines in the extended benchmark and improved downstream balanced accuracy.
Converged negative-binomial reruns showed that lower training objectives can
worsen state recovery, motivating validation-aware model selection.

The manuscript also tests the package on public fibroblast count matrices.
Public mouse fibroblast matrices support raw-count execution and bidirectional
transfer. Public human datasets support transfer, expected fibroblast
marker-gradient recovery, and directional scar-module shifts, while the
manuscript explicitly preserves the boundary that these scar-module contrasts
are not FDR-significant and do not establish diagnostic, causal, therapeutic, or
lineage-tracing claims. Each main figure is linked to source data, figure
manifests, and claim boundaries.

We believe this method-first, source-backed study is a strong fit for iScience
because it combines computational model development, single-cell fibroblast
biology, public data reuse, and reproducible research packaging. FibroDynMix
will be of interest to readers who use single-cell data to quantify cell-state
plasticity across disease contexts and need auditable methods that separate
model fit, state recovery, and biological interpretation.

This manuscript is original, is not under consideration elsewhere, and all
authors will approve the final submitted version. The study reanalyzes public,
previously generated datasets and did not generate new human participants,
animals, cell lines, or biological specimens. A final declaration of interests,
funding statement, author contributions, and data/code availability links will
be supplied after author confirmation.

Sincerely,

[Corresponding author name]  
[Institution]  
[Email]

