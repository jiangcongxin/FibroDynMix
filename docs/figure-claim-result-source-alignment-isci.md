# Figure Claim-Result-Source Alignment Check for FibroDynMix

## Basic Information

- Paper title: FibroDynMix models fibroblast state plasticity as raw-count latent mixtures
- Scope of check: Main Figures 1-6 plus manuscript-facing supplementary figure packages
- Date of check: 2026-06-15

## Consistency Issues List

- Key points summary: None found.
- Full issue list: `docs/figure-claim-result-source-alignment-issues.csv`
- The manuscript draft has no PDF page markers. The current CSV contains only
  the header because no open figure-text-source issues remain after assigning
  Supplementary Figures S1-S4.

## Summary

- Presence of serious inconsistencies: No
- Overview of main issues: None found. The main and supplementary figure spine
  is source-backed and cited from the manuscript draft.
- Note: This report compares manuscript Results text, figure manifests, figure
  legends where available, and panel source-data manifests. It does not inspect
  the rendered image content.

## Main Figure Alignment Matrix

| Figure | Manifest claim | Results paragraph alignment | Source-data manifest | Status | Notes |
|---|---|---|---|---|---|
| Figure 1 | Raw-count generative mixture architecture with latent state simplex, study/donor effects, uncertainty, and transition flow. | Results section "FibroDynMix estimates simplex fibroblast states from raw UMI counts" cites Figure 1 and describes the count matrix, simplex state weights, NB likelihood, weak marker orientation, effects, uncertainty, transfer, and transition outputs. | `figures/figure1/panel_source_data_manifest.tsv` maps panels A-D to model nodes/edges, likelihood terms, evidence map, and figure inventory. | Pass | Architecture boundary is preserved; Figure 1 is not used as standalone biological validation. |
| Figure 2 | Public mouse raw-count fibroblast subset supports real-data execution, state composition, FPI, and transition-flow readouts. | Results section "Public raw-count matrices test package execution and bidirectional transfer" reports 600 cells, 1,217 genes, NB objective improvement 0.0900, and transition convergence. | `figures/figure2/panel_source_data_manifest.tsv` maps panels A-G to condition counts, NB objective, state composition, FPI, transition flow, marker coverage, and cell weights. | Pass | Numeric values align with `figures/figure2/figure_manifest.tsv` and source data. |
| Figure 3 | Simulation benchmark quantifies state-mixture recovery, VI behavior, confounding robustness, rare-transition detection, and marker-program recovery. | Results now states that Figure 3 visualizes core simulation readouts and separately describes the extended nonvisual NMF/topic benchmark. | `figures/figure3/panel_source_data_manifest.tsv` maps panels A-F to true composition, benchmark metrics, entropy calibration, rare-transition scores, and marker recovery. | Pass | This avoids incorrectly assigning extended benchmark numbers to Figure 3 panels. |
| Figure 4 | Ridge-penalized study effect can absorb cohort-level expression shifts while retaining state-mixture recovery within the selected penalty range. | Results section "Study-effect terms reduce simulated cohort-shift sensitivity" reports study L2 = 5 and marker L2 = 0.05. | `figures/figure4/panel_source_data_manifest.tsv` maps panels A-D to RMSE/objective tradeoff, study-effect magnitude, and selected penalty. | Pass | Optional final polish can add the selection-score sentence from the legend if space allows. |
| Figure 5 | External public datasets support bounded reproducibility and biological utility across GSE167339, GSE156326, and GSE181316. | Results section "Human fibroblast datasets support transfer and bounded biological utility" reports GSE167339 leave-donor convergence, GSE156326/GSE181316 transfer convergence, scar-module deltas, AUCs, and q values. | `figures/figure5/panel_source_data_manifest.tsv` maps panels A-D to donor composition, leave-donor transfer, external state composition, and scar module validation. | Pass | Boundary is consistent: Human3 is hash-unknown robustness evidence; scar module q values are not FDR-significant. |
| Figure 6 | Cross-sectional state flow and FPI are inferred from simulated state mixtures, not lineage tracing. | Results section "Cross-sectional transition flow links state mixtures to FPI" reports entropy-regularized optimal transport, flow entropy 2.1602, expected cost 15.3812, convergence, and FPI interpretation. | `figures/figure6/panel_source_data_manifest.tsv` maps panels A-D to condition state composition, transition flow, FPI, and rare-transition labels. | Pass | Dynamic/lineage boundary is concentrated in Limitations and figure manifest. |

## Supplementary Figure Package Alignment

| Figure package | Current manuscript role | Source-data manifest | Status | Action before submission |
|---|---|---|---|---|
| `figures/supplementary_uncertainty` | Supplementary Figure S1 supports bootstrap uncertainty statements in Results and STAR Methods. | Present. | Pass | None. |
| `figures/public_realdata_smoke` | Supplementary Figure S2 supports the smaller public real-data smoke workflow. | Present. | Pass | None. |
| `figures/core_converged_tradeoff` | Supplementary Figure S3 supports the 16/16 objective-RMSE tradeoff claim. | Present. | Pass | None. |
| `figures/core_converged_tradeoff_scvi_r5` | Supplementary Figure S4 supports the 40/40 scVI-inclusive objective-RMSE tradeoff claim. | Present. | Pass | None. |

## Remarks

The manuscript is ready for language polishing. Main Figures 1-6 and
Supplementary Figures S1-S4 each have a figure claim boundary and panel-level
source-data manifest.
