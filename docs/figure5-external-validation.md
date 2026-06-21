# Figure 5 External Fibroblast Validation

This is the main Figure 5 package. It promotes the real-data layer from a
method demonstration to an external biological validation panel.

## Command

```bash
Rscript scripts/make_figure5.R
```

## Outputs

- `figures/figure5/exports/figure5.png`
- `figures/figure5/exports/figure5.pdf`
- `figures/figure5/exports/figure5.svg`
- `figures/figure5/exports/figure5.tiff`
- `figures/figure5/source_data/fig5a_gse167339_donor_state_composition.tsv`
- `figures/figure5/source_data/fig5b_gse167339_leave_donor_transfer.tsv`
- `figures/figure5/source_data/fig5c_scar_external_state_composition.tsv`
- `figures/figure5/source_data/fig5d_scar_module_validation.tsv`
- `figures/figure5/source_data/fig5_scar_transfer_diagnostics.tsv`

## Panel Logic

- Panel A: GSE167339 donor-level fibroblast state composition in Human1/Human2
  control versus perturbed collagen scar-system conditions.
- Panel B: GSE167339 leave-donor-out transfer convergence. Human3 is treated as
  hash-unknown donor robustness evidence.
- Panel C: GSE156326 and GSE181316 external scar cohort dominant-state
  composition after final-branch transfer.
- Panel D: GSE156326/GSE181316 disease-control scar-module score deltas from
  scRNA-seq pseudobulk validation, shown as descriptive effect-size evidence
  rather than FDR-significant testing.

## Current Readout

- GSE167339 leave-donor-out transfer convergence remains high:
  Human1 = 1.000, Human2 = 1.000, Human3 = 0.978.
- External scar transfer convergence is high for the two promoted cohorts:
  GSE156326 = 0.954 and GSE181316 = 0.960.
- Disease-control module deltas are positive in both scar cohorts but not
  FDR-significant after correction:
  GSE156326 delta = 0.028, AUC = 0.78, q = 0.38;
  GSE181316 delta = 0.085, AUC = 1.00, q = 0.061.

## Claim Boundary

This figure supports external reproducibility and biological utility across
independent public fibroblast datasets. It should not be presented as diagnostic
performance, causal mechanism, lineage tracing, or therapeutic efficacy.
GSE167339 Human3 hash groups support donor robustness only because public files
do not map hash IDs to treatment labels. GSE156326/GSE181316 module statistics
are small-sample pseudobulk validation and should be read as directional
effect-size support, not significance evidence.
