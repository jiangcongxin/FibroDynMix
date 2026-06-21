# Bioinformatics Validation

This analysis adds lightweight biological validation around the FibroDynMix
real-data results. It is designed to answer three reviewer-facing questions:

- Are the public inputs fibroblast-enriched rather than obvious contaminating
  cell types?
- Do inferred state programs map to interpretable fibroblast biology?
- Are state abundance and plasticity summaries reproducible at donor/sample
  level?

## Run

```bash
Rscript scripts/run_bioinformatics_validation.R
```

Outputs:

- `analysis/bioinformatics_validation/bioinformatics_validation_manifest.tsv`
- `analysis/bioinformatics_validation/fibroblast_purity_qc.tsv`
- `analysis/bioinformatics_validation/state_associated_top_genes.tsv`
- `analysis/bioinformatics_validation/state_program_pathway_enrichment.tsv`
- `analysis/bioinformatics_validation/donor_state_abundance_fpi.tsv`
- `analysis/bioinformatics_validation/donor_aware_state_effects.tsv`

## Current Results

Current run:

- validation datasets: 2 (`GSE246215`, `GSE167339`)
- QC rows: 15 sample/dataset rows
- minimum mean fibroblast purity margin: 0.823
- maximum low-purity fraction: 0.175
- state-associated top-gene rows: 360
- pathway-enrichment rows: 72
- enriched state/pathway pairs at q <= 0.05: 4
- donor-aware state-effect rows: 12

The strongest curated pathway checks support:

- `GSE246215` ECM-remodeling: ECM organization genes
  (`COL1A1`, `COL5A1`, `POSTN`, `COL3A1`, `MMP14`)
- `GSE246215` resident: resident matrix-homeostasis genes
  (`DCN`, `LUM`, `FBLN1`, `PDGFRA`)
- `GSE167339` ECM-remodeling: ECM organization genes
  (`COL1A1`, `COL3A1`, `COL1A2`)

For GSE167339, donor-aware normal-to-perturbed effects are directionally
consistent across Human1 and Human2:

- resident composition increases;
- inflammatory composition increases;
- antigen-presenting composition increases;
- ECM-remodeling composition decreases;
- entropy/FPI increases.

Human3 remains `hash_unknown` and is used for donor-level transfer robustness,
not treatment-specific biological interpretation.

## Claim Boundary

This is a curated marker/pathway validation layer, not exhaustive cell-type
annotation, full GO/Reactome enrichment, spatial validation, or causal disease
biology. The analysis supports biological plausibility and input quality for
the current public-data validation layer.
