# Figure 6 Transition Flow and FPI

## Figure-Level Claim

Figure 6 adds the dynamic/plasticity layer to FibroDynMix:

```text
Latent fibroblast state compositions can be coupled with state-program
transcriptional costs to infer cross-sectional condition-associated state flow,
and cell-level mixed-state uncertainty can be summarized as a Fibroblast
Plasticity Index.
```

The claim boundary is explicit:

```text
The flow is inferred from cross-sectional simulated compositions and state
program costs. It is not lineage tracing, RNA velocity, or directly observed
temporal transition.
```

## Panels

- **A. Condition-level state composition**
  - Mean latent state fractions for normal and disease simulated conditions.

- **B. Inferred normal-to-disease state flow**
  - Entropy-regularized optimal transport from normal composition to disease
    composition using transcriptional state-program cost.

- **C. Fibroblast Plasticity Index by condition**
  - Cell-level FPI summarizes latent-state entropy plus transition potential.

- **D. FPI in simulator-labeled transition cells**
  - FPI distribution is evaluated against simulator-provided rare-transition
    labels.

## Source Data

- `figures/figure6/source_data/fig6_condition_state_composition.tsv`
- `figures/figure6/source_data/fig6_transition_flow_long.tsv`
- `figures/figure6/source_data/fig6_transition_cost_long.tsv`
- `figures/figure6/source_data/fig6_transition_flow_summary.tsv`
- `figures/figure6/source_data/fig6_cell_fpi.tsv`

## Role in the Manuscript

This figure answers the reviewer question:

```text
Does FibroDynMix model state plasticity rather than only assign static states?
```

For the current manuscript target, this figure should be presented as a
simulation and methods-validation figure. A real-data version can later reuse
the same functions after a real fibroblast dataset interface is added.
