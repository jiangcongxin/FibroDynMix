# FibroDynMix Function Index

## Simulation and Benchmarks

- `simulate_fibrodynmix()` creates raw-count generative simulations with known
  latent state truth.
- `run_simulation_benchmark()` compares marker scoring, initialization, NB
  optimization, topic/NMF baselines, and VI posterior means across benchmark
  scenarios.
- `summarize_benchmark_results()` aggregates benchmark metrics by scenario and
  method.
- `run_marker_scoring_benchmark()` evaluates direct marker scoring as a baseline.
- `fit_topic_nmf_baseline()` fits a mature-package NMF/topic comparator with
  marker-based topic-to-state alignment.

## Model Fitting

- `fit_fibrodynmix_initializer()` initializes latent state weights and
  state-gene programs from raw counts.
- `fit_fibrodynmix_nb()` fits the raw-count NB optimizer with optional study and
  donor effects.
- `fit_fibrodynmix_prepared()` fits from a validated real-data object.

## Likelihood and Objectives

- `fibrodynmix_nb_loglik()` evaluates the raw-count NB likelihood.
- `fibrodynmix_nb_deviance()` computes NB deviance.
- `fibrodynmix_nb_objective()` evaluates the penalized NB optimization objective.

## Uncertainty and Posterior Inference

- `bootstrap_fibrodynmix()` estimates cell, sample, and marker-program
  uncertainty by bootstrap.
- `fit_fibrodynmix_vi()` fits a lightweight logistic-normal posterior over cell
  state logits around the NB mode.
- `evaluate_posterior_intervals()` evaluates posterior state-interval coverage
  in simulation.
- `calibrate_posterior_interval_scale()` performs simulation-only interval
  scale calibration.

## Cross-Cohort Transfer

- `fit_fibrodynmix_transfer()` freezes a fitted NB program and optimizes held-out
  cell state weights.
- `run_cross_cohort_transfer_benchmark()` runs leave-study-out simulation
  transfer benchmarks.

## Transition and Plasticity

- `compute_state_cost()` computes transcriptional state-to-state cost.
- `estimate_transition_flow()` estimates entropy-regularized state flow.
- `compute_fpi()` computes cell-level Fibroblast Plasticity Index.

## Visualization

- `plot_state_composition()` visualizes sample, donor, condition, or dataset
  state fractions.
- `plot_cell_state_heatmap()` visualizes cell-level state mixture weights.
- `plot_transfer_diagnostics()` visualizes transfer convergence diagnostics.
- `plot_transition_flow()` visualizes source-to-target state flow.
- `plot_fpi_distribution()` visualizes FPI or entropy distributions.
- `plot_benchmark_rankings()` visualizes method performance across stress
  scenarios.
- `plot_purity_qc()` visualizes fibroblast purity QC margins.
- `plot_pathway_enrichment()` visualizes state-program pathway enrichment.
- `plot_fibroblast_annotation()` visualizes FibroDynMix or user-supplied
  fibroblast annotations on Seurat embeddings.
- `plot_fibroblast_marker_dot()` visualizes fibroblast marker expression by
  annotation group.

## Real Data

- `prepare_fibrodynmix_data()` validates raw counts, metadata, marker priors, and
  filtering decisions before model fitting.

## Seurat Integration

- `prepare_fibrodynmix_seurat()` extracts Seurat assay counts and metadata into
  a `FibroDynMixData` object.
- `fit_fibrodynmix_seurat()` runs the FibroDynMix fitting workflow from a Seurat
  object and optionally returns an annotated Seurat object.
- `add_fibrodynmix_to_seurat()` attaches state weights, dominant state,
  entropy, FPI, and an optional state-weight `DimReduc` to a Seurat object.
