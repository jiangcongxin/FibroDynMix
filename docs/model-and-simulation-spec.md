# FibroDynMix Model and Simulation Spec

## Generative Model

For cell `i`, gene `g`, donor `j`, and study `s`, the simulator follows the
initial FibroDynMix count model:

```text
x_ig ~ NegativeBinomial(mu_ig, phi_g)

log(mu_ig) =
    log(library_i)
  + alpha_g
  + sum_k z_ik beta_kg
  + b_sg
  + u_jg
```

where:

- `library_i` is a cell-level size factor.
- `alpha_g` is the gene baseline expression.
- `z_i` is a latent fibroblast state mixture on the K-state simplex.
- `beta_kg` is the state-by-gene expression program.
- `b_sg` is the study-level gene effect.
- `u_jg` is the donor-level gene effect.
- `phi_g` is gene-specific negative-binomial overdispersion.

The latent mixture is generated with a logistic-normal layer:

```text
eta_i ~ Normal(m_j, sigma_state^2 I)
z_i = softmax(eta_i)
```

`m_j` is shifted by donor-level disease metadata so disease donors contain more
inflammatory, myofibroblast, and ECM-remodeling composition.

## Marker Program Generation

The simulator assigns weakly oriented state marker programs:

```text
beta_kg ~ Normal(0, tau_low^2)              otherwise
beta_kg ~ Normal(1.2, (tau_high / 4)^2)     for state-k marker genes
```

This mirrors the planned model principle: prior marker sets orient
identifiability, but state mixtures are generated and later inferred through the
count likelihood rather than fixed gene-set scores.

## Benchmark Scenarios

The first simulator supports four scenarios:

- `continuous`: cells are sampled from smooth logistic-normal state mixtures.
- `discrete`: cells are collapsed to one-hot dominant states.
- `batch_confounding`: disease is confounded with the last study, and selected
  disease-associated marker genes receive an additional study effect.
- `rare_transition`: a small fraction of cells are forced into a resident /
  inflammatory / myofibroblast mixed profile.

These scenarios are intended to benchmark FibroDynMix against AddModuleScore,
AUCell, ssGSEA, NMF, LDA/topic models, and latent-space classifier baselines.

## Initial Evaluation Targets

The first benchmark layer should measure:

- state weight RMSE against the true `z_i`.
- dominant-state accuracy.
- marker recovery AUPRC against true state marker assignments.
- calibration of posterior or bootstrap state uncertainty.
- robustness under study/disease confounding.
- rare transition detection.

## Raw-Count Objective Layer

The package exposes the first raw-count likelihood layer:

```text
fibrodynmix_nb_loglik()
fibrodynmix_nb_deviance()
fibrodynmix_nb_objective()
```

These functions evaluate:

```text
x_ig ~ NegativeBinomial(mu_ig, phi_g)

log(mu_ig) =
    log(library_i)
  + alpha_g
  + sum_k z_ik beta_kg
  + optional study/donor effects
```

This layer is the bridge from the current log-normalized initializer toward a
proper raw-count optimization or variational-inference implementation. The
current optimizer estimates ridge-penalized study- and donor-by-gene effects
when identifiers are supplied.

The first optimizer is:

```text
fit_fibrodynmix_nb()
```

It uses alternating updates:

```text
initialize z, beta from fit_fibrodynmix_initializer()
repeat:
  update alpha_g and beta_kg gene-wise under NB likelihood
  update each z_i on the simplex through softmax logits
  optionally update phi_g by method-of-moments
```

The stabilized optimizer also tracks:

```text
objective trace
best objective and best iteration
rollback flags
early-stopping reason
marker-orientation penalty
```

When `fibrodynmix_nb` is used inside `run_simulation_benchmark()`, the benchmark
table records:

```text
nb_initial_objective
nb_final_objective
nb_best_objective
nb_objective_improvement
nb_best_iteration
nb_executed_iterations
nb_stop_reason
nb_any_rollback
nb_rollback_count
vi_best_elbo
vi_interval_coverage
vi_mean_interval_width
vi_calibrated_interval_coverage
vi_calibrated_mean_interval_width
```

These fields can be summarized with:

```text
summarize_optimizer_diagnostics()
```

Study-effect penalty sensitivity is handled by:

```text
run_study_effect_sensitivity()
```

It benchmarks `study_l2` and `marker_l2` grids under the batch-confounding
scenario and reports:

```text
RMSE
dominant accuracy
NB objective
rollback diagnostics
study_effect_l2_norm
study_effect_mean_abs
```

The default penalty pair can be selected with:

```text
select_study_effect_penalty()
```

The selector ranks candidate penalties by a weighted tradeoff among state RMSE,
NB objective, study-effect magnitude, and rollback frequency, while preferring
solutions within configurable tolerance bands around the best RMSE and
objective.

Current scope:

- includes raw counts, library offset, latent simplex, gene programs, and NB
  overdispersion;
- includes first ridge-penalized study- and donor-effect optimizers through
  `fit_study_effect = TRUE` and `fit_donor_effect = TRUE`;
- is an optimization bridge, not the final hierarchical Bayesian VI model.

## Bootstrap Uncertainty Layer

The first uncertainty implementation is:

```text
bootstrap_fibrodynmix()
summarize_bootstrap_uncertainty()
```

It provides pragmatic uncertainty summaries before full Bayesian VI is
implemented:

```text
cell-level z draws and intervals
cell entropy uncertainty
sample-level composition intervals
marker-program stability from beta bootstrap draws
```

This bootstrap layer is complemented by `fit_fibrodynmix_vi()`, a lightweight
logistic-normal variational posterior over cell state logits around the fitted
NB mode. It returns posterior state draws, credible intervals, and an
ELBO-like trace. In simulation, `evaluate_posterior_intervals()` and
`calibrate_posterior_interval_scale()` quantify raw interval coverage and the
width cost needed to reach target coverage. This is not yet a full amortized or
fully hierarchical variational posterior over all model parameters.

## Transition Flow and Plasticity

The first dynamic layer is:

```text
compute_state_cost()
estimate_transition_flow()
compute_fpi()
```

It estimates condition-to-condition state flow with entropy-regularized optimal
transport:

```text
source composition -> target composition
cost = transcriptional distance between state programs
flow = Sinkhorn-scaled transport matrix
```

The Fibroblast Plasticity Index is currently:

```text
FPI_i = normalized entropy(z_i) + lambda * transition_potential_i
```

This is a cross-sectional transition-flow approximation. It should be described
as inferred state flow, not observed lineage trajectory.
