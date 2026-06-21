# Figure 3. Simulation-based validation of fibroblast state-mixture inference

**A.** Ground-truth fibroblast state composition across continuous-mixture, discrete-state, batch-confounding, and rare-transition simulation scenarios.
**B.** State-weight RMSE for marker scoring, the FibroDynMix initializer, the FibroDynMix NB optimizer, and the FibroDynMix VI posterior mean.
**C.** Dominant-state assignment accuracy across the same simulations.
**D.** Mean absolute error between inferred and true latent entropy, summarizing preservation of mixed-state structure without relying on overplotted entropy scatter.
**E.** Entropy-score distributions for simulator-labeled rare transition-like cells and other cells in the rare-transition scenario.
**F.** Mean state-wise AUPRC for recovery of true marker programs from learned initializer, NB-optimized, and VI-backed `beta_hat`.

All panels are computational benchmark readouts from raw-count generative simulations with known latent truth. They do not claim causal biological transitions in patient tissue.
