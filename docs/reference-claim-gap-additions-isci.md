# Reference Claim-Gap Additions for iScience Draft

This note records the new references added during claim-gap supplementation.
Each source is tied to a bounded manuscript claim so the reference list does not
become a general background dump.

| Claim gap | Added source | Evidence role | Manuscript-safe use |
|---|---|---|---|
| Fibroblast states recur across tissues and diseases, including inflammatory and fibrotic settings | Korsunsky et al., Med 2022 | Cross-tissue single-cell stromal atlas across inflammatory diseases | Supports shared pathological fibroblast activation states; does not validate FibroDynMix |
| Human fibrosis contains heterogeneous scar-forming fibroblast/myofibroblast populations | Kuppe et al., Nature 2021 | Human kidney fibrosis single-cell and spatial analyses | Supports fibroblast heterogeneity and myofibroblast origin context; not a scar-module validation source |
| Wound fibroblasts can occupy transition-like and plastic states | Guerrero-Juarez et al., Nat. Commun. 2019 | Mouse wound single-cell, pseudotime, RNA velocity, lineage tracing | Supports plasticity and differentiation-state language in wound repair |
| Fibrotic lung disease includes multiple fibroblast and myofibroblast populations | Valenzi et al., Ann. Rheum. Dis. 2019 | Human SSc-ILD single-cell analysis | Supports disease-associated fibroblast heterogeneity |
| Single cells may contain mixtures of identity and activity gene programs | Kotliar et al., eLife 2019 | cNMF benchmark and single-cell gene-program usage analysis | Supports mixture-program and topic-model baseline framing |
| Cross-dataset correction is commonly handled by aligning cells in expression or embedding spaces | Korsunsky et al., Nat. Methods 2019; Haghverdi et al., Nat. Biotechnol. 2018 | Harmony and MNN single-cell integration methods | Supports contrast with FibroDynMix study/donor terms; does not imply identical modeling target |
| Single-cell expression snapshots can be linked with optimal transport for trajectory-style analysis | Schiebinger et al., Cell 2019 | Waddington-OT applied to reprogramming time courses | Supports OT precedent; does not make FPI a temporal lineage measure |
| Keloid/scar ECM can show disease-associated matrix remodeling and mis-differentiation signals | Barallobre-Barreiro et al., Matrix Biology Plus 2019 | Proteomic comparison of normal skin, normal scar, and keloid ECM | Supports scar/fibrosis biological plausibility; not diagnostic validation |
| FibroDynMix state weights live on a simplex and should be interpreted as compositional quantities | Aitchison and Shen, Biometrika 1980; Aitchison, J. R. Stat. Soc. Series B 1982 | Logistic-normal and compositional-data foundations | Supports simplex/logistic-normal state-weight wording; does not validate fibroblast-state biology |
| Bootstrap intervals summarize refitting-based uncertainty rather than posterior probability | Efron, Ann. Stat. 1979 | Bootstrap resampling foundation | Supports the resampling uncertainty procedure; not evidence that intervals are calibrated in every dataset |
| The lightweight VI posterior is an approximation over fitted state logits | Blei et al., J. Am. Stat. Assoc. 2017 | Variational-inference review | Supports approximation language; does not imply exact Bayesian posterior inference |
| Raw UMI count modeling is preferable to treating normalized expression as the primary response for model fitting | Townes et al., Genome Biology 2019 | Single-cell count-modeling and dimension-reduction method | Supports direct count-response framing; does not independently validate the FibroDynMix NB likelihood |
| Small-sample scar-module checks should preserve sample-level grouping and be interpreted cautiously | Soneson and Robinson, Nat. Methods 2018; Crowell et al., Nat. Commun. 2020 | Single-cell differential-expression benchmarking and multi-sample differential-state modeling | Supports pseudobulk/sample-level boundary language; not proof of significant scar-module differences |
| STAR Methods should document reproducible single-cell workflows and analysis conventions | Luecken and Theis, Mol. Syst. Biol. 2019 | Single-cell analysis best-practice tutorial | Supports methods-convention wording; does not add new empirical validation |
