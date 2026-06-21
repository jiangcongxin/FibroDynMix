# Figure 4. Study-effect penalty calibration under batch confounding

**A.** Mean state-weight RMSE across study-effect ridge penalties and marker-orientation penalties. Dashed line indicates the no-study-effect FibroDynMix NB baseline.
**B.** Mean raw-count NB objective across the same penalty grid. Dashed line indicates the no-study-effect baseline.
**C.** Mean L2 norm of the fitted study-by-gene effect, showing stronger shrinkage at larger `study_l2`.
**D.** Weighted tradeoff score used to select the default penalty pair. Circled point marks the selected setting.

Selected study_l2=5 and marker_l2=0.05 using weighted tradeoff score 0.111; RMSE 0.1892 is within the 5.0% tolerance around best RMSE 0.1892, and NB objective 1.2691 is within the 3.0% tolerance around best objective 1.2648.

All panels are simulation-calibration readouts under batch confounding and do not claim biological validation in patient tissue or a universal hyperparameter setting.
