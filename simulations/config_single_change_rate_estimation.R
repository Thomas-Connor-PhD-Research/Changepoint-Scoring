# Simulation Control Parameters
simulation_params <- list(
  simulation_name = "compare_rate_estimate_t_oracle_tests",
  n_reps = 20000,
  simulation_function = "simulate_single_changepoint_rate_estimation",
  seed = 100
)

# Data Generation Parameters
data_params <- list(
  n = c(200, 400, 800, 1600, 3200, 6400, 12800),
  changepoint_spec = list(
    lambda = 0.5,
    delta = function(n){return(3/n^(0.25))}
  ),
  noise_dist = "t",
  noise_params = list(df=3, scale = 1),
  noise_mu = 0
)

# Estimator Parameters
# A named list where each estimator has its own parameter list.
estimator_params <- list(
  
  # # # CUSUM takes no parameters
  "cusum" = list(estimator_type = "cusum"),
  
  "oracle" = list(estimator_type = "oracle"),

  "score_test_fisher" = list(estimator_type = "score_test",
                             use_fisher_update = TRUE),
  
  "bayes" = list(estimator_type = "oracle_bayes",
                 use_fisher_update = TRUE)
  
  # 
  # "score_test_median" = list(estimator_type = "score_test",
  #                            use_fisher_update = FALSE)
  
)
  
  # "cusum_chosen_variance_0" = list(estimator_type = "cusum_chosen_variance",
  #                                    alpha = 0),
  # 
  # "cusum_chosen_variance_0-25" = list(estimator_type = "cusum_chosen_variance",
  #                                  alpha = 0.25),
  # 
  # "cusum_chosen_variance_0-5" = list(estimator_type = "cusum_chosen_variance",
  #                                     alpha = 0.5),
  # 
  # "cusum_chosen_variance_0-75" = list(estimator_type = "cusum_chosen_variance",
  #                                     alpha = 0.75)
  # 
  # 
  # # Oracle takes no parameters (it gets them from data_params)
  # "oracle" = list(estimator_type = "oracle"),

  # "sign" = list(estimator_type = "sign_of_medians"),
  # 
  # "oracle_sign" = list(estimator_type = "oracle_sign"),
  # 
  # "cusum_no_scaling" = list(estimator_type = "oracle_trimmed_cusum",
  #                                  trim = 0.00), # cusum no scaling
  # 
  # "oracle_trimmed_cusum_2%" = list(estimator_type = "oracle_trimmed_cusum",
  #                                  trim = 0.02),
  # 
  # "oracle_trimmed_cusum_10%" = list(estimator_type = "oracle_trimmed_cusum",
  #                                  trim = 0.10),
  # 
  # "oracle_trimmed_cusum_50%" = list(estimator_type = "oracle_trimmed_cusum",
  #                                  trim = 0.50)
  # 
  
  
  
  # "ECDF" = list()

# Hinkley [65, 66]

# Non-parametric estimators:
# Rank-based statistics
# Mann-Whitney statistic
# Kolomorgogorv-Smirnov statistics