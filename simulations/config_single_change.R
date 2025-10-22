# Simulation Control Parameters
simulation_params <- list(
  simulation_name = "single_change_t_dist_refactored",
  n_reps = 10,
  # This simulation function name is now just a label,
  # the R script directly calls simulate_single_changepoint
  simulation_function = "simulate_single_changepoint" 
)

# Data Generation Parameters
data_params <- list(
  n = 1000,
  burnin = 40, # Note: 'burnin' isn't used in your generate_timeseries
  changepoint_spec = list(
    tau = 500,
    delta = 0.3
  ),
  noise_dist = "t",
  noise_params = list(scale = 0.5, df = 3)
)

# Estimator Parameters
# A named list where each estimator has its own parameter list.
estimator_params <- list(
  
  # CUSUM takes no parameters
  # "cusum" = list(), 
  
  # Oracle takes no parameters (it gets them from data_params)
  "oracle" = list(), 
  
  # Scoring (1-iteration)
  "scoring" = list(
    scoring_method = "spline_df_1se", 
    iterations = 1
  ),
  
  # Iterative Scoring (3-iterations)
  "iterative_scoring" = list(
    scoring_method = "spline_df_min",
    iterations = 3
  ),
  
  # Scoring with known noise
  "scoring_known_noise" = list(
    scoring_method = "spline_df_1se" # Spline-based estimation
  )
)

