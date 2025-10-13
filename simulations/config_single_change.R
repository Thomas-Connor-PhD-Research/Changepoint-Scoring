sim_params <- list(
  simulation_name = "single_change_cauchy_dist", # used to label ouput file
  n_reps = 1000,

  # Simulation type to run
  simulation_function = "simulate_single_changepoint", 
  
  # Simulation specific parameters
  n = 1000,
  burnin = 40, # typically floor(n/25)
  changepoint_spec = list(
    tau = 500, # typically floor(n/2)
    delta = 1
  ),

  noise_dist = "cauchy",
  noise_params = list(scale = 1),
  # Estimators to compute
  # To add an estimator: new_estimator = list() and write a calculate_new_estim-
  # ator function in R/simulation_functions.R
  # Ensure estimators are listed AFTER their dependencies (eg: scoring with 
  # use_cusum_tau = TRUE after cusum)
  estimators = list("cusum", "oracle", "scoring", "scoring_known_noise", "iterative_scoring") 
)