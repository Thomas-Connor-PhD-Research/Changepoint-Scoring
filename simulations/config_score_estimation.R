#' Config file to for the "calculate_score_mse" simulation function
# sim_params <- list(
#   simulation_name = "score_estimation_normal_dist", # used to label ouput file
#   n_reps = 100,
#   simulation_function = "calculate_score_mse", 
#   # simulation specific parameters
#   
#   estimators = c("asm", "spline_df_min", "spline_df_1se")
# )

simulation_params <- list(
  simulation_name = "score_est_t_multiple_se",
  n_reps = 10,
  simulation_function = "simulate_score_mse",
  seed = 123
)

data_params <- list(
  n_values = c(250, 500, 1000, 2000),
  noise_dist = "t",
  noise_params = list(scale = 0.5, df = 3)
)

estimator_params <- list(
  "asm" = list(),
  "spline_df_min" = list(),
  "spline_df_1se" = list(),
  "spline_df_2se" = list(),
  "spline_df_3se" = list()
)