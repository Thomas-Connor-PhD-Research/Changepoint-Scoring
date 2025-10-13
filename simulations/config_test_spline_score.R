#' Config file to for the "calculate_score_mse" simulation function
sim_params <- list(
  simulation_name = "spline_score_t_dist", # used to label ouput file
  n_reps = 500,
  simulation_function = "calculate_score_mse", 
  # simulation specific parameters
  n_values = c(100, 200, 500, 1000, 2000),
  noise_dist = "t",
  noise_params = list(scale = 1, df = 3)
)