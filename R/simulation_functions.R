simulate_score_mse <- function(data_params, estimator_params, seed) {

  dist_args <- list(
    dist_name = data_params$noise_dist,
    params = data_params$noise_params
  )
  
  # Create the true noise functions
  true_quantile_fn <- do.call(create_quantile_function, dist_args)
  true_score_fn <- do.call(create_score_function, dist_args)
  true_density_fn <- do.call(create_density_function, dist_args)
  
  # Define integration grid / pre-calculate true values
  x_grid <- seq(true_quantile_fn(0.001), true_quantile_fn(0.999), length.out = 2000)
  grid_spacing <- diff(x_grid)[1]
  
  true_scores_on_grid <- true_score_fn(x_grid)
  true_density_on_grid <- true_density_fn(x_grid)
  
  # Initialise result list
  
  estimators <- names(estimator_params)
  results_list <- list()
  for (estimator_name in estimators) {
    results_list[[estimator_name]] <- list()
  }
  
  # Main calculation loop
  for (n in data_params$n_values) {
    
    # Generate the data sample ONCE for the current n
    noise_sample <- do.call(sample_from_distribution, list(
      n = n,
      dist_name = data_params$noise_dist,
      params = data_params$noise_params
    ))
    
    estimated_score_fn_list <- score_estimation(noise_sample, estimators)
    
    # INNER loop through estimators for the SAME data sample
    for (estimator_name in estimators) {
      

      estimated_score_fn <- estimated_score_fn_list[[estimator_name]]
      # Calculate density-weighted squared error (approximate integral)
      estimated_scores_on_grid <- estimated_score_fn(x_grid)
      squared_errors <- (true_scores_on_grid - estimated_scores_on_grid)^2
      integral_approx <- sum(squared_errors * true_density_on_grid) * grid_spacing
      
      results_list[[estimator_name]][[as.character(n)]] <- integral_approx
    }
  }
  
  return(unlist(results_list))
}


