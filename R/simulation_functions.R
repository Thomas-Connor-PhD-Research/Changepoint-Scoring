simulate_score_mse <- function(data_params, estimator_params) {

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
  results_list <- list()
  for (estimator_name in names(estimator_params)) {
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
    
    # INNER loop through estimators for the SAME data sample
    for (estimator_name in names(estimator_params)) {
      

      if (estimator_name %in% c("spline_df_min", "spline_df_1se", "asm")) {
        estimated_score_fn <- score_estimation(noise_sample, estimator_name)
      } else {
        stop(paste("Unknown estimator specified:", estimator_name))
      } # Change to a do.call in future if necessary
      
      # Calculate density-weighted squared error (approximate integral)
      estimated_scores_on_grid <- estimated_score_fn(x_grid)
      squared_errors <- (true_scores_on_grid - estimated_scores_on_grid)^2
      integral_approx <- sum(squared_errors * true_density_on_grid) * grid_spacing
      
      results_list[[estimator_name]][[as.character(n)]] <- integral_approx
    }
  }
  
  return(unlist(results_list))
}


