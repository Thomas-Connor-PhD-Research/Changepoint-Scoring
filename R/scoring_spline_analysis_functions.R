#' @description Given a set of params with n_values, noise_dist, noise_params,
#' calculates density weighted MSE for estimated score v true score at each n in 
#' n_values

calculate_score_mse <- function(params) {
  results_list <- list()
  
  true_quantile_fn <- do.call(create_quantile_function,list(dist_name = params$noise_dist, 
                                                                 params =  params$noise_params))
  
  x_grid <- seq(true_quantile_fn(0.001), true_quantile_fn(0.999), length.out = 2000)
  
  grid_spacing <- diff(x_grid)[1]
  
  true_score_fn <- do.call(create_score_function, list(dist_name = params$noise_dist, 
                                                       params =  params$noise_params))
  
  true_density_fn <- do.call(create_density_function, list(dist_name = params$noise_dist, 
                                                           params =  params$noise_params))
  
  true_scores_on_grid <- true_score_fn(x_grid)
  true_density_on_grid <- true_density_fn(x_grid)
  
  for (n in params$n_values) {
    
      noise_sample <- do.call(sample_from_distribution, list(n=n, 
                                                             dist_name = params$noise_dist, 
                                                             params =  params$noise_params))
                              
      estimated_score_fn <- spline_score(noise_sample, df = cv_spline_score(noise_sample)$df_min)$rho
      
      # 5. Calculate the density-weighted squared error via numerical integration
      estimated_scores_on_grid <- estimated_score_fn(x_grid)
      
      # Calculate the squared error at each grid point
      squared_errors <- (true_scores_on_grid - estimated_scores_on_grid)^2
      
      integral_approx <- sum(squared_errors * true_density_on_grid) * grid_spacing
      
      results_list[[as.character(n)]] <- integral_approx
    }
    
  return(results_list)
}


  