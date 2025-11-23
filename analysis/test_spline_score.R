# source("R/noise_distributions.R")
# 
# plot_single_estimated_versus_true <- function(n, noise_dist, ...) {
#   noise_params <- list(...)
#   noise_sample <- do.call(sample_from_distribution, list(n = n, 
#                                                          dist = noise_dist, 
#                                                          params = noise_params))
#   
#   true_score_fn <- do.call(create_score_function, list(
#     dist = noise_dist, 
#     params = noise_params)
#   )
#   
#   cv_df <- cv_spline_score(noise_sample)$df_min
#   estimated_score_fn <- spline_score(noise_sample, df = cv_df)$rho
#   
#   # 3. Set up a grid and calculate scores for plotting
#   x_grid <- seq(min(noise_sample), max(noise_sample), length.out = 1000)
#   true_scores <- true_score_fn(x_grid)
#   estimated_scores <- estimated_score_fn(x_grid)
#   
#   # 4. Create the main plot
#   plot(x_grid, true_scores, type = "l", col = "blue", lwd = 2,
#        main = "Comparison of True vs. Spline-Estimated Score Function",
#        xlab = "x", ylab = "Score rho(x)",
#        ylim = range(c(true_scores, estimated_scores), na.rm = TRUE))
#   
#   # 5. Add the estimated score and a rug plot for data density
#   lines(x_grid, estimated_scores, col = "red", lty = 2, lwd = 2)
#   
#   # Add ticks to the x-axis to show the density of the sampled points
#   rug(noise_sample, col = scales::alpha("black", 0.3))
#   
#   # 6. Add legend and grid
#   legend("topright", 
#          legend = c("True Score", "Spline-Estimated Score", "Data Points"),
#          col = c("blue", "red", "darkgrey"), 
#          lwd = c(2, 2, 2), 
#          lty = c(1, 2, 1),
#          bty = "n") # bty = "n" removes the box
#   
#   grid()
#   
#   return(invisible(NULL))
# }
# 
# plot_multiple_estimated_versus_true <- function(n, noise_dist, n_reps = 20, ...) {
#   
#   # 1. Capture noise parameters and get the true score function
#   noise_params <- list(...)
#   
#   true_score_fn <- do.call(create_score_function, list(dist_name = noise_dist, 
#                                                        params = noise_params))
#   
#   # Generate one sample just to define a reasonable x-axis grid
#   initial_sample <- do.call(sample_from_distribution, list(n = n,
#                                                            dist_name = noise_dist, 
#                                                            params = noise_params))
#   x_grid <- seq(min(initial_sample), max(initial_sample), length.out = 1000)
#   true_scores <- true_score_fn(x_grid)
#   
#   # 2. Pre-calculate all estimated scores to determine the plot's y-range
#   all_estimated_scores <- list()
#   for (i in 1:n_reps) {
#     noise_sample <- do.call(sample_from_distribution, list(n = n,
#                                                            dist_name = noise_dist, 
#                                                            params = noise_params))
#     
#     cv_df <- cv_spline_score(noise_sample)$df_min
#     estimated_score_fn <- spline_score(noise_sample, df = cv_df)$rho
#     all_estimated_scores[[i]] <- estimated_score_fn(x_grid)
#   }
#   
#   # 3. Create the main plot with the true score function
#   plot(x_grid, true_scores, type = "l", col = "blue", lwd = 3,
#        main = "True Score vs. Multiple Spline-Estimated Scores",
#        xlab = "x", ylab = "Score rho(x)",
#        ylim = range(c(true_scores, unlist(all_estimated_scores)), na.rm = TRUE))
#   
#   # 4. Add all the estimated score functions as semi-transparent grey lines
#   for (i in 1:n_reps) {
#     lines(x_grid, all_estimated_scores[[i]], col = scales::alpha("grey", 0.6), lwd = 1.5)
#   }
#   
#   # Re-plot the true score on top to ensure it's not obscured
#   lines(x_grid, true_scores, col = "blue", lwd = 3)
#   
#   # 5. Add a rug plot for the density of the *first* sample for reference
#   rug(initial_sample, col = scales::alpha("black", 0.3))
#   
#   # 6. Add legend and grid
#   legend("topright", 
#          legend = c("True Score", "Estimated Score Replicates"),
#          col = c("blue", "grey"), 
#          lwd = c(3, 2), 
#          bty = "n")
#   
#   grid()
#   
#   return(invisible(NULL))
# }
# 
# generate_QQ_plot <- function(n, noise_dist, ...) {
#   noise_params <- list(...)
#   
#   noise_sample <- do.call(sample_from_distribution, list(n = n,
#                                                          dist_name = noise_dist, 
#                                                          params = noise_params))
#   
#   cv_df <- cv_spline_score(noise_sample)$df_min
#   estimated_score_fn <- spline_score(noise_sample, df = cv_df)$rho
#   
#   # The dist of rho(eps) under normal model is N(0, 1/sigma^2)
#   if (noise_dist == "normal"){
#     sigma <- noise_params$sd
#     scaled_estimated_score <- estimated_score_fn(noise_sample) * sigma
#     plot_title <- paste("Q-Q Plot for", noise_dist, "Noise (n =", n, ", df =", round(cv_df, 2), ")")
#     qqnorm(scaled_estimated_score,
#            main = plot_title,
#            xlab = "Theoretical N(0,1) Quantiles",
#            ylab = "Sample Quantiles of Scaled Score")
#     qqline(scaled_estimated_score, col = "steelblue", lwd = 2)
#     
#   }
# }
# 
# plot_single_estimated_versus_true(1000, "normal")

source("R/noise_distributions.R")
# We assume the `score_estimation` function from the previous prompt 
# is also sourced and available.

plot_single_estimated_versus_true <- function(n, Y = NULL, noise_dist, noise_params,
                                              scoring_type = "spline_df_min") {
  
  # 1. Generate data
  if (is.null(Y)){
    noise_sample <- do.call(sample_from_distribution, list(n = n, 
                                                           dist_name = noise_dist, 
                                                           params = noise_params))
  }else{
     noise_sample <- Y
  }
  
  
  # 2. Get true score function
  true_score_fn <- do.call(create_score_function, list(
    dist_name = noise_dist, 
    params = noise_params)
  )
  
  # 3. Get estimated score function using the new generalized function
  # This returns a named list, e.g., list(spline_df_min = <fn>)
  estimated_fn_list <- score_estimation(noise_sample, scoring_type)
  estimated_score_fn <- estimated_fn_list[[1]] # Extract the function
  
  # 4. Set up a grid and calculate scores for plotting
  x_grid <- seq(min(noise_sample), max(noise_sample), length.out = 1000)
  true_scores <- true_score_fn(x_grid)
  estimated_scores <- estimated_score_fn(x_grid)
  
  # 5. Create the main plot
  plot_title <- paste("True vs.", scoring_type, "Estimated Score")
  plot(x_grid, true_scores, type = "l", col = "blue", lwd = 2,
       main = plot_title,
       xlab = "x", ylab = "Score rho(x)",
       ylim = range(c(true_scores, estimated_scores), na.rm = TRUE))
  
  # 6. Add the estimated score and a rug plot
  lines(x_grid, estimated_scores, col = "red", lty = 2, lwd = 2)
  rug(noise_sample, col = scales::alpha("black", 0.3))
  
  # 7. Add legend and grid
  legend_text <- c("True Score", paste(scoring_type, "Estimate"), "Data Points")
  legend("topright", 
         legend = legend_text,
         col = c("blue", "red", "darkgrey"), 
         lwd = c(2, 2, 2), 
         lty = c(1, 2, 1),
         bty = "n") # bty = "n" removes the box
  
  grid()
  
  return(invisible(NULL))
}

plot_multiple_estimated_versus_true <- function(n, noise_dist, scoring_type = "spline_df_min", n_reps = 20, ...) {
  
  # 1. Capture noise parameters and get the true score function
  noise_params <- list(...)
  
  true_score_fn <- do.call(create_score_function, list(dist_name = noise_dist, 
                                                       params = noise_params))
  
  # Generate one sample just to define a reasonable x-axis grid
  initial_sample <- do.call(sample_from_distribution, list(n = n,
                                                           dist_name = noise_dist, 
                                                           params = noise_params,
                                                           seed = NULL))
  x_grid <- seq(min(initial_sample), max(initial_sample), length.out = 1000)
  true_scores <- true_score_fn(x_grid)
  
  # 2. Pre-calculate all estimated scores
  all_estimated_scores <- list()
  for (i in 1:n_reps) {
    noise_sample <- do.call(sample_from_distribution, list(n = n,
                                                           dist_name = noise_dist, 
                                                           params = noise_params))
    
    # Use the generalized score_estimation function
    estimated_fn_list <- score_estimation(noise_sample, scoring_type)
    estimated_score_fn <- estimated_fn_list[[1]] # Extract the function
    
    all_estimated_scores[[i]] <- estimated_score_fn(x_grid)
  }
  
  # 3. Create the main plot with the true score function
  plot_title <- paste("True Score vs. Multiple", scoring_type, "Estimates")
  plot(x_grid, true_scores, type = "l", col = "blue", lwd = 3,
       main = plot_title,
       xlab = "x", ylab = "Score rho(x)",
       ylim = range(c(true_scores, unlist(all_estimated_scores)), na.rm = TRUE))
  
  # 4. Add all the estimated score functions
  for (i in 1:n_reps) {
    lines(x_grid, all_estimated_scores[[i]], col = scales::alpha("grey", 0.6), lwd = 1.5)
  }
  
  # Re-plot the true score on top
  lines(x_grid, true_scores, col = "blue", lwd = 3)
  
  # 5. Add a rug plot for the density
  rug(initial_sample, col = scales::alpha("black", 0.3))
  
  # 6. Add legend and grid
  legend_text <- c("True Score", paste(scoring_type, "Replicates"))
  legend("topright", 
         legend = legend_text,
         col = c("blue", "grey"), 
         lwd = c(3, 2), 
         bty = "n")
  
  grid()
  
  return(invisible(NULL))
}

generate_QQ_plot <- function(n, noise_dist, scoring_type = "spline_df_min", ...) {
  noise_params <- list(...)
  
  noise_sample <- do.call(sample_from_distribution, list(n = n,
                                                         dist_name = noise_dist, 
                                                         params = noise_params))
  
  # Use the generalized score_estimation function
  estimated_fn_list <- score_estimation(noise_sample, scoring_type)
  estimated_score_fn <- estimated_fn_list[[1]] # Extract the function
  
  # The dist of rho(eps) under normal model is N(0, 1/sigma^2)
  if (noise_dist == "normal"){
    sigma <- noise_params$sd
    scaled_estimated_score <- estimated_score_fn(noise_sample) * sigma
    
    # Updated plot title to be dynamic
    plot_title <- paste("Q-Q Plot for", scoring_type, "on", noise_dist, "Noise (n =", n, ")")
    
    qqnorm(scaled_estimated_score,
           main = plot_title,
           xlab = "Theoretical N(0,1) Quantiles",
           ylab = "Sample Quantiles of Scaled Score")
    qqline(scaled_estimated_score, col = "steelblue", lwd = 2)
  }
}

