# This scratch script generates multiple realizations of a changepoint problem
# and plots all of the resulting scoring statistic paths on a single graph
# to visualize the estimator's variability.

source("R/noise_distributions.R")


calculate_cusum_path <- function(Y, params) {
  burnin <- params$burnin
  n <- params$n
  tau_grid <- (burnin + 1):(n - burnin)
  
  denom <- (tau_grid) * (1 - (tau_grid) / n)
  unscaled_cusum <- cumsum(Y - mean(Y))[tau_grid]
  statistic_values <- unscaled_cusum^2 / denom
  
  return(data.frame(tau = tau_grid, statistic_value = statistic_values))
}

calculate_scoring_path <- function(Y, params) {
  burnin <- params$burnin
  n <- params$n
  
  cusum_path <- calculate_cusum_path(Y, params)
  tau_hat_initial <- cusum_path$tau[which.max(cusum_path$statistic_value)]
  
  residuals <- Y
  mean_pre <- mean(Y[1:tau_hat_initial])
  mean_post <- mean(Y[(tau_hat_initial + 1):n])
  residuals[(tau_hat_initial + 1):n] <- residuals[(tau_hat_initial + 1):n] - (mean_post - mean_pre)
  
  score_fn <- spline_score(residuals, df = cv_spline_score(residuals)$df_min)$rho
  Y_scored <- score_fn(Y)
  
  tau_grid <- (burnin + 1):(n - burnin)
  denom <- (tau_grid) * (1 - (tau_grid) / n)
  unscaled_scoring <- cumsum(Y_scored - mean(Y_scored))[tau_grid]
  statistic_values <- unscaled_scoring^2 / denom
  
  return(data.frame(tau = tau_grid, statistic_value = statistic_values))
}


# --- 3. Set Simulation Parameters ---
sim_params <- list(
  n_reps = 50, # NEW: Number of realizations to plot
  n = 1000,
  burnin = 40,
  changepoint_spec = list(tau = 500, delta = 0.8),
  noise_dist = "normal",
  noise_params = list(sd = 1)
)


# --- 4. Generate Data and Compute Paths for Multiple Realizations ---
set.seed(123) # for reproducibility
all_paths <- list() # A list to store the data frame for each path

for (i in 1:sim_params$n_reps) {
  # a. Generate a new dataset for each replication
  noise <-  do.call(sample_from_distribution, list(n = sim_params$n, 
                                                   dist_name = sim_params$noise_dist, 
                                                   params = sim_params$noise_params))
  Y <- noise
  Y[(sim_params$changepoint_spec$tau + 1):sim_params$n] <- Y[(sim_params$changepoint_spec$tau + 1):sim_params$n] + sim_params$changepoint_spec$delta
  
  # b. Calculate the scoring path for this realization
  path_df <- calculate_scoring_path(Y, sim_params)
  path_df$replication <- i # Add an ID for this specific run
  
  all_paths[[i]] <- path_df
}

# Combine all the individual paths into one large data frame
plot_data <- bind_rows(all_paths)


# --- 5. Plot the Results ---
# We wrap the ggplot object in print() to ensure it renders when run as a script.
print(
  ggplot(plot_data, aes(x = tau, y = statistic_value, group = replication)) +
    # Plot all the individual paths in semi-transparent grey
    geom_line(color = "grey80") +
    
    # Add a vertical line for the true changepoint
    geom_vline(xintercept = sim_params$changepoint_spec$tau, 
               linetype = "dashed", 
               color = "black", 
               linewidth = 1) +
    
    # Add a text label for the true changepoint
    annotate("text", 
             x = sim_params$changepoint_spec$tau, 
             y = max(plot_data$statistic_value) * 0.9, 
             label = paste("True Tau =", sim_params$changepoint_spec$tau), 
             hjust = 1.1, 
             color = "black") +
    
    labs(
      title = "Multiple Realizations of the Scoring Test Statistic Path",
      subtitle = paste(sim_params$n_reps, "replications with", sim_params$noise_dist, "noise"),
      x = "Potential Changepoint (tau)",
      y = "Test Statistic Value"
    ) +
    theme_minimal()
)