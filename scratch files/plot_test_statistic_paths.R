# This scratch script generates a single realization of a changepoint problem
# and plots the entire test statistic path for CUSUM and a scoring-based
# estimator to visualize their behavior.
source("R/noise_distributions.R")
source("R/score_estimation2.R")

calculate_cusum_path <- function(Y, params) {
  burnin <- params$burnin
  n <- params$n
  
  tau_grid <- (burnin + 1):(n - burnin)
  
  denom <- (tau_grid) * (1 - (tau_grid) / n)
  unscaled_cusum <- cumsum(Y - mean(Y))[tau_grid]
  statistic_values <- unscaled_cusum^2 / denom
  
  # Return a data frame of the path
  return(data.frame(tau = tau_grid, statistic_value = statistic_values))
}

#' Calculates the full scoring test statistic path.
calculate_scoring_path <- function(Y, params) {
  burnin <- params$burnin
  n <- params$n
  
  # a. Get an initial CUSUM estimate to calculate residuals
  cusum_path <- calculate_cusum_path(Y, params)
  tau_hat_initial <- cusum_path$tau[which.max(cusum_path$statistic_value)]
  
  # b. Calculate residuals based on the initial estimate
  residuals <- Y
  mean_pre <- mean(Y[1:tau_hat_initial])
  mean_post <- mean(Y[(tau_hat_initial + 1):n])
  residuals[(tau_hat_initial + 1):n] <- residuals[(tau_hat_initial + 1):n] - (mean_post - mean_pre)
  
  # c. Learn the score function from the residuals
  score_fn <- spline_score(residuals, df = cv_spline_score(residuals)$df_min)$rho
  Y_scored <- score_fn(Y)
  
  # d. Calculate the test statistic on the SCORED data
  tau_grid <- (burnin + 1):(n - burnin)
  denom <- (tau_grid) * (1 - (tau_grid) / n)
  unscaled_scoring <- cumsum(Y_scored - mean(Y_scored))[tau_grid]
  statistic_values <- unscaled_scoring^2 / denom
  
  # Return a data frame of the path
  return(data.frame(tau = tau_grid, statistic_value = statistic_values))
}


# --- 3. Set Simulation Parameters (for a single run) ---
sim_params <- list(
  n = 1000,      # Sample size
  burnin = 1,
  changepoint_spec = list(
    tau = 500,   # True changepoint location
    delta = 0.8
  ),
  
  # Noise Specification (Gaussian)
  noise_dist = "normal",
  noise_params = list(sd = 1)
)


# --- 4. Generate Data and Compute Paths --- # for reproducibility

x_limits <- c(1, 999)

# a. Generate a single dataset with a changepoint
noise <- do.call(sample_from_distribution, list(n = sim_params$n, 
                                             dist = sim_params$noise_dist, 
                                             params = sim_params$noise_params))
Y <- noise
Y[(sim_params$changepoint_spec$tau + 1):sim_params$n] <- Y[(sim_params$changepoint_spec$tau + 1):sim_params$n] + sim_params$changepoint_spec$delta

# b. Calculate the full path for each estimator
cusum_path_df <- calculate_cusum_path(Y, sim_params)
scoring_path_df <- calculate_scoring_path(Y, sim_params)

tau_hat_cusum <- cusum_path_df$tau[which.max(cusum_path_df$statistic_value)]
tau_hat_scoring <- scoring_path_df$tau[which.max(scoring_path_df$statistic_value)]


# c. Combine into a single data frame for plotting
cusum_path_df$estimator <- "CUSUM"
scoring_path_df$estimator <- "Scoring"
plot_data <- rbind(cusum_path_df, scoring_path_df)


# --- 5. Plot the Results ---
print(
  ggplot(plot_data, aes(x = tau, y = statistic_value, color = estimator)) +
    geom_line(linewidth = 1.1) +
    
    # Add vertical lines for the TRUE and ESTIMATED changepoints
    geom_vline(xintercept = sim_params$changepoint_spec$tau, linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = tau_hat_cusum, linetype = "dotted", color = "dodgerblue", linewidth = 1) +
    geom_vline(xintercept = tau_hat_scoring, linetype = "dotted", color = "firebrick", linewidth = 1) +
    
    # Add labels for the changepoints
    annotate("text", x = sim_params$changepoint_spec$tau, y = 0, label = paste("True Tau =", sim_params$changepoint_spec$tau), hjust = 1.1, vjust = -0.5) +
    annotate("text", x = tau_hat_cusum, y = 0, label = paste("CUSUM Tau =", tau_hat_cusum), hjust = -0.1, vjust = -2.5, color = "dodgerblue") +
    annotate("text", x = tau_hat_scoring, y = 0, label = paste("Scoring Tau =", tau_hat_scoring), hjust = -0.1, vjust = -0.5, color = "firebrick") +
    
    scale_color_manual(values = c("CUSUM" = "dodgerblue", "Scoring" = "firebrick")) +
    
    # NEW: Add coordinate limits if they are defined
    coord_cartesian(xlim = x_limits) +
    
    labs(
      title = "Test Statistic Path vs. Potential Changepoint",
      subtitle = paste("Gaussian Noise, delta =", sim_params$changepoint_spec$delta),
      x = "Potential Changepoint (tau)",
      y = "Test Statistic Value",
      color = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
)


true_score_fn <- do.call(create_score_function, list(
  dist = sim_params$noise_dist, 
  params = noise_params)
)

cv_df <- cv_spline_score(noise_sample)$df_min
estimated_score_fn <- spline_score(noise_sample, df = cv_df)$rho

# 3. Set up a grid and calculate scores for plotting
x_grid <- seq(min(noise_sample), max(noise_sample), length.out = 1000)
true_scores <- true_score_fn(x_grid)
estimated_scores <- estimated_score_fn(x_grid)

# 4. Create the main plot
plot(x_grid, true_scores, type = "l", col = "blue", lwd = 2,
     main = "Comparison of True vs. Spline-Estimated Score Function",
     xlab = "x", ylab = "Score rho(x)",
     ylim = range(c(true_scores, estimated_scores), na.rm = TRUE))

# 5. Add the estimated score and a rug plot for data density
lines(x_grid, estimated_scores, col = "red", lty = 2, lwd = 2)

# Add ticks to the x-axis to show the density of the sampled points
rug(noise_sample, col = scales::alpha("black", 0.3))

# 6. Add legend and grid
legend("topright", 
       legend = c("True Score", "Spline-Estimated Score", "Data Points"),
       col = c("blue", "red", "darkgrey"), 
       lwd = c(2, 2, 2), 
       lty = c(1, 2, 1),
       bty = "n") # bty = "n" removes the box

grid()

return(invisible(NULL))

