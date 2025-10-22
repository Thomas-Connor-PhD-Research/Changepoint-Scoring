# This script runs a self-contained simulation to compare a naive vs. a
# score-based delta estimator for a single changepoint in Cauchy noise.

# --- 1. Load All Necessary Functions and Packages ---
# This ensures the script can run independently.
source("R/simulation_functions.R")
source("R/noise_distributions.R")
source("R/score_estimation.R")

library(foreach)
library(progressr)
library(future)
library(doFuture)
library(doRNG)


# --- 2. Define a NEW Score-Based Delta Estimator ---
# This function will be used in our simulation. It estimates delta by
# taking the difference in the mean of the *scored* data.

calculate_scoringdelta_estimator <- function(Y, params, all_results) {
  
  # This estimator depends on a prior CUSUM run to get tau
  if (!"cusum" %in% names(all_results)) {
    stop("The 'scoringdelta' estimator requires 'cusum' to be run first.", call. = FALSE)
  }
  tau_hat <- all_results$cusum$cusum.tau
  
  # a. Calculate residuals based on the CUSUM tau estimate
  residuals <- Y
  mean_pre <- mean(Y[1:tau_hat])
  mean_post <- mean(Y[(tau_hat + 1):params$n])
  residuals[(tau_hat + 1):params$n] <- residuals[(tau_hat + 1):params$n] - (mean_post - mean_pre)
  
  # b. Learn the score function from the residuals
  score_fn <- spline_score(residuals, df = cv_spline_score(residuals)$df_min)$rho
  Y_scored <- score_fn(Y)
  
  delta_hat_scoring <- mean(Y_scored[(tau_hat + 1):params$n]) - mean(Y_scored[1:tau_hat])
  delta_hat <- mean()
  # We only care about delta for this estimator
  results <- list(scoringnaive.delta = delta_hat_scoring,
                  scoringclever.delta = )
  names(results) <- paste0("scoringdelta.", names(results))
  
  return(results)
}


# --- 3. Define a NEW Simulation Function ---
# This function is a modified version of your main simulation function.
# It's tailored to calculate only the necessary delta estimates.

simulate_delta_comparison <- function(params) {
  # a. Generate data with a single changepoint
  noise <- sample_from_distribution(params$n, params$noise_dist, !!!params$noise_params)
  Y <- noise
  Y[(params$changepoint_spec$tau + 1):params$n] <- Y[(params$changepoint_spec$tau + 1):params$n] + params$changepoint_spec$delta
  
  # b. Run the estimators sequentially
  all_results <- list()
  
  # The CUSUM estimator provides the naive delta and the tau for the scoring method
  all_results[["cusum"]] <- calculate_cusum_estimator(Y, params)
  
  # The new scoring-based delta estimator
  all_results[["scoringdelta"]] <- calculate_scoringdelta_estimator(Y, params, all_results)
  
  # c. Combine and return results
  return(unlist(all_results))
}


# --- 4. Set Simulation Parameters ---
# Configuration is defined directly in this script for clarity.
sim_params_cauchy <- list(
  simulation_name = "cauchy_delta_comparison",
  n_reps = 500, # Number of simulation repetitions
  n = 1000,     # Sample size
  burnin = 40,
  
  # Changepoint Specification
  changepoint_spec = list(
    tau = 500,
    delta = 1
  ),
  
  # Noise Specification (Cauchy)
  noise_dist = "cauchy",
  noise_params = list(scale = 1)
)


# --- 5. Execute the Simulation in Parallel ---
plan(
  multisession,
  workers = parallel::detectCores() - 1
)
registerDoFuture()
handlers(global = FALSE)

print("--- Starting Cauchy Delta Comparison Simulation ---")

with_progress({
  p <- progressor(steps = sim_params_cauchy$n_reps)
  
  results_df <- foreach(i = 1:sim_params_cauchy$n_reps, .combine = rbind) %dorng% {
    p()
    source("R/simulation_functions.R")
    source("R/noise_distributions.R")
    source("R/score_estimation.R")
    simulate_delta_comparison(params = sim_params_cauchy)
  }
})

# Shut down parallel workers
plan(sequential)
print("--- Simulation Complete ---")


# --- 6. Save and Analyze Results ---
output_data <- list(
  parameters = sim_params_cauchy,
  results = as.data.frame(results_df)
)

# Save the final output
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
filename <- paste0("results/cauchy_delta_sim_", timestamp, ".rds")
saveRDS(output_data, file = filename)

# Use your existing analysis tools to view the results
source("analysis/analysis_functions.R")
analysis_tools <- analyze_simulation_results(output_data)

print(paste("Results saved to:", filename))

# a. Calculate error metrics for the delta estimates
delta_errors <- analysis_tools$delta_error()
print("--- Delta Estimator Error Metrics (Ranked by MSE) ---")
print(delta_errors)

# b. Plot the distributions of the delta estimates
print("--- Plotting Delta Distributions ---")
analysis_tools$plot_delta()