# --- 0. LIBRARIES AND PLACEHOLDER FUNCTIONS ---
# Make sure to install these packages
# install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)



# --- 1. CORE HELPER FUNCTIONS ---

#' Core Logic: Calculate all intermediate vectors
detect_cp_and_get_all_internals <- function(Y, rho) {
  n <- length(Y)
  k <- 1:(n - 1)
  
  centered_rho_y <- cumsum(rho(Y) - mean(rho(Y)))
  delta_scaling <- n / (k * (n - k))
  possible_deltas <- delta_scaling * centered_rho_y[k]
  
  constant_variance_scaling <- sqrt((k * (n - k)) / n)
  possible_deltas[is.nan(possible_deltas)] <- 0
  
  cusum_statistic_un_squared <- constant_variance_scaling * possible_deltas
  test_statistics_squared <- cusum_statistic_un_squared^2
  
  tau_hat <- which.max(test_statistics_squared)
  if (length(tau_hat) > 1) tau_hat <- tau_hat[1] 
  
  test_stat_max <- test_statistics_squared[tau_hat]
  delta_hat <- possible_deltas[tau_hat]
  
  return(list(
    estimates = list(
      tau = tau_hat,
      delta = delta_hat,
      max = test_stat_max
    ),
    possible_deltas = possible_deltas,
    cusum_statistic_un_squared = cusum_statistic_un_squared,
    test_statistics_squared = test_statistics_squared
  ))
}


#' Core Wrapper: Create the interactive object
#'
#' This (UPDATED) function no longer returns plotters.
#' It now includes an 'object_name' for use in legends.
create_investigation_object <- function(object_name, # <-- NEW
                                        Y,
                                        all_internals,
                                        rho_fn,
                                        true_score_fn = NULL,
                                        oracle_fi = NULL,
                                        est_score_fn = NULL,
                                        est_fi = NULL,
                                        noise_estimate = NULL) {
  force(object_name)
  force(Y)
  force(all_internals)
  force(rho_fn)
  force(true_score_fn)
  force(oracle_fi)
  force(est_score_fn)
  force(est_fi)
  force(noise_estimate)
  
  list(
    # --- Accessor Functions ---
    get_name = function() object_name, # <-- NEW
    get_Y = function() Y,
    get_rho = function() rho_fn,
    get_true_score_fn = function() true_score_fn,
    get_estimated_score_fn = function() est_score_fn,
    get_oracle_fi = function() oracle_fi,
    get_estimated_fi = function() est_fi,
    get_noise_estimate = function() noise_estimate,
    
    # --- Data Vectors ---
    get_possible_deltas = function() all_internals$possible_deltas,
    get_cusum_statistic = function() all_internals$cusum_statistic_un_squared,
    get_test_statistics = function() all_internals$test_statistics_squared,
    
    # --- Final Estimates ---
    get_estimates = function() all_internals$estimates,
    
    # --- Comparison Functions ---
    compare_fi = function() {
      message(paste("--- FI Comparison for:", object_name, "---"))
      if (!is.null(oracle_fi)) {
        message(paste("Oracle FI (True):   ", format(oracle_fi, digits = 5)))
      }
      if (!is.null(est_fi)) {
        message(paste("Estimated FI:       ", format(est_fi, digits = 5)))
      }
      if (!is.null(oracle_fi) && !is.null(est_fi)) {
        message(paste("Difference (Est - True):", format(est_fi - oracle_fi, digits = 5)))
      }
    }
    
    # --- Plotting functions have been REMOVED ---
  )
}


# --- 2. MAIN DEBUG FUNCTIONS (Updated) ---

#' Create investigation object for CUSUM
debug_cusum_estimator <- function(Y, noise_dist, noise_params) {
  
  rho <- function(x) -x
  true_score_fn <- create_score_function(noise_dist, noise_params)
  oracle_fi <- get_stein_fisher_info(noise_dist, noise_params)
  all_internals <- detect_cp_and_get_all_internals(Y, rho)
  
  return(create_investigation_object(
    object_name = "cusum", # <-- Pass name
    Y = Y,
    all_internals = all_internals,
    rho_fn = rho,
    true_score_fn = true_score_fn,
    oracle_fi = oracle_fi
  ))
}


#' Create investigation object for Oracle
debug_oracle_estimator <- function(Y, noise_dist, noise_params) {
  
  oracle_score_fn <- create_score_function(noise_dist, noise_params)
  oracle_fi <- get_stein_fisher_info(noise_dist, noise_params)
  rho <- function(x) oracle_score_fn(x) / oracle_fi
  
  all_internals <- detect_cp_and_get_all_internals(Y, rho)
  
  return(create_investigation_object(
    object_name = "oracle", # <-- Pass name
    Y = Y,
    all_internals = all_internals,
    rho_fn = rho,
    true_score_fn = oracle_score_fn,
    oracle_fi = oracle_fi
  ))
}


#' Create investigation object for Scoring with Known Noise
debug_scoring_known_noise_estimator <- function(Y, noise, scoring_method,
                                                noise_dist, noise_params) {
  
  score_fn_estimate <- score_estimation(noise, scoring_method)[[scoring_method]]
  fi_estimate <- var(score_fn_estimate(noise))
  rho <- function(x) score_fn_estimate(x) / fi_estimate
  
  true_score_fn <- create_score_function(noise_dist, noise_params)
  oracle_fi <- get_stein_fisher_info(noise_dist, noise_params)
  
  all_internals <- detect_cp_and_get_all_internals(Y, rho)
  
  return(create_investigation_object(
    object_name = paste0("scoring_known_noise (", scoring_method, ")"), # <-- Pass name
    Y = Y,
    all_internals = all_internals,
    rho_fn = rho,
    true_score_fn = true_score_fn,
    oracle_fi = oracle_fi,
    est_score_fn = score_fn_estimate,
    est_fi = fi_estimate,
    noise_estimate = noise
  ))
}

debug_scoring_estimator <- function(Y,
                                    scoring_method,
                                    noise_dist, noise_params,
                                    iterations = 1,
                                    rho_initial = NULL) {
  
  # --- NEW: Function factory to fix lazy evaluation ---
  # This creates a *new* function, "locking in" the values of
  # est_fn and est_fi so they don't change later.
  make_rho_fn <- function(est_fn, est_fi) {
    force(est_fn)
    force(est_fi)
    function(x) {
      est_fn(x) / est_fi
    }
  }
  # --- END NEW ---
  
  n <- length(Y)
  investigation_steps <- list()
  true_score_fn <- create_score_function(noise_dist, noise_params)
  oracle_fi <- get_stein_fisher_info(noise_dist, noise_params)
  
  current_noise_estimate <- Y
  if (!is.null(rho_initial)) {
    current_rho <- rho_initial
    # We still estimate from Y for the iter_0 object
    current_est_score_fn <- score_estimation(Y, scoring_method)[[scoring_method]]
    current_est_fi <- var(current_est_score_fn(Y))
  } else {
    current_est_score_fn <- score_estimation(Y, scoring_method)[[scoring_method]]
    current_est_fi <- var(current_est_score_fn(Y))
    # --- FIX: Use the factory ---
    current_rho <- make_rho_fn(current_est_score_fn, current_est_fi)
  }
  
  internals <- detect_cp_and_get_all_internals(Y, current_rho)
  
  investigation_steps[["iter_0"]] <- create_investigation_object(
    object_name = "iter_0",
    Y = Y,
    all_internals = internals,
    rho_fn = current_rho,
    true_score_fn = true_score_fn,
    oracle_fi = oracle_fi,
    est_score_fn = current_est_score_fn,
    est_fi = current_est_fi,
    noise_estimate = current_noise_estimate
  )
  
  if (iterations > 0) {
    for (i in 1:iterations) {
      prev_estimates <- investigation_steps[[paste0("iter_", i - 1)]]$get_estimates()
      tau_hat <- prev_estimates$tau
      delta_hat <- prev_estimates$delta
      
      current_noise_estimate <- numeric(n)
      current_noise_estimate[1:tau_hat] <- Y[1:tau_hat]
      current_noise_estimate[(tau_hat + 1):n] <- Y[(tau_hat + 1):n] - delta_hat
      
      current_est_score_fn <- score_estimation(current_noise_estimate, scoring_method)[[scoring_method]]
      current_est_fi <- var(current_est_score_fn(current_noise_estimate))
      
      # --- FIX: Use the factory ---
      current_rho <- make_rho_fn(current_est_score_fn, current_est_fi)
      
      internals <- detect_cp_and_get_all_internals(Y, current_rho)
      
      investigation_steps[[paste0("iter_", i)]] <- create_investigation_object(
        object_name = paste0("iter_", i),
        Y = Y,
        all_internals = internals,
        rho_fn = current_rho,
        true_score_fn = true_score_fn,
        oracle_fi = oracle_fi,
        est_score_fn = current_est_score_fn,
        est_fi = current_est_fi,
        noise_estimate = current_noise_estimate
      )
    }
  }
  
  return(investigation_steps)
}


# --- 3. NEW EXTERNAL PLOTTING FUNCTIONS ---

#' Plot statistic vectors from multiple investigation objects
#'
#' @param ... One or more investigation objects.
#' @param statistic_to_plot A string: "test_statistics" (default),
#'   "cusum_statistic", or "possible_deltas".
#' @param title An optional plot title.
plot_statistic_vectors <- function(..., statistic_to_plot = "test_statistics", title = NULL) {
  
  # 1. Capture all investigation objects
  objects <- list(...)
  if (length(objects) == 0) {
    stop("You must provide at least one investigation object.")
  }
  
  # 2. Map the statistic name to the correct getter function
  getter_fn_name <- switch(statistic_to_plot,
                           "test_statistics" = "get_test_statistics",
                           "cusum_statistic" = "get_cusum_statistic",
                           "possible_deltas" = "get_possible_deltas",
                           stop("Invalid 'statistic_to_plot'. Must be one of: 'test_statistics', 'cusum_statistic', 'possible_deltas'")
  )
  
  # 3. Build the data frame
  all_paths <- list()
  for (obj in objects) {
    vec <- obj[[getter_fn_name]]()
    all_paths[[obj$get_name()]] <- data.frame(
      k = 1:length(vec),
      value = vec,
      name = obj$get_name(),
      tau_hat = obj$get_estimates()$tau
    )
  }
  plot_data <- bind_rows(all_paths)
  
  # 4. Set default title
  if (is.null(title)) {
    title <- paste("Comparison of:", statistic_to_plot)
  }
  
  # 5. Create plot
  ggplot(plot_data, aes(x = k, y = value, color = name)) +
    geom_line(linewidth = 0.8) +
    # Add vertical lines for each estimator's guess
    geom_vline(aes(xintercept = tau_hat, color = name), linetype = "dotted", linewidth = 0.8) +
    labs(
      title = title,
      x = "Potential Change-Point (k)",
      y = "Statistic Value",
      color = "Estimator"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


#' Plot score functions from multiple investigation objects
#'
#' (UPDATED with fix for NULL est_fn)
#'
#' @param ... One or more investigation objects.
#' @param plot_rho Boolean. Plot the obj$get_rho() function? (Default TRUE)
#' @param plot_est_score Boolean. Plot the obj$get_estimated_score_fn()? (Default TRUE)
#' @param plot_true_score Boolean. Plot the obj$get_true_score_fn()? (Default TRUE)
#' @param x_grid Numeric vector. An optional grid to evaluate functions on.
plot_score_functions <- function(..., plot_rho = TRUE, plot_est_score = TRUE, plot_true_score = TRUE, x_grid = NULL) {
  
  # 1. Capture all investigation objects
  objects <- list(...)
  if (length(objects) == 0) {
    stop("You must provide at least one investigation object.")
  }
  
  # 2. Determine plot grid if not provided
  if (is.null(x_grid)) {
    Y <- objects[[1]]$get_Y()
    x_grid <- seq(min(Y, na.rm = TRUE), max(Y, na.rm = TRUE), length.out = 200)
  }
  
  # 3. Build the data frame
  all_scores <- list()
  
  # Add the True Score (only once)
  if (plot_true_score) {
    true_fn <- objects[[1]]$get_true_score_fn()
    if (!is.null(true_fn)) {
      true_values <- true_fn(x_grid)
      all_scores[["true_score"]] <- data.frame(
        x = x_grid,
        value = true_values,
        name = "True Score"
      )
    }
  }
  
  # Loop through objects for rho and est_score
  for (obj in objects) {
    obj_name <- obj$get_name()
    
    # Add rho function
    if (plot_rho) {
      rho_fn <- obj$get_rho()
      rho_values <- rho_fn(x_grid)
      all_scores[[paste0(obj_name, "_rho")]] <- data.frame(
        x = x_grid,
        value = rho_values,
        name = paste0(obj_name, " (rho)")
      )
    }
    
    # Add estimated score function
    if (plot_est_score) {
      est_fn <- obj$get_estimated_score_fn()
      # --- FIX: Check for NULL *before* calling the function ---
      if (!is.null(est_fn)) {
        est_values <- est_fn(x_grid) 
        all_scores[[paste0(obj_name, "_est_score")]] <- data.frame(
          x = x_grid,
          value = est_values,
          name = paste0(obj_name, " (Estimated Score)")
        )
      }
    }
  }
  
  plot_data <- bind_rows(all_scores)
  
  # 4. Create plot
  ggplot(plot_data, aes(x = x, y = value, color = name)) +
    geom_line(linewidth = 0.8) +
    labs(
      title = "Score Function Comparison",
      x = "x",
      y = "Score Value",
      color = "Function"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


# --- 4. EXAMPLE USAGE ---

# --- Setup ---
seed <- 123
stream_to_run <- 65426 # Specify the iteration you want to debug

data_params <- list(
  n = 1000,
  changepoint_spec = list(
    tau = 500,
    delta = 0.3
  ),
  noise_dist = "t",
  noise_params = list(scale = 0.5, df = 3)
)

# Set the seed for this specific run
.Random.seed <- get_dorng_stream(seed, stream_to_run)

# Generate the data *once*
data <- generate_timeseries(
  data_params$n,
  data_params$changepoint_spec$tau,
  data_params$changepoint_spec$delta,
  data_params$noise_dist,
  data_params$noise_params
)
Y <- data$Y
noise <- data$noise

# --- Example 1: Comparing CUSUM and Oracle ---
message("--- Running CUSUM and Oracle debuggers ---")
cusum_obj <- debug_cusum_estimator(Y, data_params$noise_dist, data_params$noise_params)
oracle_obj <- debug_oracle_estimator(Y, data_params$noise_dist, data_params$noise_params)

# Plot their FINAL test statistics on the same graph
# To see the plot, run this line interactively:
# print(plot_statistic_vectors(cusum_obj, oracle_obj, statistic_to_plot = "test_statistics"))

# Plot their CUSUM (un-squared) statistics on the same graph
# To see the plot, run this line interactively:
# print(plot_statistic_vectors(cusum_obj, oracle_obj, statistic_to_plot = "cusum_statistic"))

# Plot their rho functions against the true score
# To see the plot, run this line interactively:
# print(plot_score_functions(cusum_obj, oracle_obj, plot_est_score = FALSE))


# --- Example 2: Comparing Iterations of the Scoring Estimator ---
message("\n--- Running debug_scoring_estimator(iterations = 2) ---")
debug_obj_list <- debug_scoring_estimator(
  Y = Y,
  scoring_method = "spline_df_min",
  noise_dist = data_params$noise_dist,
  noise_params = data_params$noise_params,
  iterations = 1
)

# Compare the test statistics of iter_0, iter_1, and iter_2
# To see the plot, run this line interactively:
# print(plot_statistic_vectors(
#   debug_obj_list$iter_0,
#   debug_obj_list$iter_1,
#   debug_obj_list$iter_2,
#   statistic_to_plot = "test_statistics"
# ))

# Compare the estimated score functions from all 3 steps
# To see the plot, run this line interactively:
print(plot_score_functions(
  debug_obj_list$iter_0,
  debug_obj_list$iter_1,
  plot_rho = FALSE,
  plot_est_score = TRUE# Don't plot rho, just the estimated scores
))

# Compare the FI of the first and last step
message("\n--- Comparing FI for first and last steps ---")
debug_obj_list$iter_0$compare_fi()

# Get final estimates
# message("\n--- Final estimates from iter_2 ---")
# print(debug_obj_list$iter_2$get_estimates())