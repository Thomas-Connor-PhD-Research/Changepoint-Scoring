library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

#' Calculate test statistic paths for multiple estimators
#'
#' @param data_params List. Parameters for data generation.
#' @param estimator_params List. Configuration for each estimator.
#' @param seed Integer. A master seed to set for reproducibility.
#'
#' @return A list containing a 'plot_data' data.frame and 'true_tau'.
calculate_statistic_paths <- function(data_params, estimator_params, seed = NULL){
  
  # 1. Set seed for reproducibility if provided
  if (!is.null(seed)) {
    # Use the correct method to set the stream for this *single* run
    .Random.seed <- get_dorng_stream(seed, 1) # Assumes iteration 1
    # Or, if you're just running this once, set.seed(seed) is fine.
    # Let's use the stream method as in your example.
    # If the seed is just a simple seed, not a stream:
    # set.seed(seed) 
  }
  
  # 2. Generate data
  data <- generate_timeseries(data_params$n, 
                              data_params$changepoint_spec$tau, 
                              data_params$changepoint_spec$delta, 
                              data_params$noise_dist, 
                              data_params$noise_params)
  Y <- data$Y
  noise <- data$noise
  
  # 3. Initialize a list to store results
  all_paths <- list()
  
  # 4. Loop through the estimator configurations
  for (estimator in names(estimator_params)) {
    
    config <- estimator_params[[estimator]]
    results <- NULL
    
    tryCatch({
      results <- switch(estimator,
                        "cusum" = {
                          calculate_cusum_estimator(Y, return_path = TRUE)
                        },
                        "oracle" = {
                          calculate_oracle_estimator(Y, 
                                                     data_params$noise_dist, 
                                                     data_params$noise_params,
                                                     return_path = TRUE)
                        },
                        "scoring_known_noise" = {
                          calculate_scoring_known_noise_estimator(Y, noise, 
                                                                  config$scoring_method, 
                                                                  return_path = TRUE)
                        },
                        "scoring" = {
                          calculate_scoring_estimator(Y,
                                                      scoring_method = config$scoring_method,
                                                      iterations = ifelse(is.null(config$iterations), 1, config$iterations),
                                                      rho_initial = config$rho_initial, 
                                                      return_path = TRUE
                          )
                        },
                        {
                          stop(paste("Unknown estimator type:", estimator))
                        }
      )
    }, error = function(e) {
      warning(paste("Failed to run estimator", estimator, ":", e$message))
    })
    
    if (!is.null(results) && !is.null(results$path)) {
      all_paths[[estimator]] <- data.frame(
        k = 1:length(results$path),
        test_statistic = results$path,
        estimator = estimator,
        tau_hat = results$tau
      )
    }
  }
  
  # 5. Combine all data frames into one
  if (length(all_paths) == 0) {
    stop("No estimator paths were successfully generated.")
  }
  
  return(list(
    timeseries = data,
    plot_data = bind_rows(all_paths),
    true_tau = data_params$changepoint_spec$tau
  ))
}

#' Plot absolute test statistic paths
#'
#' @param path_results List. The output from 'calculate_statistic_paths'.
#' @param x_min Numeric. Optional minimum x-axis limit.
#' @param x_max Numeric. Optional maximum x-axis limit.
#'
#' @return A ggplot object.
plot_absolute_paths <- function(path_results, x_min = NULL, x_max = NULL) {
  
  plot_data <- path_results$plot_data
  true_tau <- path_results$true_tau
  
  p <- ggplot(plot_data, aes(x = k, y = test_statistic, color = estimator)) +
    geom_line(linewidth = 0.8, aes(linetype = estimator)) + 
    
    # Add a vertical line for the true CP
    geom_vline(xintercept = true_tau, linetype = "dashed", color = "black", linewidth = 1) +
    
    # Add vertical lines for each estimator's guess (will be colored)
    geom_vline(aes(xintercept = tau_hat, color = estimator), linetype = "dotted", linewidth = 0.8) +
    
    labs(
      title = "Test Statistic Paths for Change-Point Estimators",
      subtitle = paste("True Change-Point (tau) =", true_tau, "(dashed black line)"),
      x = "Potential Change-Point (k)",
      y = "Test Statistic Value"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(x_min, x_max))
  
  return(p)
}

#' Plot the difference between two test statistic paths
#'
#' @param path_results List. The output from 'calculate_statistic_paths'.
#'   This *must* contain exactly two estimators.
#' @param x_min Numeric. Optional minimum x-axis limit.
#' @param x_max Numeric. Optional maximum x-axis limit.
#'
#' @return A ggplot object.
plot_difference_path <- function(path_results, x_min = NULL, x_max = NULL) {
  
  # 1. Check for exactly two estimators
  estimators <- unique(path_results$plot_data$estimator)
  if (length(estimators) != 2) {
    stop(paste(
      "plot_difference_path requires *exactly* 2 estimators, but found", 
      length(estimators), ": (", paste(estimators, collapse = ", "), ")"
    ))
  }
  
  # 2. Get data and true tau
  plot_data <- path_results$plot_data
  true_tau <- path_results$true_tau
  
  # 3. Calculate the difference
  # Use tidyr::pivot_wider to get estimators as columns
  wide_data <- plot_data %>%
    select(k, estimator, test_statistic) %>%
    tidyr::pivot_wider(names_from = estimator, values_from = test_statistic)
  
  # Create the difference column. We use .data[[]] to access columns by string name
  diff_data <- wide_data %>%
    mutate(difference = .data[[estimators[1]]] - .data[[estimators[2]]])
  
  # Get the tau_hat values for the dotted lines
  tau_hats <- plot_data %>% 
    select(estimator, tau_hat) %>% 
    distinct()
  
  # 4. Create plot
  plot_title <- paste("Difference Plot:", estimators[1], "-", estimators[2])
  
  p <- ggplot(diff_data, aes(x = k, y = difference)) +
    geom_line(linewidth = 0.8, color = "black") +
    
    # Add a vertical line for the true CP
    geom_vline(xintercept = true_tau, linetype = "dashed", color = "black", linewidth = 1) +
    
    # Add vertical lines for *both* estimators' guesses
    geom_vline(data = tau_hats, 
               aes(xintercept = tau_hat, color = estimator), 
               linetype = "dotted", linewidth = 0.8) +
    
    # Add a horizontal line at 0
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    
    labs(
      title = plot_title,
      subtitle = paste("True Change-Point (tau) =", true_tau, "(dashed black line)"),
      x = "Potential Change-Point (k)",
      y = "Difference in Test Statistic"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = c(x_min, x_max))
  
  return(p)
}

# --- 4. EXAMPLE USAGE ---

seed = 123
stream_to_run = 101 # Specify the iteration you want to debug

# Data Generation Parameters
data_params <- list(
  n = 1000,
  changepoint_spec = list(
    tau = 500,
    delta = 0.6
  ),
  noise_dist = "t",
  noise_params = list(scale = 1, df = 3)
)

# Estimator Parameters
estimator_params_for_diff <- list(
  "oracle" = list(), 
  "scoring" = list(
    scoring_method = "spline_df_min", 
    iterations = 1
  )
)




# --- 1. Calculate the paths ONCE ---
# We pass the master seed and the specific stream number
.Random.seed <- get_dorng_stream(seed, stream_to_run)
path_results <- calculate_statistic_paths(data_params, estimator_params_for_diff)


# --- 2. Create the plots ---

# Plot 1: Absolute Paths
plot_abs <- plot_absolute_paths(path_results)
print(plot_abs)

# Plot 2: Difference Path (zoomed in)
plot_diff <- plot_difference_path(path_results, x_min = NULL, x_max = NULL)
print(plot_diff)

Y <- path_results$timeseries$Y
plot(1:1000, Y)