# This script analyzes and visualises the output from single change-point
# simulation files.

# --- Helper Functions ---

#' Internal helper function to plot distributions for a given metric.
.plot_distributions <- function(data, column_suffix, title,
                                estimator_names = NULL, min_val = NULL,
                                max_val = NULL, plot_hist = TRUE, plot_density = TRUE) {
  
  all_cols <- grep(paste0("\\.", column_suffix, "$"), colnames(data), value = TRUE)
  
  if (!is.null(estimator_names)) {
    cols_to_plot <- paste0(estimator_names, ".", column_suffix)
    cols_to_plot <- intersect(cols_to_plot, all_cols)
  } else {
    cols_to_plot <- all_cols
  }
  
  if (length(cols_to_plot) == 0) {
    message("No columns found to plot after filtering for suffix '", column_suffix, "'.")
    return(invisible(NULL))
  }
  
  # --- Plot Setup ---
  legend_names <- sub(paste0("\\.", column_suffix, "$"), "", cols_to_plot)
  colors <- RColorBrewer::brewer.pal(max(3, length(cols_to_plot)), "Set1")
  
  x_range <- range(data[, cols_to_plot], na.rm = TRUE)
  if (!is.null(min_val)) x_range[1] <- min_val
  if (!is.null(max_val)) x_range[2] <- max_val
  
  max_y <- 0
  for (col_name in cols_to_plot) {
    if (plot_hist) {
      h <- hist(data[[col_name]], plot = FALSE, breaks = "FD")
      max_y <- max(max_y, max(h$density, na.rm = TRUE))
    }
    if (plot_density){
      d <- density(data[[col_name]], na.rm = TRUE)
      max_y <- max(max_y, max(d$y, na.rm = TRUE))
    }
   
  }
  max_y <- max_y * 1.1
  
  plot(NULL, xlim = x_range, ylim = c(0, max_y), main = title, xlab = "Value", ylab = "Density")
  
  for (i in 1:length(cols_to_plot)) {
    col_name <- cols_to_plot[i]
    if (plot_hist) {
      hist_color <- scales::alpha(colors[i], 0.3)
      hist(data[[col_name]], freq = FALSE, add = TRUE, col = hist_color, border = "grey", breaks = "FD")
    }
    if (plot_density){
      lines(density(data[[col_name]], na.rm = TRUE), col = colors[i], lwd = 2)
    }
  }
  
  legend("topright", legend = legend_names, col = colors, lwd = 2, bty = "n")
}

#' Internal helper function to calculate error metrics.
.calculate_error_metrics <- function(data, column_suffix, true_value, estimator_names = NULL) {
  if (missing(true_value) || !is.numeric(true_value) || length(true_value) != 1) {
    stop("You must provide a single numeric value for 'true_value'.")
  }
  
  all_cols <- grep(paste0("\\.", column_suffix, "$"), colnames(data), value = TRUE)
  
  if (!is.null(estimator_names)) {
    cols_to_calculate <- paste0(estimator_names, ".", column_suffix)
    cols_to_calculate <- intersect(cols_to_calculate, all_cols)
  } else {
    cols_to_calculate <- all_cols
  }
  
  if (length(cols_to_calculate) == 0) {
    message("No columns found to calculate errors for suffix '", column_suffix, "'.")
    return(NULL)
  }
  
  clean_estimator_names <- sub(paste0("\\.", column_suffix, "$"), "", cols_to_calculate)
  
  error_df <- data.frame(
    Estimator = clean_estimator_names,
    BIAS = sapply(data[, cols_to_calculate], function(est) mean(est - true_value, na.rm = TRUE)),
    MSE = sapply(data[, cols_to_calculate], function(est) mean((est - true_value)^2, na.rm = TRUE)),
    MAE = sapply(data[, cols_to_calculate], function(est) mean(abs(est - true_value), na.rm = TRUE)),
    MaxAE = sapply(data[, cols_to_calculate], function(est) max(abs(est - true_value), na.rm = TRUE)),
    row.names = NULL
  )

  error_df <- error_df[order(error_df$MSE), ]
  
  return(error_df)
}

# --- Main Analysis Function ---
analyze_simulation_results <- function(simulation_output) {
  results_df <- as.data.frame(simulation_output$results)
  data_params <- simulation_output$data_params
  
  plot_tau_distributions <- function(...) {
    .plot_distributions(results_df, "tau", "Distribution of Tau Estimators", ...)
  }
  plot_max_stat_distributions <- function(...) {
    .plot_distributions(results_df, "max", "Distribution of Max Test Statistics", ...)
  }
  plot_delta_distributions <- function(...) {
    .plot_distributions(results_df, "delta", "Distribution of Delta Estimators", ...)
  }
  
  calculate_tau_error <- function(...) {
    .calculate_error_metrics(results_df, "tau", data_params$changepoint_spec$tau, ...)
  }
  
  calculate_delta_error <- function(...){
    .calculate_error_metrics(results_df, "delta", data_params$changepoint_spec$delta, ...)
  }
  
  return(list(
    plot_tau = plot_tau_distributions,
    plot_max = plot_max_stat_distributions,
    plot_delta = plot_delta_distributions,
    tau_error = calculate_tau_error,
    delta_error = calculate_delta_error
  ))
}



results_filename <- "results/single_change_t_dist_refactored_2025-10-21_18-05-49.rds"
sim_output <- readRDS(results_filename)
analysis_tools <- analyze_simulation_results(sim_output)

# --- Use Examples ---


error_table <- analysis_tools$tau_error()
print("--- Error Metrics for Tau ---")
print(error_table)

 
# analysis_tools$plot_tau()
# analysis_tools$plot_tau(estimator_names = c("cusum", "scoring"), plot_hist = FALSE)