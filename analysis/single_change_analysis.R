# This script provides functions to analyse and visualise the output
# from sing_change2-2.R: plot histograms / kernel density estimates of tau (and
# corresponding lambda) and max. test statistic
#

# install.packages(c("RColorBrewer", "scales"))

#' A helper function to plot distributions for a given metric.
#'
#' This internal function generates a plot with histograms and overlaid
#' kernel density estimates for all columns in a dataframe that match a
#' specific suffix (e.g., ".tau" or ".max").
#'
#' @param data The results dataframe from the simulation.
#' @param column_suffix The suffix of the columns to plot (e.g., "tau").
#' @param title The main title for the plot.
#' @return Invisibly returns NULL. The function's purpose is to draw a plot.
plot_distributions <- function(data, column_suffix, title, estimator_names = NULL) {
  # Find cols that match a suffix
  cols_to_plot <- grep(paste0("\\.", column_suffix, "$"), colnames(data), value = TRUE)
  
  if (length(cols_to_plot) == 0) {
    message("No columns found with suffix '.", column_suffix, "' to plot.")
    return(invisible(NULL))
  }
  # Extract clean estimator names if not provided (check functionality)
  if (is.null(estimator_names)){
  estimator_names <- sub(paste0("\\.", column_suffix, "$"), "", cols_to_plot)
  }
  
  # 2. Set up colors for plotting and calculate plot boundaries
  colors <- RColorBrewer::brewer.pal(max(3, length(cols_to_plot)), "Set1")
  x_range <- range(data[, cols_to_plot], na.rm = TRUE)
  
  # Calculate the maximum y-value needed to fit all histograms and density curves
  max_y <- 0
  for (col_name in cols_to_plot) {
    # Use "FD" algorithm for a reasonable number of histogram bins
    h <- hist(data[[col_name]], plot = FALSE, breaks = "FD")
    d <- density(data[[col_name]], na.rm = TRUE)
    max_y <- max(max_y, max(h$density), max(d$y), na.rm = TRUE)
  }
  max_y <- max_y * 1.1 # Add a 10% buffer to the top of the plot
  
  # 3. Create an empty plot with the calculated dimensions
  plot(NULL, xlim = x_range, ylim = c(0, max_y),
       main = title,
       xlab = "Value",
       ylab = "Density")
  
  # 4. Loop through each estimator to draw its histogram and density curve
  for (i in 1:length(cols_to_plot)) {
    col_name <- cols_to_plot[i]
    # Use semi-transparent colors for histogram bars to see overlaps
    hist_color <- scales::alpha(colors[i], 0.3)
    hist(data[[col_name]], freq = FALSE, add = TRUE, col = hist_color, border = "grey", breaks = "FD")
    lines(density(data[[col_name]], na.rm = TRUE), col = colors[i], lwd = 2)
  }
  
  # 5. Add a legend to identify the estimators
  legend("topright",
         legend = estimator_names,
         col = colors,
         lwd = 2,
         bty = "n") # bty="n" removes the box around the legend
}


#' Main analysis function to generate plotting tools.
#'
#' This is the primary function you will call. It takes a simulation results
#' dataframe and returns a list of functions that you can use to generate plots.
#'
#' @param results_df The dataframe containing your simulation results.
#' @return A list containing two functions: `plot_tau()` and `plot_max()`.
analyze_simulation_results <- function(results_df) {
  results_df <- as.data.frame(results_df)
  
  # Function for plotting the dist of tau estimates
  plot_tau_distributions <- function() {
    plot_distributions(results_df, "tau", "Distribution of Tau Estimators")
  }
  
  # Function for plotting max statistic dist
  plot_max_stat_distributions <- function() {
    plot_distributions(results_df, "max", "Distribution of Max Test Statistics")
  }
  
  # Function for plotting the dist of the lambda estimates
  plot_lambda_distributions <- function(n) {
    lambda_cols_exist <- any(grepl("\\.lambda$", colnames(results_df)))
    
    # If lambda columns don't exist, create and append them
    if (!lambda_cols_exist) {
      if (missing(n) || !is.numeric(n) || length(n) != 1 || n <= 0) {
        stop("First time calling plot_lambda(): You must provide a single positive number for the sample size 'n'.")
      }
      
      # Find all columns ending in ".tau"
      tau_cols <- grep("\\.tau$", colnames(results_df), value = TRUE)
      if (length(tau_cols) == 0) {
        message("No '.tau' columns found to convert to lambda.")
        return(invisible(NULL))
      }
      
      # Create a new dataframe containing only the new lambda columns
      lambda_df_new_cols <- results_df[, tau_cols, drop = FALSE] / n
      colnames(lambda_df_new_cols) <- sub("\\.tau$", ".lambda", colnames(lambda_df_new_cols))
      
      # Append lambda cols for future use
      results_df <<- cbind(results_df, lambda_df_new_cols)
      message("Created and appended .lambda columns to the results dataframe.")
    }
    
    plot_distributions(results_df, "lambda", "Distribution of Lambda Estimators (tau / n)")
  }
  
  
  calculate_error_metrics <- function(true_tau) {
    if (missing(true_tau) || !is.numeric(true_tau) || length(true_tau) != 1) {
      stop("You must provide a single numeric value for 'true_tau'.")
    }
    
    tau_cols <- grep("\\.tau$", colnames(results_df), value = TRUE)
    if (length(tau_cols) == 0) {
      message("No '.tau' columns found to calculate errors.")
      return(invisible(NULL))
    }
    estimator_names <- sub("\\.tau$", "", tau_cols)
    
    # Calculate errors versus the true tau 
    errors_vs_true <- data.frame(
      Estimator = estimator_names,
      MSE = sapply(results_df[, tau_cols], function(est) mean((est - true_tau)^2, na.rm = TRUE)),
      MAE = sapply(results_df[, tau_cols], function(est) mean(abs(est - true_tau), na.rm = TRUE)),
      MaxAE = sapply(results_df[, tau_cols], function(est) max(abs(est - true_tau), na.rm = TRUE)),
      row.names = NULL
    )
    
    # Calculate errors versus the CUSUM estimates (if available)
    errors_vs_cusum <- NULL
    if ("cusum.tau" %in% tau_cols) {
      cusum_estimates <- results_df[["cusum.tau"]]
      other_cols <- tau_cols[tau_cols != "cusum.tau"]
      
      if (length(other_cols) > 0) {
        errors_vs_cusum <- data.frame(
          Estimator = sub("\\.tau$", "", other_cols),
          MSE_vs_CUSUM = sapply(results_df[, other_cols], function(est) mean((est - cusum_estimates)^2, na.rm = TRUE)),
          MAE_vs_CUSUM = sapply(results_df[, other_cols], function(est) mean(abs(est - cusum_estimates), na.rm = TRUE)),
          MaxAE_vs_CUSUM = sapply(results_df[, other_cols], function(est) max(abs(est - cusum_estimates), na.rm = TRUE)),
          row.names = NULL
        )
      }
    }
    
    # Return a list containing the two results dataframes
    return(list(
      vs_True_Tau = errors_vs_true,
      vs_CUSUM = errors_vs_cusum
    ))
  }
  # Return the callable functions in a list
  return(list(
    plot_tau = plot_tau_distributions,
    plot_max = plot_max_stat_distributions,
    plot_lambda = plot_lambda_distributions,
    calculate_errors = calculate_error_metrics
  ))
}


sim_results <- readRDS("t_1000_0-5_3_0-5_0-3_sim_results.rds")
analysis_tools <- analyze_simulation_results(sim_results)

error_results <- analysis_tools$calculate_errors(true_tau = 500)
print(error_results$vs_True_Tau)

analysis_tools$plot_tau()
# analysis_tools$plot_lambda(n = 1000)
