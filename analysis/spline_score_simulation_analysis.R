library(dplyr)
library(tibble) 
library(ggplot2)

get_summary_data <- function(results_data) {

  sample_sizes <- as.numeric(names(results_df))

  # Use sapply to iterate through the list. For each sample size, it will:
  # 1. unlist() the inner list of results to create a simple numeric vector.
  # 2. Compute the mean of that vector.
  mean_mses <- sapply(results_df, function(results_for_n) {
    mean(unlist(results_for_n))
  })

  # Do the same to calculate the standard deviation for each sample size.
  sd_mses <- sapply(results_df, function(results_for_n) {
    sd(unlist(results_for_n))
  })

  # Combine the calculated statistics into a new data frame for plotting
  summary_df <- data.frame(
    n = sample_sizes,
    MSE = mean_mses,
    SD = sd_mses
  )}

plot_mse_vs_n <- function(summary_data, use_log_scale = TRUE, ...) {
    # Check for required columns
    required_cols <- c("n", "MSE", "SD")
    if (!all(required_cols %in% colnames(summary_data))) {
      stop("Input data must contain 'n', 'MSE', and 'SD' columns.", call. = FALSE)
    }

    # Set up plotting arguments
    plot_args <- list(
      x = summary_data$n,
      y = summary_data$MSE,
      type = "b",
      pch = 16,
      col = "darkgreen",
      main = "Mean MSE vs. Sample Size",
      xlab = "Sample Size (n)",
      ylab = "Mean MSE",
      ylim = c(min(summary_data$MSE - summary_data$SD), max(summary_data$MSE + summary_data$SD)),
      ...
    )

    # Adjust for log scale if requested
    if (use_log_scale) {
      plot_args$log <- "xy"
      plot_args$main <- paste(plot_args$main, "(Log-Log Scale)")
      # ylim must be handled differently for log scale
      plot_args$ylim <- NULL
    }

    # Create the plot
    do.call(plot, plot_args)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

    # Add error bars for Standard Deviation
    arrows(
      x0 = summary_data$n,
      y0 = summary_data$MSE - summary_data$SD,
      x1 = summary_data$n,
      y1 = summary_data$MSE + summary_data$SD,
      length = 0.05,
      angle = 90,
      code = 3,
      col = scales::alpha("forestgreen", 0.6)
    )

    return(invisible(NULL))
  }

plot_mse_histogram <- function(results_data, n_to_plot) {
    # Convert the numeric n to a character to subset the list
    n_char <- as.character(n_to_plot)

    # Check if the requested n exists in the data
    if (!n_char %in% names(results_data)) {
      stop(paste("Data for n =", n_to_plot, "not found."), call. = FALSE)
    }

    # Extract the MSE values for the specified sample size
    mse_values <- unlist(results_data[[n_char]])

    # Create the histogram with a higher number of breaks for smaller bins
    hist(
      mse_values,
      main = paste("Histogram of MSEs for n =", n_char),
      xlab = "MSE",
      col = "lightblue",
      border = "darkblue",
      breaks = 100,
      xlim = c(0,1)# Increased number of breaks for smaller bins
    )

    return(invisible(NULL))
  }


#' @title Get summary statistics using a parameter list
#' @param results_data A data frame where rows are simulation runs
#'   (e.g., "result.1") and columns are "method.n" (e.g., "asm.100").
#' @param sim_params The configuration list used to run the simulation.
#'   This function uses `sim_params$estimators` and `sim_params$n_values`.
#' @return A tidy data frame (tibble) with summary statistics
#'   grouped by method and sample size.
get_long_data <- function(results_data, data_params, estimator_params) {
  
  all_results_list <- list()
  
  for (est_name in names(estimator_params)) {
    for (n_val in data_params$n_values) {
      
      col_name <- paste(est_name, n_val, sep = ".")
      
      # 4. Check if this column exists in the data frame
      if (col_name %in% names(results_data)) {
        
        # 5. Extract the vector of MSEs
        mse_vector <- results_data[[col_name]]
        
        # 6. Create a long-format data frame for this single combination
        temp_df <- data.frame(
          method = est_name,
          n = n_val,
          mse = mse_vector
        )
        
        # 7. Add this data frame to our list
        all_results_list[[col_name]] <- temp_df
      } else {
        # Optional: Warn if a column from params is missing
        warning(paste("Column '", col_name, "' not found in results_data."))
      }
    }
  }
  
  
  data_long <- do.call(rbind, all_results_list)
  
  return(data_long)
}

get_summary_df <- function(long_data){
  
  summary_df <- long_data %>%
    dplyr::group_by(method, n) %>%
    dplyr::summarise(
      mean_mse = mean(mse, na.rm = TRUE),
      sd_mse = sd(mse, na.rm = TRUE),
      min_mse = min(mse, na.rm = TRUE),
      lower_quartile_mse = quantile(mse, 0.25, na.rm = TRUE),
      median_mse = median(mse, na.rm = TRUE), # Median (Q2)
      upper_quartile_mse = quantile(mse, 0.75, na.rm = TRUE),
      max_mse = max(mse, na.rm = TRUE),
      .groups = 'drop' # Ungroup the resulting data frame
    ) %>%
    dplyr::arrange(method, n)
  return(summary_df)
}

generate_mse_box_plot <- function(long_data){
  
  mse_plot <- ggplot(long_data, aes(x = as.factor(n), y = mse, color = method)) +
    
    # Add boxplots
    # position_dodge(width = 0.8) places the boxplots side-by-side
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    
    # Add labels and titles
    labs(
      title = "MSE Distribution by Method and Sample Size",
      x = "Sample Size (n)",
      y = "Mean Squared Error (MSE)",
      color = "Estimator Method"
    ) +
    
    # Apply a clean theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  return(mse_plot)
}


results_address <- "results/score_est_t_multiple_se_2025-10-22_12-02-22.rds"
results_address <- "results/score_est_normal_multiple_se_2025-10-22_11-55-46.rds"
results_file <- readRDS(results_address)

data_params <- results_file$data_params
estimator_params <- results_file$estimator_params
results <- results_file$results

long_data <- get_long_data(results, data_params, estimator_params)
summary_data <- get_summary_df(long_data)
mse_plot <- generate_mse_box_plot(long_data)
print(mse_plot)
print(summary_data)
# results_df <- as.data.frame(results_data$results)
# summary_df <- get_summary_data(results_data)
# print(summary_df)
# 
# plot_mse_histogram(results_data$results, n_to_plot = 1000)


