# ============================================================
# Convergence & Error Analysis for CP Rate Estimation
# ============================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)

# ------------------------------------------------------------
# Helper: Riemann zeta(3) constant (approximate)
#   Used for the theoretical Bayes RMSE:
#     4 * sqrt(zeta(3)) / r^2
# ------------------------------------------------------------
zeta_3_value <- function() {
  # Numeric value of Riemann zeta(3)
  1.2020569031595942854
}

# ------------------------------------------------------------
# Helper: prepare long-format lambda results
#
# - Takes a wide results data.frame with columns like:
#     estimator.metric_col.n   (e.g. "method.tau.1000")
# - Produces a long data.frame with:
#     estimator, metric_col, n, tau_hat, lambda_hat, error
# - Optionally trims lambda_hat outside [trim_lambda, 1 - trim_lambda]
# ------------------------------------------------------------
prepare_lambda_results <- function(results_df,
                                   true_lambda,
                                   trim_lambda = 0) {
  # Long format with columns: estimator, metric_col, n, tau_hat
  results_long <- results_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "col_name",
      values_to = "tau_hat"
    ) %>%
    # Expected col_name pattern: estimator.metric_col.n
    separate(
      col = col_name,
      into = c("estimator", "metric_col", "n"),
      sep = "\\.",
      convert = FALSE,
      extra = "merge",
      fill = "right"
    )
  
  lambda_results <- results_long %>%
    mutate(
      tau_hat    = as.numeric(tau_hat),
      n          = as.numeric(n),
      lambda_hat = tau_hat / n,
      error      = lambda_hat - true_lambda
    )
  
  # Basic input checks for trimming
  if (!is.numeric(trim_lambda) || trim_lambda < 0 || trim_lambda >= 0.5) {
    stop("trim_lambda must be numeric in [0, 0.5).")
  }
  
  # Optional trimming by lambda_hat
  if (trim_lambda > 0) {
    cat(sprintf(
      "\n--- Trimming: keeping lambda_hat in [%.4f, %.4f] ---\n",
      trim_lambda, 1 - trim_lambda
    ))
    
    lambda_results <- lambda_results %>%
      filter(
        !is.na(lambda_hat),
        lambda_hat >= trim_lambda,
        lambda_hat <= (1 - trim_lambda)
      )
  }
  
  lambda_results
}

# ============================================================
# 1) Plot Estimator Convergence Metric
#    + Theoretical MLE & Bayes RMSE curves (for RMSE only)
# ============================================================
# Arguments:
#   - results_df:   data frame of simulation results (wide format)
#   - true_lambda:  true lambda used in the data-generating process
#   - metric:       "rmse" or "mae"
#   - scale_target: "lambda" or "tau" (which parameter is the error on)
#   - log_scale:    if TRUE, use log-log axes
#   - trim_lambda:  trimming in prepare_lambda_results()
#   - estimators:   optional subset of estimators to keep
#   - alpha:        multiplies metric by n^alpha
#   - delta:        changepoint magnitude (can be:
#                     * scalar numeric
#                     * function(n) returning numeric
#                     * numeric vector aligned with unique n values)
#   - fisher_info:  I(f), Fisher information (usually scalar numeric)
# Notes:
#   - Theoretical curves are added only when metric == "rmse"
#   - Theoretical tau-RMSE:
#       MLE:   sqrt(26) / r^2
#       Bayes: 4 * sqrt(zeta(3)) / r^2
#       where r = delta * |fisher_info|
#   - For lambda-scale, theory is divided by n and then scaled by n^alpha
# ============================================================
plot_convergence_metric <- function(results_df,
                                    true_lambda,
                                    metric = "rmse",
                                    scale_target = "lambda",
                                    log_scale = TRUE,
                                    trim_lambda = 0,
                                    estimators = NULL,
                                    alpha = 0,
                                    delta,
                                    fisher_info) {
  
  # -----------------------------------
  # Prepare data
  # -----------------------------------
  lambda_results <- prepare_lambda_results(
    results_df   = results_df,
    true_lambda  = true_lambda,
    trim_lambda  = trim_lambda
  )
  
  # Optional estimator filter (applied early to save work)
  if (!is.null(estimators)) {
    lambda_results <- lambda_results %>%
      filter(estimator %in% estimators)
  }
  
  # -----------------------------------
  # Argument checks
  # -----------------------------------
  metric <- tolower(metric)
  if (!metric %in% c("rmse", "mae")) {
    stop("Invalid metric. Choose 'rmse' or 'mae'.")
  }
  
  scale_target <- tolower(scale_target)
  if (!scale_target %in% c("lambda", "tau")) {
    stop("Invalid scale_target. Choose 'lambda' or 'tau'.")
  }
  
  # -----------------------------------
  # Choose error column: lambda or tau
  # -----------------------------------
  if (scale_target == "lambda") {
    error_col    <- "error"
    y_label_base <- "Lambda"
  } else { # scale_target == "tau"
    lambda_results <- lambda_results %>%
      mutate(tau_error = (lambda_hat - true_lambda) * n)
    error_col    <- "tau_error"
    y_label_base <- "Tau"
  }
  
  # -----------------------------------
  # Aggregate error into RMSE or MAE by (estimator, n)
  # -----------------------------------
  if (metric == "rmse") {
    error_data <- lambda_results %>%
      filter(!is.na(.data[[error_col]])) %>%
      group_by(estimator, n) %>%
      summarise(
        Metric_Value = sqrt(mean((.data[[error_col]])^2, na.rm = TRUE)),
        .groups = "drop"
      )
    y_label <- paste0("RMSE of ", y_label_base, " Estimate")
  } else { # metric == "mae"
    error_data <- lambda_results %>%
      filter(!is.na(.data[[error_col]])) %>%
      group_by(estimator, n) %>%
      summarise(
        Metric_Value = mean(abs(.data[[error_col]]), na.rm = TRUE),
        .groups = "drop"
      )
    y_label <- paste0("MAE of ", y_label_base, " Estimate")
  }
  
  # Apply estimator filter again (for safety)
  if (!is.null(estimators)) {
    error_data <- error_data %>%
      filter(estimator %in% estimators)
  }
  
  # -----------------------------------
  # Scale metric by n^alpha (if requested)
  # -----------------------------------
  error_data <- error_data %>%
    mutate(Metric_Value = Metric_Value * n^alpha)
  
  if (alpha != 0) {
    y_label <- paste0(y_label, " × n^", alpha)
  }
  
  trim_info <- if (trim_lambda > 0) {
    paste0(" (Trim ", trim_lambda * 100, "%)")
  } else {
    ""
  }
  
  main_title <- paste0(
    "Estimator Convergence - ",
    toupper(metric), " (", y_label_base, ")", trim_info
  )
  
  # -----------------------------------
  # Base plot: estimator curves
  # -----------------------------------
  p <- ggplot(
    data = error_data,
    aes(
      x     = n,
      y     = Metric_Value,
      color = estimator,
      shape = estimator,
      group = estimator
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    labs(
      title = main_title,
      x     = "Sample Size (n)",
      y     = y_label
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # -----------------------------------
  # (Optional) log-log scaling
  # -----------------------------------
  if (log_scale) {
    p <- p +
      scale_x_log10(breaks = sort(unique(error_data$n))) +
      scale_y_log10() +
      labs(title = paste0(main_title, " (Log-Log Scale)"))
  }
  
  # -----------------------------------
  # Theoretical MLE & Bayes RMSE curves (ONLY for RMSE)
  # -----------------------------------
  if (metric == "rmse") {
    n_vals <- sort(unique(error_data$n))
    
    # Handle delta that might be scalar, function(n), or vector
    delta_n <- NULL
    
    if (is.function(delta)) {
      # delta is a function of n
      delta_n <- vapply(n_vals, delta, numeric(1))
    } else if (is.numeric(delta)) {
      if (length(delta) == 1L) {
        # scalar: same delta for all n
        delta_n <- rep(delta, length(n_vals))
      } else if (length(delta) == length(n_vals)) {
        # vector aligned to n_vals by position
        delta_n <- delta
      } else {
        warning(
          "Length of 'delta' does not match number of unique n. ",
          "Using the first element of delta for all n."
        )
        delta_n <- rep(delta[1L], length(n_vals))
      }
    } else {
      stop("delta must be numeric (scalar/vector) or a function of n.")
    }
    
    # In most settings fisher_info is a scalar; if vector, align similarly
    if (!is.numeric(fisher_info) || length(fisher_info) < 1L) {
      stop("fisher_info must be a numeric scalar or vector.")
    }
    if (length(fisher_info) == 1L) {
      fisher_n <- rep(fisher_info, length(n_vals))
    } else if (length(fisher_info) == length(n_vals)) {
      fisher_n <- fisher_info
    } else {
      warning(
        "Length of 'fisher_info' does not match number of unique n. ",
        "Using the first element for all n."
      )
      fisher_n <- rep(fisher_info[1L], length(n_vals))
    }
    
    # r = delta(n) * |I(f)|
    r_n <- delta_n * sqrt(abs(fisher_n))
    
    # Theoretical RMSE for tau
    z3 <- zeta_3_value()
    tau_rmse_mle   <- sqrt(26)       / (r_n^2)
    tau_rmse_bayes <- 4 * sqrt(z3)   / (r_n^2)
    
    theory_df <- tibble(n = n_vals)
    
    if (scale_target == "lambda") {
      # Convert tau RMSE to lambda RMSE and apply n^alpha scaling
      theory_df <- theory_df %>%
        mutate(
          Theory_MLE   = (tau_rmse_mle   / n) * n^alpha,
          Theory_Bayes = (tau_rmse_bayes / n) * n^alpha
        )
    } else { # scale_target == "tau"
      theory_df <- theory_df %>%
        mutate(
          Theory_MLE   = tau_rmse_mle   * n^alpha,
          Theory_Bayes = tau_rmse_bayes * n^alpha
        )
    }
    
    # Add dashed/dotted black theoretical curves (no legend entries)
    p <- p +
      geom_line(
        data        = theory_df,
        aes(x = n, y = Theory_MLE),
        inherit.aes = FALSE,
        linetype    = "dashed",
        linewidth   = 0.9,
        colour      = "black"
      ) +
      geom_line(
        data        = theory_df,
        aes(x = n, y = Theory_Bayes),
        inherit.aes = FALSE,
        linetype    = "dotted",
        linewidth   = 0.9,
        colour      = "black"
      )
  }
  
  print(p)
  invisible(p)
}

# ============================================================
# 2) Estimate Convergence Rate and Constant
#    via linear regression in log-log space
# ============================================================
estimate_convergence_rate <- function(results_df,
                                      true_lambda,
                                      metric      = "rmse",
                                      min_n       = NULL,
                                      max_n       = NULL,
                                      trim_lambda = 0) {
  lambda_results <- prepare_lambda_results(
    results_df  = results_df,
    true_lambda = true_lambda,
    trim_lambda = trim_lambda
  )
  
  metric <- tolower(metric)
  if (!metric %in% c("rmse", "mae")) {
    stop("Invalid metric. Choose 'rmse' or 'mae'.")
  }
  
  if (metric == "rmse") {
    error_data <- lambda_results %>%
      filter(!is.na(error)) %>%
      group_by(estimator, n) %>%
      summarise(
        Metric_Value = sqrt(mean(error^2, na.rm = TRUE)),
        .groups = "drop"
      )
  } else { # metric == "mae"
    error_data <- lambda_results %>%
      filter(!is.na(error)) %>%
      group_by(estimator, n) %>%
      summarise(
        Metric_Value = mean(abs(error), na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Filter by n if requested
  if (!is.null(min_n)) {
    error_data <- error_data %>% filter(n >= min_n)
  }
  if (!is.null(max_n)) {
    error_data <- error_data %>% filter(n <= max_n)
  }
  if (nrow(error_data) == 0) {
    stop("No observations after filtering by n and trim_lambda.")
  }
  
  # Fit separate linear models in log-log space per estimator:
  #   log(Metric_Value) ~ log(n)
  fit_one <- function(df) {
    # Need at least two distinct n
    if (n_distinct(df$n) < 2) {
      return(tibble(
        term     = c("(Intercept)", "log(n)"),
        estimate = NA_real_,
        std.error = NA_real_
      ))
    }
    # Robustify with tryCatch
    tryCatch(
      {
        tidy(lm(log(Metric_Value) ~ log(n), data = df))
      },
      error = function(e) {
        tibble(
          term     = c("(Intercept)", "log(n)"),
          estimate = NA_real_,
          std.error = NA_real_
        )
      }
    )
  }
  
  rate_models <- error_data %>%
    group_by(estimator) %>%
    group_modify(~ fit_one(.x)) %>%
    ungroup()
  
  rate_estimates <- rate_models %>%
    select(estimator, term, estimate, std.error) %>%
    pivot_wider(
      names_from  = term,
      values_from = c(estimate, std.error)
    ) %>%
    rename(
      intercept    = `estimate_(Intercept)`,
      rate         = `estimate_log(n)`,
      intercept_se = `std.error_(Intercept)`,
      rate_se      = `std.error_log(n)`
    ) %>%
    mutate(
      constant    = exp(intercept),
      constant_se = ifelse(
        !is.na(intercept_se),
        constant * intercept_se,
        NA_real_
      )
    ) %>%
    select(estimator, rate, rate_se, constant, constant_se) %>%
    arrange(rate)
  
  metric_str <- toupper(metric)
  if (trim_lambda > 0) {
    metric_str <- paste0(metric_str, " (Trim ", trim_lambda * 100, "%)")
  }
  
  cat(paste(
    "\n--- Estimated Convergence Rates (",
    metric_str,
    " ≈ C * n^rate) ---\n"
  ))
  print(as.data.frame(rate_estimates))
  rate_estimates
}

# ============================================================
# 3) Plot Trim Sensitivity of the Estimated Convergence Rate
# ============================================================
plot_trim_sensitivity <- function(results_df,
                                  true_lambda,
                                  metric   = "rmse",
                                  min_n    = NULL,
                                  max_n    = NULL,
                                  max_trim = 0.05,
                                  n_steps  = 20) {
  trim_values        <- seq(0, max_trim, length.out = n_steps)
  all_rate_estimates <- vector("list", length(trim_values))
  
  cat(paste0(
    "--- Running sensitivity analysis for ",
    toupper(metric),
    " ---\n"
  ))
  
  for (i in seq_along(trim_values)) {
    trim_val <- trim_values[i]
    
    # Suppress printed output from estimate_convergence_rate
    rate_table <- suppressMessages(
      tryCatch(
        estimate_convergence_rate(
          results_df  = results_df,
          true_lambda = true_lambda,
          metric      = metric,
          min_n       = min_n,
          max_n       = max_n,
          trim_lambda = trim_val
        ),
        error = function(e) {
          # Return empty table with expected columns if estimation fails
          tibble(
            estimator   = character(0),
            rate        = numeric(0),
            rate_se     = numeric(0),
            constant    = numeric(0),
            constant_se = numeric(0)
          )
        }
      )
    )
    
    if (nrow(rate_table) > 0) {
      rate_table$trim_lambda <- trim_val
      all_rate_estimates[[i]] <- rate_table
    } else {
      all_rate_estimates[[i]] <- tibble()
    }
  }
  
  results <- bind_rows(all_rate_estimates)
  if (nrow(results) == 0) {
    stop("No rate estimates produced for any trim values.")
  }
  
  p <- ggplot(
    results,
    aes(
      x     = trim_lambda,
      y     = rate,
      color = estimator,
      group = estimator
    )
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    labs(
      title    = paste(
        "Convergence Rate Sensitivity to Lambda Trimming (",
        toupper(metric), ")"
      ),
      x        = "Trim Value (exclude lambda_hat outside [x, 1-x])",
      y        = "Estimated Rate (slope of log(Metric) vs log(n))",
      subtitle = "Look for stability (an 'elbow') where trimming no longer changes the rate much."
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  invisible(p)
}

# ============================================================
# 4) Error Distribution for a Fixed n
# ============================================================
plot_error_distribution_for_n <- function(results_df,
                                          true_lambda,
                                          n_value,
                                          estimators  = NULL,
                                          trim_lambda = 0,
                                          bins        = 30,
                                          show_density = TRUE,
                                          facet        = TRUE) {
  lambda_results <- prepare_lambda_results(
    results_df  = results_df,
    true_lambda = true_lambda,
    trim_lambda = trim_lambda
  )
  
  df_n <- lambda_results %>%
    filter(n == n_value)
  
  if (!is.null(estimators)) {
    df_n <- df_n %>% filter(estimator %in% estimators)
  }
  if (nrow(df_n) == 0) {
    stop("No data for the requested n/estimators after trimming.")
  }
  
  p <- ggplot(
    df_n,
    aes(
      x     = error,
      fill  = estimator,
      color = estimator
    )
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins    = bins,
      alpha   = 0.35,
      position = "identity"
    ) +
    labs(
      title = paste0(
        "Distribution of Error (lambda_hat - lambda) at n = ",
        n_value
      ),
      x     = expression(lambda - hat(lambda)),
      y     = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (show_density) {
    p <- p + geom_density(alpha = 0.0, linewidth = 1)
  }
  if (facet) {
    p <- p + facet_wrap(~ estimator, scales = "free_y")
  }
  
  print(p)
  invisible(p)
}

# ============================================================
# 5) Tail Probabilities for a Fixed n over a Grid of k
#    - P(|error| > k) and P(error > k) for each estimator
# ============================================================
plot_tail_probabilities_for_n <- function(results_df,
                                          true_lambda,
                                          n_value,
                                          k_values   = seq(0, 0.1, length.out = 50),
                                          estimators = NULL,
                                          trim_lambda = 0) {
  lambda_results <- prepare_lambda_results(
    results_df  = results_df,
    true_lambda = true_lambda,
    trim_lambda = trim_lambda
  )
  
  df_n <- lambda_results %>%
    filter(n == n_value)
  
  if (!is.null(estimators)) {
    df_n <- df_n %>% filter(estimator %in% estimators)
  }
  if (nrow(df_n) == 0) {
    stop("No data for the requested n/estimators after trimming.")
  }
  
  grid <- expand.grid(
    estimator = unique(df_n$estimator),
    k         = k_values,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  tail_probs <- grid %>%
    rowwise() %>%
    mutate(
      abs_prob   = mean(abs(df_n$error[df_n$estimator == estimator]) > k, na.rm = TRUE),
      upper_prob = mean(df_n$error[df_n$estimator == estimator] > k, na.rm = TRUE)
    ) %>%
    ungroup()
  
  p1 <- ggplot(
    tail_probs,
    aes(x = k, y = abs_prob, color = estimator)
  ) +
    geom_line(linewidth = 1) +
    labs(
      title = paste0(
        "Tail: P(|lambda - lambda_hat| > k) at n = ",
        n_value
      ),
      x     = "k (absolute error threshold)",
      y     = "P(|error| > k)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(
    tail_probs,
    aes(x = k, y = upper_prob, color = estimator)
  ) +
    geom_line(linewidth = 1) +
    labs(
      title = paste0(
        "Upper Tail: P(lambda - lambda_hat > k) at n = ",
        n_value
      ),
      x     = "k (one-sided error threshold)",
      y     = expression(P(lambda - hat(lambda) > k))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p1)
  print(p2)
  
  invisible(list(
    abs_tail   = p1,
    upper_tail = p2,
    data       = tail_probs
  ))
}

# ============================================================
# 6) Robust Dispersion Metrics (IQR, MAD) for all n
# ============================================================
compute_robust_dispersion <- function(results_df,
                                      true_lambda,
                                      trim_lambda  = 0,
                                      mad_constant = 1) {
  lambda_results <- prepare_lambda_results(
    results_df  = results_df,
    true_lambda = true_lambda,
    trim_lambda = trim_lambda
  )
  
  dispersion <- lambda_results %>%
    filter(!is.na(error)) %>%
    group_by(estimator, n) %>%
    summarise(
      IQR          = IQR(error, na.rm = TRUE),
      MAD          = median(
        abs(error - median(error, na.rm = TRUE)),
        na.rm = TRUE
      ) * mad_constant,
      median_error = median(error, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    arrange(estimator, n)
  
  dispersion
}

# ============================================================
# 7) Plot Robust Dispersion vs n (IQR & MAD)
# ============================================================
plot_robust_dispersion_vs_n <- function(results_df,
                                        true_lambda,
                                        trim_lambda  = 0,
                                        mad_constant = 1,
                                        log_scale    = TRUE) {
  disp <- compute_robust_dispersion(
    results_df   = results_df,
    true_lambda  = true_lambda,
    trim_lambda  = trim_lambda,
    mad_constant = mad_constant
  )
  
  if (nrow(disp) == 0) {
    stop("No dispersion data available.")
  }
  
  disp_long <- disp %>%
    pivot_longer(
      cols      = c(IQR, MAD),
      names_to  = "metric",
      values_to = "value"
    )
  
  p <- ggplot(
    disp_long,
    aes(
      x     = n,
      y     = value,
      color = estimator,
      group = estimator
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(
      title = "Robust Dispersion vs n (IQR and MAD)",
      x     = "Sample Size (n)",
      y     = "Dispersion"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  print(p)
  invisible(list(plot = p, data = disp))
}

# ============================================================
# Example usage (adapt paths and calls to your setup)
# ============================================================
## Load your results (uncomment and adapt)
results_filename <- "results/compare_rate_estimate_t_oracle_tests_2025-11-13_16-15-25.rds"

sim_output <- readRDS(results_filename)
df          <- as.data.frame(sim_output$results)
true_lambda <- sim_output$data_params$changepoint_spec$lambda

# # Extract delta and Fisher information
delta        <- sim_output$data_params$changepoint_spec$delta
noise_params <- sim_output$data_params$noise_params
noise_dist   <- sim_output$data_params$noise_dist
fisher_info  <- get_stein_fisher_info(noise_dist, noise_params)
# 
# # Basic convergence plot (RMSE, with theoretical curves)
plot_convergence_metric(
  results_df   = df,
  true_lambda  = true_lambda,
  metric       = "rmse",
  scale_target = "lambda",
  log_scale    = TRUE,
  trim_lambda  = 0.00,
  alpha        = 0.5,
  delta        = delta,
  fisher_info  = fisher_info,
  estimators = c("bayes", "oracle", "score_test_fisher")
)
# 
# # Basic convergence plot (MAE, no theory curves):
# # plot_convergence_metric(
# #   df, true_lambda, metric = "mae",
# #   scale_target = "lambda", log_scale = TRUE,
# #   trim_lambda  = 0.00,
# #   alpha        = 0,
# #   delta        = delta,
# #   fisher_info  = fisher_info
# # )
# 
# # Estimate rates for a range of n
# # estimate_convergence_rate(
# #   df, true_lambda,
# #   metric      = "rmse",
# #   min_n       = 500,
# #   max_n       = 2000,
# #   trim_lambda = 0.05
# # )
# 
# # Trim sensitivity
# # plot_trim_sensitivity(
# #   df, true_lambda,
# #   metric   = "rmse",
# #   min_n    = 500,
# #   max_n    = 2000,
# #   max_trim = 0.1,
# #   n_steps  = 20
# # )
# 
# # Fixed-n distribution & tails:
# # plot_error_distribution_for_n(
# #   df, true_lambda,
# #   n_value    = 1000,
# #   trim_lambda = 0.05
# # )
# # plot_tail_probabilities_for_n(
# #   df, true_lambda,
# #   n_value    = 1000,
# #   k_values   = seq(-0.05, 0.05, length.out = 40),
# #   trim_lambda = 0.05
# # )
# 
# # Robust dispersion across n:
# # disp <- compute_robust_dispersion(
# #   df, true_lambda,
# #   trim_lambda = 0.05
# # )
# # plot_robust_dispersion_vs_n(
# #   df, true_lambda,
# #   trim_lambda = 0.05,
# #   log_scale   = TRUE
# # )
