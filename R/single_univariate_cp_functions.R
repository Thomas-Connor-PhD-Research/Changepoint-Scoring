#' Generate a uni-variate time series with a single mean-shift change-point (CP)
#' 
#' @param n Integer. Length of time series.
#' @param tau Integer. Location of the CP (last point before the change).
#'  1 <= tau <= n (if tau == n, then the series does not contain a CP)
#' @param delta Numeric. Size of the mean shift (may be +ve or -ve)
#' @param noise_dist  String. Name of the family of the noise distribution 
#'  (see R/noise_distributions.R for complete list).
#' @param noise_params List. Family specific noise distribution parameters.
#' 
#' @return List containing time series 'Y' and true underlying noise 'noise'.

generate_timeseries <- function(n, tau, delta, noise_dist, noise_params, noise_mu){
  noise <- sample_from_distribution(n, noise_dist, noise_params)
  
  # Start from pure noise on noise_mu, then add a shift after tau
  Y <- noise + noise_mu
  if (tau < n){
    Y[(tau+1L):n] <- Y[(tau+1L):n] + delta
  }
  
  return(list(Y = Y, noise = noise))
}

# ------------------------------------------------------------------------------
# CUSUM-STYLE ESTIMATORS
# ------------------------------------------------------------------------------

#' Detect CP from (normalised) score estimate
#' 
#' Find the CP by maximising the CUSUM test statistic based on the transformation
#' function 'rho' (for accurate estimation of 'delta', 'rho' should be normalised
#' to have variance 1 on the noise (or noise estimates))
#' 
#' @param Y. Numeric vector. Univariate time series.
#' @param rho. Function. Normalised score estimate.
#' @param return_path. Bool. If true, returns 
#' 
#' @return List with estimate of tau 'tau_hat', delta 'delta_hat'
#' and the maximum test statistic 'test_stat_max'
#' 
detect_cp_from_rho <- function(Y, rho, mu_hat,
                               return_delta = TRUE,
                               return_max = TRUE, 
                               return_path = FALSE){
  n <- length(Y)
  k <- 1:(n-1)
  
  # Calculate delta estimates at each t
  centered_rho_y <- cumsum(rho(Y-mu_hat) - mean(rho(Y-mu_hat)))
  delta_scaling <- n/(k * (n-k))
  possible_deltas <- delta_scaling*centered_rho_y[k]
  
  # Generate vector of test statistics
  constant_variance_scaling <- sqrt((k * (n-k))/n)
  test_statistics <- (constant_variance_scaling * possible_deltas)^2
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- possible_deltas[tau_hat]
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  
  
  return(results)
}

#' Calculate the CUSUM estimator (rho(x) = x, mu_hat = N/A [since rho linear])
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @param return_path. Bool. If true, returns test statistics too.
#' @return List with CP estimates.
#' 
calculate_cusum_estimator <- function(Y, ...){
  rho <- function(x){-x}
  other_params <- list(...)
  return(do.call("detect_cp_from_rho", c(list(Y = Y, 
                                              rho = rho,
                                              mu_hat = 0), 
                                         other_params)))
}

# We need to implement some variance scaling for rho
# calculate_trimmed_cusum_estimator <- function(Y, mu_hat, trim, ...){
#   neg_id <- function(x){-x}
#   C <- qnorm(1-trim/2)
#   rho <- trim_function(neg_id, C)
#   other_params <- list(...)
# return(do.call("detect_cp_from_rho", c(list(Y = Y, 
#                                             rho = rho,
#                                             mu_hat = mu_hat), other_params)))
# }

#' Calculate the oracle estimator
#' 
#' Use true noise distribution to build optimal 'rho' based off known score 
#' and known Fisher info
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @param noise_dist  String. Name of the family of the noise distribution 
#'  (see R/noise_distributions.R for complete list).
#' @param noise_params List. Family specific noise distribution parameters.
#' @param return_path. Bool. If true, returns test statistics too.
#' @return List with CP estimates.
calculate_oracle_estimator <- function(Y, noise_dist, noise_params, noise_mu,
                                       ...){
  # Get oracle rho
  oracle_score_fn <- create_score_function(noise_dist, noise_params)
  
  true_fi <- get_stein_fisher_info(noise_dist, noise_params)
  
  if (!is.na(true_fi)){
    oracle_fi <- true_fi
  } else{
    oracle_fi <- var(oracle_score_fn(Y))
  }
    
  rho <- function(x){oracle_score_fn(x)/oracle_fi}
  
  other_params <- list(...)
  
  return(do.call("detect_cp_from_rho", c(list(Y = Y, 
                                              rho = rho,
                                              mu = noise_mu
                                              ), other_params)))
  oracle_score_fn <- create_score_function(noise_dist, noise_params)
  
  
}


calculate_score_test_estimator <- function(Y, noise_dist, noise_params,
                                           use_fisher_update,
                                           return_delta = TRUE,
                                           return_max = TRUE,
                                           return_path = FALSE
                                           ){
  # Get oracle rho
  score_fn <- create_score_function(noise_dist, noise_params)

  fi <- get_stein_fisher_info(noise_dist, noise_params)


  # global_loc <- median(Y)

  if (use_fisher_update){
    global_loc <- median(Y) - mean(score_fn(Y - median(Y)))* (1/fi) #One step-Fisher update
  } else {
    global_loc <- median(Y)
  }
  
  n <- length(Y)
  k <- 1:(n-1)

  # Calculate delta estimates at each t
  S <- cumsum(score_fn(Y-global_loc) - mean(score_fn(Y-global_loc)))
  test_statistics <- (n / (k * (n-k))) * (1 / fi) * (S[k]^2)

  possible_deltas <- (n / (k * (n-k))) * (S[k] / fi)
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- possible_deltas[tau_hat]

  results <- list(tau = tau_hat)

  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }




return(results)}


calculate_oracle_bayes_estimator <- function(Y, noise_dist, noise_params,
                                           use_fisher_update,
                                           return_delta = TRUE,
                                           return_max = TRUE,
                                           return_path = FALSE
){
  # Get oracle rho
  score_fn <- create_score_function(noise_dist, noise_params)
  
  fi <- get_stein_fisher_info(noise_dist, noise_params)
  
  
  # global_loc <- median(Y)
  
  if (use_fisher_update){
    global_loc <- median(Y) - mean(score_fn(Y - median(Y)))* (1/fi) #One step-Fisher update
  } else {
    global_loc <- median(Y)
  }
  
  n <- length(Y)
  k <- 1:(n-1)
  
  w = k*(n-k)/n
  s = cumsum(score_fn(Y-global_loc) - mean(score_fn(Y-global_loc)))
  
  m = (fi*w[k])^(-0.5) * exp(s[k]^2 / (2 * fi * w))
  
  tau_hat <- sum(k*m)/sum(m)
  
  results <- list(tau = round(tau_hat))
  
  # Return desired estimators
  if (return_delta){
    results$delta <- 0
  }
  if (return_max){
    results$max <- 0
  }
  if (return_path) {
    results$path <- 0
  }
  
  
  
  
  return(results)}



#' Calculate estimator using a score estimated from KNOWN noise
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @param noise. Numeric vector. True underlying noise
#' @param scoring_method. String. Method used in 'score_estimation' fn
#' @param return_path. Bool. If true, returns test statistics too.
#' @return List with CP estimates.
calculate_scoring_known_noise_estimator <- function(Y, noise, scoring_method, ...){
  # Get rho estimate from noise
  score_fn_estimate <- score_estimation(noise, scoring_method)[[scoring_method]]
  fi_estimate <- var(score_fn_estimate(noise))
  rho <- function(x){score_fn_estimate(x)/fi_estimate}
  
  other_params <- list(...)
  
  return(do.call("detect_cp_from_rho", c(list(Y = Y, 
                                              rho = rho, 
                                              mu = 0), # Mu is estimated implicitly in the rho via scoring
                                         other_params)))
}



calculate_scoring_estimator <- function(Y, 
                                        scoring_method, 
                                        iterations = 1, 
                                        rho_initial = NULL,
                                        ...){
  other_params <- list(...)
  
  n <- length(Y)
  if (!is.null(rho_initial)){
    rho <- rho_initial
  } else{
    # Get score estimate from Y (which potentially includes CP)
    score_fn_estimate <- score_estimation(Y, scoring_method)[[scoring_method]]
    fi_estimate <- var(score_fn_estimate(Y))
    rho <- function(x){score_fn_estimate(x)/fi_estimate}
  }
  noise_estimate <- numeric(n)
  
  i <- 0
  
  while (i < iterations){
    int_results <- detect_cp_from_rho(Y, rho)
    
    tau_hat <- int_results$tau
    delta_hat <- int_results$delta
    
    # Update noise estimate
    noise_estimate[1:tau_hat] <- Y[1:tau_hat]
    noise_estimate[(tau_hat+1):n] <- Y[(tau_hat+1):n] - delta_hat
    
    # Update score estimate
    score_fn_estimate <- score_estimation(noise_estimate, scoring_method)[[scoring_method]]
    fi_estimate <- var(score_fn_estimate(noise_estimate))
    rho <- function(x){score_fn_estimate(x)/fi_estimate}
  
    i <- i + 1
  }
  return(do.call("detect_cp_from_rho", c(list(Y = Y, 
                                              rho = rho,
                                              mu = 0), # Mu is estimated implicitly in the rho via scoring
                                         other_params)))
}

# ------------------------------------------------------------------------------
# NON-PARAMETRIC ESTIMATION METHODS
# ------------------------------------------------------------------------------

calculate_wilcoxon_estimator <- function(Y, return_delta = TRUE,
                                         return_max = TRUE, return_path = FALSE){
  n <- length(Y)
  k <- 1:(n-1)
  
  # Convert data to ranks
  ranks <- rank(Y)[k]
  
  # Sum or ranks (centered so that min at 0)
  sum_of_ranks <- cumsum(ranks) - k*(k+1)/2 

  # Centre and variance scale under H0 (no change)
  centered_sum_of_ranks <- sum_of_ranks  - k*(n-k)/2
  variance_scaling <- sqrt(k*(n-k)*(n+1) / 12)
  
  test_statistics <- abs(centered_sum_of_ranks/variance_scaling)
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- median(Y[(tau_hat + 1):n]) - median(Y[1:tau_hat]) # Robust
  
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  return(results)
}
  

calculate_sign_estimator <- function(Y, return_delta = TRUE,
                                       return_max = TRUE, return_path = FALSE){
  n <- length(Y)
  k <- 1:(n-1)
  
  centre <- median(Y)
  
  signs <- sign(Y - centre)
  cumsum_of_signs <- cumsum(signs - mean(signs))
  
  test_statistics <- abs(cumsum_of_signs[k])
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- median(Y[(tau_hat + 1):n]) - median(Y[1:tau_hat]) # Robust
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  return(results)
}

calculate_oracle_sign_estimator <- function(Y, noise_mu, 
                                            delta,
                                            return_delta = TRUE,
                                            return_max = TRUE,
                                            return_path = FALSE){
  n <- length(Y)
  k <- 1:(n-1)
  
  centre <- noise_mu + delta/2
  
  signs <- sign(Y - centre)
  cumsum_of_signs <- cumsum(signs - mean(signs))
  
  test_statistics <- abs(cumsum_of_signs[k])
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- median(Y[(tau_hat + 1):n]) - median(Y[1:tau_hat]) # Robust
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  return(results)
}

calculate_oracle_trimmed_cusum_estimator <- function(Y, trim,
                                                     noise_dist, 
                                                     noise_params,
                                                     noise_mu, 
                                            delta,
                                            return_delta = TRUE,
                                            return_max = TRUE,
                                            return_path = FALSE){
  
  noise_variance <- get_variance(noise_dist, noise_params)
  
  n <- length(Y)
  k <- 1:(n-1)
  
  neg_id <- function(x){-x}
  C <- noise_variance*qnorm(1-trim/2) + abs(delta) # Expected Trim% of 
                                                   # observations in each regime are "Winsorized"
  rho <- trim_function(neg_id, C)
  
  centre <- noise_mu + delta/2
  
  centered_rho_y <- cumsum(rho(Y-centre) - mean(rho(Y-centre)))
  
  test_statistics <- abs(centered_rho_y[k])
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- median(Y[(tau_hat + 1):n]) - median(Y[1:tau_hat]) # Robust
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    results$delta <- delta_hat
  }
  if (return_max){
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  return(results)
}
  
calculate_cusum_chosen_variance_estimator <- function(Y, alpha,
                                                      return_delta = TRUE,
                                                      return_max = TRUE,
                                                      return_path = FALSE){

  n <- length(Y)
  k <- 1:(n-1)
  
  neg_id <- function(x){-x}
  
  centered_rho_y <- cumsum(neg_id(Y) - mean(neg_id(Y)))
  
  variance_scaling <- (k*(1-k/n))^(-alpha/2)
  
  test_statistics <- abs(variance_scaling*centered_rho_y[k])
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  
  results <- list(tau = tau_hat)
  
  # Return desired estimators
  if (return_delta){
    delta_hat <- median(Y[(tau_hat + 1):n]) - median(Y[1:tau_hat]) # Robust
    results$delta <- delta_hat
  }
  if (return_max){
    test_stat_max <- test_statistics[tau_hat]
    results$max <- test_stat_max
  }
  if (return_path) {
    results$path <- test_statistics
  }
  return(results)
}
                                            
# ------------------------------------------------------------------------------
# SIMULATION FUNCTIONS
# ------------------------------------------------------------------------------

#' Run a single simulation for univariate single CP detection
#' 
#' @param data_params List. Parameters for 'generate_timeseries'
#' @param estimator_params Named list. List where each element is an estimator
#' to run, and contains its own params.
#' 
#' @return List containing results from all estimators 
#'  (eg: "cusum.tau_hat", "oracle.delta_hat")
                                        
simulate_single_changepoint <- function(data_params, estimator_params) {
  # Generate time series data
  ts_data <- generate_timeseries(
    n = data_params$n,
    tau = data_params$changepoint_spec$tau,
    delta = data_params$changepoint_spec$delta,
    noise_dist = data_params$noise_dist,
    noise_params = data_params$noise_params
  )
  Y <- ts_data$Y
  noise <- ts_data$noise
  
  results_list <- list()
  
  # Loop through all requested estimators
  for (estimator_name in names(estimator_params)) {
    params <- estimator_params[[estimator_name]]
    
    result <- switch(
      estimator_name,
      
      "cusum" = calculate_cusum_estimator(Y),
      
      "cusum_trimmed" = calculate_trimmed_cusum_estimator(Y, trim),
      
      "oracle" = calculate_oracle_estimator(
        Y,
        data_params$noise_dist,
        data_params$noise_params
      ),
      
      "scoring_known_noise" = calculate_scoring_known_noise_estimator(
        Y,
        noise,
        params$scoring_method
      ),
      
      "scoring" = calculate_scoring_estimator(
        Y,
        params$scoring_method,
        params$iterations
      ),
      "iterative_scoring" = calculate_scoring_estimator(
        Y,
        params$scoring_method,
        params$iterations
      ),
      "scoring_initial_rho" = calculate_scoring_estimator(
        Y,
        params$scoring_method,
        params$iterations,
        params$rho_initial
      ),
      "sign_of_medians" = calculate_median_estimator(
        Y
      ),
      
      "wilcoxon_ranks" = calculate_wilcoxon_estimator(
        Y),
      
      "ECDF" = calculate_ecdf_estimator(
        Y
      ),
      
      
     # Add estimators (if needed)
      
      # Stop if invalid estimator added
      stop(paste("Unknown estimator in estimator_params:", estimator_name))
    )
    
    # Name the results vector (e.g., "cusum.tau_hat", "cusum.delta_hat")
    # names(result) <- paste(estimator_name, names(result), sep = ".")
    results_list[[estimator_name]] <- result
  }
  
  # Return a single flat vector of all results
  return(unlist(results_list))
}
  
#' Run a single simulation for multiple sample sizes 'n' (No Burn-in)
#'
#' This function iterates over a vector of 'n' values, generates a dataset
#' for each 'n', and runs a suite of estimators. It returns a single named
#' list where each result is suffixed with its corresponding 'n' value.
#' The changepoint search is conducted over the full range 1:(n-1).
#'
#' @param data_params List. Parameters for 'generate_timeseries'.
#'   `data_params$n` MUST be a vector of sample sizes (e.g., c(200, 500, 1000)).
#'   `data_params$changepoint_spec$tau` MUST be a fraction (0 < tau < 1).
#' @param estimator_params Named list. List of estimators to run.
#'
#' @return A single named list of all results, e.g., "cusum.tau.200", "oracle.tau.500"

simulate_single_changepoint_rate_estimation <- function(data_params, estimator_params) {
  
  all_n_results_list <- list()
  
  # --- 1. Outer loop over the vector of n values ---
  for (current_n in data_params$n) {
    current_data_params <- data_params
    
    current_data_params$n <- current_n
    current_data_params$changepoint_spec$tau <- floor(data_params$changepoint_spec$lambda * current_n)
    
    delta_spec <- data_params$changepoint_spec$delta
    
    if (is.function(delta_spec)) {
      # If it's a function, call it with current_n
      current_data_params$changepoint_spec$delta <- delta_spec(current_n)
    } else {
      # Otherwise, use it as a static value (backward-compatible)
      current_data_params$changepoint_spec$delta <- delta_spec
    }
    
    
    ts_data <- generate_timeseries(
      n = current_data_params$n,
      tau = current_data_params$changepoint_spec$tau,
      delta = current_data_params$changepoint_spec$delta,
      noise_dist = current_data_params$noise_dist,
      noise_params = current_data_params$noise_params,
      noise_mu = current_data_params$noise_mu
    )

    Y <- ts_data$Y
    noise <- ts_data$noise
    
    results_list_for_n <- list()
    
    # --- 2. Loop through all requested estimators ---
    for (estimator_name in names(estimator_params)) {
      params <- estimator_params[[estimator_name]]
      estimator_type <- params$estimator_type

      # 3. Use switch to call the correct function
      estimator_func_name <- switch(
        estimator_type, 
        "cusum" = "calculate_cusum_estimator",
        "oracle" = "calculate_oracle_estimator",
        "score_test" = "calculate_score_test_estimator",
        "oracle_bayes" = "calculate_oracle_bayes_estimator",
        "scoring_known_noise" = "calculate_scoring_known_noise_estimator",
        "scoring" = "calculate_scoring_estimator",
        "iterative_scoring" = "calculate_scoring_estimator",
        "scoring_initial_rho" = "calculate_scoring_estimator",
        "wilcoxon_ranks" = "calculate_wilcoxon_estimator",
        "ECDF" = "calculate_ecdf_estimator",
        "sign_of_medians" = "calculate_sign_estimator",
        "oracle_sign" = "calculate_oracle_sign_estimator",
        "oracle_trimmed_cusum" = "calculate_oracle_trimmed_cusum_estimator",
        "cusum_chosen_variance" = "calculate_cusum_chosen_variance_estimator",
        
        stop(paste("Unknown estimator:", estimator_type))
      )
      
      # 4. Build the argument list for do.call
      args_list <- list(
        Y = Y,
        return_delta = FALSE,
        return_max = FALSE
      )
      
      # Add estimator-specific arguments
      if (estimator_type == "oracle") {
        args_list$noise_dist <- current_data_params$noise_dist
        args_list$noise_params <- current_data_params$noise_params
        args_list$noise_mu <- current_data_params$noise_mu
      } else if (estimator_type == "scoring_known_noise") {
        args_list$noise <- noise
        args_list$scoring_method <- params$scoring_method
      } else if (grepl("scoring", estimator_type)) {
        args_list$scoring_method <- params$scoring_method
        args_list$iterations <- params$iterations
        if (estimator_type == "scoring_initial_rho") {
          args_list$rho_initial <- params$rho_initial
        }
      }else if (estimator_type == "oracle_sign"){
        args_list$noise_mu <- current_data_params$noise_mu
        args_list$delta <- current_data_params$changepoint_spec$delta
      } else if (estimator_type == "oracle_trimmed_cusum"){
        args_list$trim <- params$trim
        args_list$noise_dist <- current_data_params$noise_dist
        args_list$noise_params <- current_data_params$noise_params
        args_list$noise_mu <- current_data_params$noise_mu
        args_list$delta <- current_data_params$changepoint_spec$delta
      } else if (estimator_type == "cusum_chosen_variance"){
        args_list$alpha <- params$alpha
      } else if (estimator_type == "score_test"){
        args_list$use_fisher_update <- params$use_fisher_update
        args_list$noise_dist <- current_data_params$noise_dist
        args_list$noise_params <- current_data_params$noise_params
      } else if (estimator_type == "oracle_bayes"){
        args_list$use_fisher_update <- params$use_fisher_update
        args_list$noise_dist <- current_data_params$noise_dist
        args_list$noise_params <- current_data_params$noise_params
      }
      
      
      # 5. Call the estimator function
      result <- do.call(estimator_func_name, args_list)
      
      # 6. Name the results with the .n suffix
      names(result) <- paste(estimator_name, names(result), current_n, sep = ".")
      
      results_list_for_n <- c(results_list_for_n, result)
    } # end estimator loop
    
    all_n_results_list <- c(all_n_results_list, results_list_for_n)
    
  } # end n loop
  
  # 7. Return a single flat list of all results
  return(all_n_results_list)
}