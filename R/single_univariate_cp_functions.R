#' Generate a uni-variate time series with a single mean-shift change-point (CP)
#' 
#' @param n Integer. Length of time series.
#' @param tau Integer. Location of the CP (last point before the change).
#'  1 <= tau <= n-1
#' @param delta Numeric. Size of the mean shift (may be +ve or -ve)
#' @param noise_dist  String. Name of the family of the noise distribution 
#'  (see R/noise_distributions.R for complete list).
#' @param noise_params List. Family specific noise distribution parameters.
#' 
#' @return List containing time series 'Y' and true underlying noise 'noise'.
generate_timeseries <- function(n, tau, delta, noise_dist, noise_params){
  noise <- sample_from_distribution(n, noise_dist, noise_params)
  
  Y <- numeric(n)
  Y[1:tau] <- noise[1:tau]
  Y[(tau+1):n] <- noise[(tau+1):n] + delta
  
  return(list(Y = Y, noise = noise))
}

#' Detect CP from (normalised) score estimate
#' 
#' Find the CP by maximising the CUSUM test statistic based on the transformation
#' function 'rho' (for accurate estimation of 'delta', 'rho' should be normalised
#' to have variance 1 on the noise (or noise estimates))
#' 
#' @param Y. Numeric vector. Univariate time series.
#' @param rho. Function. Normalised score estimate.
#' 
#' @return List with estimate of tau 'tau_hat', delta 'delta_hat'
#' and the maximum test statistic 'test_stat_max'
#' 
detect_cp_from_rho <- function(Y, rho){
  n <- length(Y)
  k <- 1:(n-1)
  
  # Calculate delta estimates at each t
  centered_rho_y <- -cumsum(rho(Y) - mean(rho(Y)))
  delta_scaling <- n/(k * (n-k))
  possible_deltas <- delta_scaling*centered_rho_y[k]
  
  # Generate vector of test statistics
  constant_variance_scaling <- sqrt((k * (n-k))/n)
  test_statistics <- (constant_variance_scaling * possible_deltas)^2
  
  # Obtain estimates
  tau_hat <- which.max(test_statistics)
  test_stat_max <- test_statistics[tau_hat]
  delta_hat <- possible_deltas[tau_hat]
  
  return(list(tau = tau_hat, 
           delta = delta_hat, 
           max = test_stat_max
           ))
}

#' Calculate the CUSUM estimator (rho(x) = x)
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @return List with CP estimates.
#' 
calculate_cusum_estimator <- function(Y){
  rho <- function(x){x}
  return(detect_cp_from_rho(Y, rho))
}

#' Calculate the oracle estimator
#' 
#' Use true noise distribution to build optimal 'rho' based off known score 
#' and known Fisher info
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @param noise_dist  String. Name of the family of the noise distribution 
#'  (see R/noise_distributions.R for complete list).
#' @param noise_params List. Family specific noise distribution parameters.
#' @return List with CP estimates.
calculate_oracle_estimator <- function(Y, noise_dist, noise_params){
  # Get oracle rho
  oracle_score_fn <- create_score_function(noise_dist, noise_params)
  oracle_fi <- get_stein_fisher_info(noise_dist, noise_params)
  rho <- function(x){oracle_score_fn(x)/oracle_fi}
  
  return(detect_cp_from_rho(Y, rho))
}

#' Calculate estimator using a score estimated from KNOWN noise
#' 
#' @param Y Numeric vector. Univariate timeseries.
#' @param noise. Numeric vector. True underlying noise
#' @param scoring_method. String. Method used in 'score_estimation' fn
#' @return List with CP estimates.
calculate_scoring_known_noise_estimator <- function(Y, noise, scoring_method){
  # Get rho estimate from noise
  score_fn_estimate <- score_estimation(noise, scoring_method)
  fi_estimate <- var(score_fn_estimate(noise))
  rho <- function(x){score_fn_estimate(x)/fi_estimate}
  
  return(detect_cp_from_rho(Y, rho))
}


calculate_scoring_estimator <- function(Y, 
                                        scoring_method, 
                                        iterations = 1, 
                                        rho_initial = NULL){
  n <- length(Y)
  
  if (!is.null(rho_initial)){
    rho <- rho_initial
  } else{
    # Get score estimate from Y (which potentially includesa CP)
    score_fn_estimate <- score_estimation(Y, scoring_method)
    fi_estimate <- var(score_fn_estimate(Y))
    rho <- function(x){score_fn_estimate(x)/fi_estimate}
  }
  
  noise_estimate <- numeric(n)
  
  for (i in 0:iterations){
    cp_results <- detect_cp_from_rho(Y, rho)
    
    # If not final iteration, update rho
    if (i < iterations){
      tau_hat <- cp_results$tau
      delta_hat <- cp_results$delta
      
      # Update noise estimate
      noise_estimate[1:tau_hat] <- Y[1:tau_hat]
      noise_estimate[(tau_hat+1):n] <- Y[(tau_hat+1):n] - delta_hat
      
      # Update score estimate
      score_fn_estimate <- score_estimation(noise_estimate, scoring_method)
      fi_estimate <- var(score_fn_estimate(noise_estimate))
      rho <- function(x){score_fn_estimate(x)/fi_estimate}
      
    }
  }
  
  return(cp_results)
}

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
  print(Y[1])
  results_list <- list()
  
  # Loop through all requested estimators
  for (estimator_name in names(estimator_params)) {
    params <- estimator_params[[estimator_name]]
    
    result <- switch(
      estimator_name,
      
      "cusum" = calculate_cusum_estimator(Y),
      
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
  
  