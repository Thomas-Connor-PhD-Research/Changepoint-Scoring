# ---- ESTIMATOR CALCULATION ----
calculate_cusum_estimator <- function(Y, params) {
  transform_fn <- function(y) y
  test_stat <- generate_test_statistic(Y, transform_fn, params$burnin, params$n)
  
  results <- get_max(test_stat, params$burnin)
  results$delta <- mean(Y[(results$tau + 1):params$n]) - mean(Y[1:results$tau])
  return(results)
}

calculate_oracle_estimator <- function(Y, params) {
  oracle_score_fn <- create_score_function(params$noise_dist, params$noise_params)
  test_stat <- generate_test_statistic(Y, oracle_score_fn, params$burnin, params$n)
  
  results <- get_max(test_stat, params$burnin)
  results$delta <- mean(Y[(results$tau + 1):params$n]) - mean(Y[1:results$tau])
  return(results)
}

calculate_scoring_estimator <- function(Y, params, tau_hat = NULL) {
  if (is.null(tau_hat)) {
    tau_hat <- calculate_cusum_estimator(Y, params)$tau
  }
  
  noise_est <- Y
  noise_est[(tau_hat + 1):params$n] <- noise_est[(tau_hat + 1):params$n] - (mean(Y[(tau_hat + 1):params$n]) - mean(Y[1:tau_hat]))
  scoring_score_fn <- spline_score(noise_est, df = cv_spline_score(noise_est)$df_min)$rho
  test_stat <- generate_test_statistic(Y, scoring_score_fn, params$burnin, params$n)
  
  results <- get_max(test_stat, params$burnin)
  results$delta <- mean(Y[(results$tau + 1):params$n]) - mean(Y[1:results$tau])
  return(results)
}


calculate_scoring_known_noise_estimator <- function(Y, noise, params){
  scoring_score_fn <- spline_score(noise, df = cv_spline_score(noise)$df_min)$rho
  test_stat <- generate_test_statistic(Y, scoring_score_fn, params$burnin, params$n)
  
  results <- get_max(test_stat, params$burnin)
  results$delta <- mean(Y[(results$tau + 1):params$n]) - mean(Y[1:results$tau])
  return(results)
}

calculate_iterative_scoring_estimator <- function(Y, params, tau_hat = NULL){
  results_initial <- calculate_scoring_estimator(Y, params, tau_hat)
  results <- calculate_scoring_estimator(Y, params, results_initial$tau)
  results$delta <- mean(Y[(results$tau + 1):params$n]) - mean(Y[1:results$tau])
  return(results)
}


# TO BE IMPLEMENTED
calculate_scoring_sample_split_estimator <- function(Y, noise, params){
  odd_indicies <- seq(1, params$n, by=2)
  even_indicies <- seq(1, params$n, by=2)
  
  Y_odd <- Y[odd_indicies]
  Y_even <- Y[even_indicies]
  

  tau_hat <- calculate_cusum_estimator(Y, params)$tau
  noise_est <- Y
  noise_est[(tau_hat + 1):params$n] <- noise_est[(tau_hat + 1):params$n] - (mean(Y[(tau_hat + 1):params$n]) - mean(Y[1:tau_hat]))
}

# TO BE IMPLEMENTED
calculate_trimmed_scoring_estimator <- function(Y, params, tau_hat){
  if (is.null(tau_hat)) {
    tau_hat <- calculate_cusum_estimator(Y, params)$tau
  }
  
  noise_est <- Y
  noise_est[(tau_hat + 1):params$n] <- noise_est[(tau_hat + 1):params$n] - (mean(Y[(tau_hat + 1):params$n]) - mean(Y[1:tau_hat]))
  scoring_score_fn <- spline_score(noise_est, df = cv_spline_score(noise_est)$df_min)$rho
  test_stat <- generate_test_statistic(Y, scoring_score_fn, params$burnin, params$n)
  
  results <- get_max(test_stat, params$burnin)
  return(results)
}

# Generates generic generalised-CUSUM test statistic
generate_test_statistic <- function(Y, transform_fn, burnin, n) {
  transformed_Y <- transform_fn(Y)
  denom <- ((burnin + 1):(n - burnin)) * (1 - ((burnin + 1):(n - burnin)) / n)
  unscaled_stat <- cumsum(transformed_Y - mean(transformed_Y))[(burnin + 1):(n - burnin)]
  return(unscaled_stat^2 / denom)
}

get_max <- function(teststat, burnin) {
  list(tau = which.max(teststat) + burnin, max = max(teststat))
}


# ---- SIMULATION FUNCTION ----
simulate_single_changepoint <- function(params) {
  # Generate observations 
  noise <- sample_from_distribution(params$n, params$noise_dist, params$noise_params)
  Y <- noise
  Y[(params$changepoint_spec$tau + 1):params$n] <- Y[(params$changepoint_spec$tau + 1):params$n] + params$changepoint_spec$delta
  
  # Loop through and run the requested estimators
  # NOTE: Order should obey dependencies
  results_list <- list()
  if ("cusum" %in% params$estimators){
    results_list[["cusum"]] <- calculate_cusum_estimator(Y, params)
  } 
  
  if ("oracle" %in% params$estimators){
    results_list[["oracle"]] <- calculate_oracle_estimator(Y, params)
  }
  
  if ("scoring_known_noise" %in% params$estimators){
    results_list[["scoring_known_noise"]] <- calculate_scoring_known_noise_estimator(Y, noise, params)
  }
  
  
  if ("scoring" %in% params$estimators){
    tau_hat <- results_list[["cusum"]]$tau # NULL if does to exist

    results_list[["scoring"]] <- calculate_scoring_estimator(Y, params, tau_hat)
  }
  
  if ("iterative_scoring" %in% params$estimators){
    tau_hat <- results_list[["cusum"]]$tau # NULL if does to exist
    
    results_list[["iterative_scoring"]] <- calculate_iterative_scoring_estimator(Y, params, tau_hat)
  }
  
  
  
  # Combine results from all estimators into single vector
  return(unlist(results_list))
}

simulate_multiple_changepoints <- function(params) {
  message("Multiple changepoint simulation not yet implemented.")
  return(NULL)
}