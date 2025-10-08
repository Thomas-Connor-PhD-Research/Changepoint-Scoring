# Parameters
n <- 1000
delta <- 0.3
tau <- floor(n/2)
nreps <- 1000
burnin <- floor(n / 25)

sigma <- NULL
scale <- 0.5
noise_dist <- "t"
df <- 3

# Select which statistics to compute
# Choose a subset of c("cusum", "oracle", "oracle_score_est", "scoring")
# cusum = traditional cusum estimator
# oracle = scoring with oracle score fn for the noise dist
# oracle_score_est = scoring with score fn estimated from noise
# scoring = scoring with score fn estimated from estimated noise
stats_to_compute <- c("cusum", "scoring", "oracle")

source("noise_distributions.R")
source("score_estimation2.R")
oracle_score <- create_score_function(noise_dist, sd = sigma, 
                                      df = df, scale = scale)

# Load required packages
library(foreach)
library(progressr)
library(future)
library(doFuture)
library(doRNG)



# Set up future backend for parallel + progress
plan(multisession, workers = parallel::detectCores() - 1)
registerDoFuture()
handlers(global = TRUE)

## HELPER FUNCTIONS
get_max <- function(teststat) {
  list(
    max = max(teststat),
    tau_hat = which.max(teststat) + burnin
  )
}
  
generate_cusum <- function(Y){
  denom <- ((burnin + 1):(n - burnin)) * (1 - ((burnin + 1):(n - burnin)) / n)
  unscaled_cusum <- cumsum(Y - mean(Y))[(burnin + 1):(n - burnin)]
  return (unscaled_cusum^2 / denom)
}

generate_oracle <- function(Y, oracle_score){
  denom <- ((burnin + 1):(n - burnin)) * (1 - ((burnin + 1):(n - burnin)) / n)
  unscaled_oracle <- cumsum(oracle_score(Y) - mean(oracle_score(Y)))[(burnin + 1):(n - burnin)]
  return(unscaled_oracle^2 / denom)
}

generate_oracle_score_est <- function(Y, noise){
  denom <- ((burnin + 1):(n - burnin)) * (1 - ((burnin + 1):(n - burnin)) / n)
  score_est <- spline_score(noise, df = cv_spline_score(noise)$df_min)$rho
  unscaled_oracle_score_est <- cumsum(score_est(Y) - mean(score_est(Y)))[(burnin + 1):(n - burnin)]
  return(unscaled_oracle_score_est^2 / denom)
}


generate_scoring <- function(Y, tau_est = NULL){
  if (is.null(tau_est)) {
    tau_est <- get_max(generate_cusum(Y))$tau_hat
  } 
  
  noise_est <- Y
  noise_est[1:tau_est] <- noise_est[1:tau_est] 
  noise_est[(tau_est+1):n] <- noise_est[(tau_est+1):n] - (mean(noise_est[(tau_est+1):n]) - mean(noise_est[1:tau_est]))
  
  score_est <- spline_score(noise_est, df = cv_spline_score(noise_est)$df_min)$rho
  
  denom <- ((burnin + 1):(n - burnin)) * (1 - ((burnin + 1):(n - burnin)) / n)
  
  unscaled_score_stat <- cumsum(score_est(Y) - mean(score_est(Y)))[(burnin + 1):(n - burnin)]
  return(unscaled_score_stat^2 / denom)
}

# MAIN SIM LOOP
with_progress({
  p <- progressor(steps = nreps)
  
  out_est <- foreach(i = 1:nreps, .combine = rbind,
                     .export = c("get_max", "state_to_compute", "generate_cusum", "generate_oracle",
                                 "generate_oracle_score_est", "generate_scoring",
                                 "spline_score", "cv_spline_score", "oracle_score",
                                 "burnin", "n", "tau", "noise_dist", "sigma", "scale", "df", "delta"),
                     .packages = c("stats")) %dorng% {
    p()  # update progress
                       
    Y <- noise <- sample_from_distribution(n, noise_dist, sd = sigma, scale = scale, df = df)
    Y[(tau+1):n] <- Y[(tau+1):n] + delta
    
    # initialise list to store iteration results
    tau_estimates <- list()
    max_stats <- list()
    tau_hat_cusum <- NULL
    
    # handle dependency of "scoring" on "cusum" tau_hat estimate
    if ("scoring" %in% stats_to_compute){
      cusum_stats <- generate_cusum(Y)
      cusum_results <- get_max(cusum_stats)
      tau_hat_cusum <- cusum_results$tau_hat
      
      # if "cusum" result wanted, store now
      if ("cusum" %in% stats_to_compute){
        tau_estimates$cusum <- tau_hat_cusum
        max_stats$cusum <- cusum_results$max
      }
    }
    
    # compute wanted statistics
    for (stat_name in stats_to_compute){
      if (stat_name == "cusum" && !is.null(tau_hat_cusum)){
        next
      }
      
      current_stats <-switch(stat_name,
                             cusum = generate_cusum(Y),
                             oracle = generate_oracle(Y, oracle_score),
                             oracle_score_est = generate_oracle_score_est(Y, noise),
                             scoring = generate_scoring(Y, tau_est = tau_hat_cusum),
                             stop("Unknown stat: ", stat_name)
                             )
      
      results <- get_max(current_stats)
      tau_estimates[[stat_name]] <- results$tau_hat
      max_stats[[stat_name]] <- results$max
      
    }
    # note: this is messy and needs to be reworked
    names(tau_estimates) <- paste0(names(tau_estimates), ".tau")
    names(max_stats) <- paste0(names(max_stats), ".max")
    
    unlist(c(tau_estimates, max_stats))
    }
})

# View the results
print("Simulation complete. Resulting data frame:")
print(head(out_est))

# for posterity: save file as "dist_n_(sigmaXORscaleORdf)_lambda_delta_sim_results.rds"
saveRDS(out_est, file = "t_1000_0-5_3_0-5_0-3_sim_results.rds")
print("t_1000_0-5_3_0-5_0.3_sim_results.rds")
