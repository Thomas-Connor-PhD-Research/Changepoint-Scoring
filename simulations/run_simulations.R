
# 1. SETUP
source("simulations/config_single_change.R") # <-- Only line to be changed
# source("simulations/config_test_spline_score.R")


library(foreach)
library(progressr)
library(future)
library(doFuture)
library(doRNG)
library(rlang)
library(rngtools)

plan(multisession,
     workers = parallel::detectCores() - 1
)

 registerDoFuture()
 handlers(global = FALSE)

with_progress({
  p <- progressor(steps = sim_params$n_reps+1)

  results_df <- foreach(i = 1:sim_params$n_reps, .combine = rbind) %dorng% {
    set.seed(i) # for reproducability
    
    # Currently each worker must source all functions
    # TODO: Fix
    source("R/simulation_functions.R")
    source("R/noise_distributions.R")
    source("R/score_estimation2.R")
    source("R/scoring_spline_analysis_functions.R")
    
    p()
    
    do.call(sim_params$simulation_function, list(params = sim_params))
 }
})


output_data <- list(
  parameters = sim_params,
  results = as.data.frame(results_df)
)

plan(sequential)

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

filename <- paste0("results/",
                   sim_params$simulation_name,
                   "_",
                   timestamp,
                   ".rds")

saveRDS(output_data, file = filename)

print("Simulation complete.")
print(paste("Results saved to:", filename))
print(head(output_data$results))