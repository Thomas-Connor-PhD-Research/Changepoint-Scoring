
# 1. SETUP
source("simulations/config_single_change_rate_estimation.R") # <-- Only line to be chang
library(foreach)
library(progressr)
library(future)
library(doFuture)
library(doRNG)
library(rlang)
library(rngtools)

#' Sources all .R files in a specified directory
#'
#' @param path The path to the folder containing your R scripts.
source_dir <- function(path) {
  files <- list.files(path, pattern = "\\.[rR]$", full.names = TRUE)
  for (file in files) {
    tryCatch(
      source(file, local = FALSE), # Source into the worker's global env
      error = function(e) {
        warning(paste("Error sourcing", file, ":", e$message))
      }
    )
  }
}

plan(multisession,
     workers = parallel::detectCores() - 1
)
registerDoFuture()
handlers(global = TRUE)
handlers("progress")

set.seed(simulation_params$seed)
 
with_progress({
  p <- progressor(steps = simulation_params$n_reps+1)

  results_df <- foreach(i = 1:simulation_params$n_reps, .combine = rbind) %dorng% {
    # # Currently each worker must source all functions
    # # TODO: Fix
    source_dir("R/")

    p()
     
    do.call(simulation_params$simulation_function,
            list(data_params = data_params,
                 estimator_params = estimator_params))
    
 }
})


output_data <- list(
  simulation_parameters = simulation_params,
  data_params = data_params,
  estimator_params = estimator_params,
  results = as.data.frame(results_df)
)

plan(sequential)

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

filename <- paste0("results/",
                   simulation_params$simulation_name,
                   "_",
                   timestamp,
                   ".rds")

saveRDS(output_data, file = filename)

print("Simulation complete.")
print(paste("Results saved to:", filename))
print(head(output_data$results))