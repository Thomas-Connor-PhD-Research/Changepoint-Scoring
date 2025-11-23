
# 1. SETUP
source("simulations/config_single_change.R") # <-- Only line to be changed



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

set.seed(10)
results_df <- foreach(i = 1:10, .combine = rbind) %dorng% {
  sample_from_distribution(1, "normal", list(sd=1))
}

results = as.data.frame(results_df)
plan(sequential)
print(results)

set.seed(10)

for (j in 1:10){
  .Random.seed <- get_dorng_stream(10, j)
  print(sample_from_distribution(1, "normal", list(sd=1)))
}
