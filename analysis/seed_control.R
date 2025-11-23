#' Get the specific parallel stream seed used by %dorng%
#'
#' This function reproduces the L'Ecuyer-CMRG stream generation
#' used by doRNG to find the seed for a specific iteration. Dependencies: 
#' rngtools.
#'
#' @param master_seed Numeric. The single seed (e.g., 123) provided to the
#'   entire %dorng% loop (via set.seed()).
#' @param iteration_no Numeric. The specific parallel worker iteration 
#'   (e.g., 5) you want to investigate. Must be >= 1.
#'
#' @return An RNGseed object. Pass this to setRNG()
#'   to activate the stream.
#'
get_dorng_stream <- function(master_seed, iteration_no) {
  
  # Ensure dependencies
  if (!requireNamespace("rngtools", quietly = TRUE)) {
    stop("Please install the 'rngtools' package to use this function.")
  }
  
  # Input validation
  if (iteration_no < 1 || !is.numeric(iteration_no) || (iteration_no %% 1 != 0)) {
    stop("'iteration_no' must be a positive integer (>= 1).")
  }
  
  # Set the master seed using the same method as doRNG
  set.seed(master_seed, kind = "L'Ecuyer-CMRG")
  
  current_stream_seed <- .Random.seed
  
  if (iteration_no > 1) {
    # We need to advance (iteration_no - 1) times
    for (i in 1:(iteration_no - 1)) {
      current_stream_seed <- parallel::nextRNGStream(current_stream_seed)
    }
  }
  
  # Return the final seed vector
  return(current_stream_seed)
}