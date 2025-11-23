#' Create a "trimmed" version of a function
#'
#' This is a function factory. It takes a function 'f' and a
#' threshold 'C' and returns a *new* function. This new function
#' will evaluate 'f' on inputs that have been "clipped" or "trimmed"
#' to the range [-C, C].
#'
#' @param f The original continuous function (e.g., function(x) -x).
#' @param C The (positive) threshold value for |x|.
#'
#' @return A new function, f_trimmed(x), which evaluates f
#'         on inputs clipped to [-C, C].
trim_function <- function(f, C) {
  
  # Ensure the threshold is positive
  threshold <- abs(C)
  
  # Define and return the new "trimmed" function.
  # This inner function 'closes over' f and threshold.
  f_trimmed <- function(x) {
    
    # Clip the input vector x to the range [-threshold, threshold]
    # pmax ensures no value is lower than -threshold
    # pmin ensures no value is higher than +threshold
    x_clipped <- pmin(pmax(x, -threshold), threshold)
    
    # Evaluate the original function f with the clipped inputs
    return(f(x_clipped))
  }
  
  # Return the function itself, not its evaluated result
  return(f_trimmed)
}
