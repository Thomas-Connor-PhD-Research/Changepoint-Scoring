#' @description
#' Returns useful functions (score, density, quantile) for common error 
#' distributions (with variance scaling)
#'
#' Supported distributions:
#' - "normal"
#' - "t"
#' - "cauchy"
#' - "logistic" (with location 0)
#' - "laplace"
#' - "exponential" (TBI)
#' 
#' @param distribution_name Character; name of the distribution.
#' @param sd Numeric; scale or standard deviation (default = 1).
#' @param ... Additional parameters (e.g., `df` for t-distribution, `scale`).
#'
#' @return A function of `x` (for score and density) or `p` (for quantile).
#'
#' @examples
#' score_fn <- create_score_function("t", sd = 2, df = 3)
#' density_fn <- create_density_function("t", sd = 2, df = 3)
#' quantile_fn <- create_quantile_function("t", sd = 2, df = 3)
#' 
#' x <- seq(-10, 10, length.out = 1000)
#' plot(x, score_fn(x), type = "l")

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- SCORE FUNCTION ----
create_score_function <- function(distribution_name, ...) {
  params <- list(...)
  
  if (distribution_name == "normal") {
    sd <- params$sd
    return(function(x) -x / (sd^2))
  }
  
  if (distribution_name == "t") {
    df <- params$df
    scale <- params$scale
    if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
    if (is.null(scale)) stop("Scale 'scale' must be provided for the t-distribution.")
    return(function(x) -(df + 1) * x / (df*scale^2 + x^2))
  }
  
  if (distribution_name == "laplace") {
    return(function(x) -sign(x) / scale)
  }
  
  if (distribution_name == "cauchy") {
    scale <- params$scale %||% 1
    return(function(x) -2 * x / (scale^2 + x^2))
  }
  
  if (distribution_name == "logistic") {
    scale <- params$scale %||% 1
    return(function(x) -(1 / (scale)) * tanh(x / (2 * scale)))
  }
  
  stop("Distribution not supported: ", distribution_name)
}

# ---- DENSITY FUNCTION ----
create_density_function <- function(distribution_name, ...) {
  params <- list(...)
  
  if (distribution_name == "normal") {
    sd <- params$sd
    return(function(x) dnorm(x / sd) / sd)
  }
  
  if (distribution_name == "t") {
    df <- params$df
    scale <- params$scale %||% 1
    if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
    return(function(x) dt(x / scale, df = df) / scale)
  }
  
  if (distribution_name == "laplace") {
    return(function(x) (1 / (2 * scale)) * exp(-abs(x) / s))
  }
  
  if (distribution_name == "cauchy") {
    scale <- params$scale %||% 1
    return(function(x) dcauchy(x / scale) / scale)
  }
  
  if (distribution_name == "logistic") {
    scale <- params$scale %||% 1
    return(function(x) dlogis(x / scale) /  scale)
  }
  
  stop("Distribution not supported: ", distribution_name)
}


# ---- QUANTILE FUNCTION ----
create_quantile_function <- function(distribution_name, ...) {
  params <- list(...)
  
  if (distribution_name == "normal") {
    sd <- params$sd %||% 1
    return(function(p) sd * qnorm(p))
  }
  
  if (distribution_name == "t") {
    df <- params$df
    scale <- params$scale %||% 1
    if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
    return(function(p) scale * qt(p, df = df))
  }
  
  if (distribution_name == "laplace") {
    return(function(p) {
      ifelse(p < 0.5,
             sd * log(2 * p),
             -sd * log(2 * (1 - p)))
    })
  }
  
  if (distribution_name == "cauchy") {
    scale <- params$scale %||% 1
    return(function(p)  scale * qcauchy(p))
  }
  
  if (distribution_name == "logistic") {
    scale <- params$scale %||% 1
    return(function(p) scale * qlogis(p))
  }
  
  stop("Distribution not supported: ", distribution_name)
}

# ---- SAMPLE NOISE ----
sample_from_distribution <- function(n, distribution_name, ...) {
  params <- list(...)
  
  if (distribution_name == "normal") {
    sd <- params$sd %||% 1
    return(sd * rnorm(n))
  }
  
  if (distribution_name == "t") {
    scale <- params$scale %||% 1
    df <- params$df
    if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
    return(scale * rt(n, df = df))
  }
  
  if (distribution_name == "laplace") {
    # via inverse CDF method
    u <- runif(n)
    return(scale * ifelse(u < 0.5, log(2 * u), -log(2 * (1 - u))))
  }
  
  if (distribution_name == "cauchy") {
    scale <- params$scale %||% 1
    return(scale * rcauchy(n))
  }
  
  if (distribution_name == "logistic") {
    scale <- params$scale %||% 1
    return(scale * rlogis(n))
  }
  
  stop("Distribution not supported: ", distribution_name)
}