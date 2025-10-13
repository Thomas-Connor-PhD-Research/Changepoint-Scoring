#' @description
#' Generator for creating score, density, quantile and sampling functions for
#' various noise distributions with mean 0 
#' 
#' Supported distributions:
#' - "normal" (params: sd)
#' - "t" (params: df, scale)
#' - "cauchy" (params: scale)
#' - "logistic" (params: scale)
#' - "laplace" (params: scale)

# ---- HELPER FUNCTIONS ----
# Merge user-provided params with a list of defaults
.merge_defaults <- function(user_params, default_params){ 
  final_params <- default_params
  
  for (name in names(user_params)){
    if (name %in% names(final_params)){
      final_params[[name]] <- user_params[[name]]
    }
  }
  return(final_params)
}

# ---- DISTRIBUTION DEFINITIONS ----
# Contains all logic for each distribution
.distributions <- list(
  normal = list(
    params = c("sd"),
    defaults  = list(sd = 1),
    score = function(params){
      function(x) -x / (params$sd^2)
    },
    density = function(params){
      function(x) dnorm(x, mean = 0, sd = params$sd)
    },
    quantile = function(params){
      function(p) qnorm(p, mean = 0, sd = params$sd)
    },
    sampler = function(params){
      function(n) rnorm(n, mean = 0, params$sd)
    }
  ),
  t = list(
    params = c("df","scale"),
    defaults  = list(df = 1, scale = 1),
    score = function(params){
      function(x) -(params$df + 1) * x / (params$df*params$scale^2 + x^2)
    },
    density = function(params){
      function(x) 1/params$scale*dt(x/params$scale, df = params$df)
    },
    quantile = function(params){
      function(p) params$scale*qt(p, df = params$df)
    },
    sampler = function(params){
      function(n) params$scale*rt(n, df = params$df)
    }
  ),
  laplace = list(
    params = c("scale"),
    defaults = list(scale = 1),
    score = function(params) {
      function(x) -sign(x) / params$scale
    },
    density = function(params) {
      function(x) (1 / (2 * params$scale)) * exp(-abs(x) / params$scale)
    },
    quantile = function(params) {
      function(p) ifelse(p < 0.5, params$scale * log(2 * p), -params$scale * log(2 * (1 - p)))
    },
    sampler = function(params) {
      function(n) {
        u <- runif(n)
        ifelse(u < 0.5, params$scale * log(2 * u), -params$scale * log(2 * (1 - u)))
      }
    }
  ),
  cauchy = list(
    params = c("scale"),
    defaults = list(scale = 1),
    score = function(params) {
      function(x) -2 * x / (params$scale^2 + x^2)
    },
    density = function(params) {
      function(x) dcauchy(x, location = 0, scale = params$scale)
    },
    quantile = function(params) {
      function(p) qcauchy(p, location = 0, scale = params$scale)
    },
    sampler = function(params) {
      function(n) rcauchy(n, location = 0, scale = params$scale)
    }
  ),
  logistic = list(
    params = c("scale"),
    defaults = list(scale = 1),
    score = function(params) {
      function(x) -(1 / params$scale) * tanh(x / (2 * params$scale))
    },
    density = function(params) {
      function(x) dlogis(x, location = 0, scale = params$scale)
    },
    quantile = function(params) {
      function(p) qlogis(p, location = 0, scale = params$scale)
    },
    sampler = function(params) {
      function(n) rlogis(n, location = 0, scale = params$scale)
    }
  )
)

# --- MAIN HANDLER ----
# points to the correct function in .distributions
.create_handler <- function(dist_name, fn_type, user_params){
  if (!dist_name %in% names(.distributions)){
    stop("Distribution not supported: ", dist_name, call. = FALSE)
  }
  dist_def <- .distributions[[dist_name]]

  params <- .merge_defaults(user_params, dist_def$defaults)
  
  fn_shell <- dist_def[[fn_type]]
  return(fn_shell(params))
}

# ---- PUBLIC FUNCTIONS ----
create_score_function <- function(dist_name, params){
  .create_handler(dist_name, "score", params)
}

create_density_function <- function(dist_name, params){
  .create_handler(dist_name, "density", params)
}

create_quantile_function <- function(dist_name, params){
  .create_handler(dist_name, "quantile", params)
}

sample_from_distribution <- function(n, dist_name, params){
  sampler_shell <- .create_handler(dist_name, "sampler", params)
  sampler_shell(n)
}


# 
# 
# 
# # ---- SCORE FUNCTION ----
# create_score_function <- function(distribution_name, ...) {
#   params <- list(...)
#   
#   if (distribution_name == "normal") {
#     sd <- params$sd
#     return(function(x) -x / (sd^2))
#   }
#   
#   if (distribution_name == "t") {
#     df <- params$df
#     scale <- params$scale
#     if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
#     if (is.null(scale)) stop("Scale 'scale' must be provided for the t-distribution.")
#     return(function(x) -(df + 1) * x / (df*scale^2 + x^2))
#   }
#   
#   if (distribution_name == "laplace") {
#     return(function(x) -sign(x) / scale)
#   }
#   
#   if (distribution_name == "cauchy") {
#     scale <- params$scale %||% 1
#     return(function(x) -2 * x / (scale^2 + x^2))
#   }
#   
#   if (distribution_name == "logistic") {
#     scale <- params$scale %||% 1
#     return(function(x) -(1 / (scale)) * tanh(x / (2 * scale)))
#   }
#   
#   stop("Distribution not supported: ", distribution_name)
# }
# 
# # ---- DENSITY FUNCTION ----
# create_density_function <- function(distribution_name, ...) {
#   params <- list(...)
#   
#   if (distribution_name == "normal") {
#     sd <- params$sd
#     return(function(x) dnorm(x / sd) / sd)
#   }
#   
#   if (distribution_name == "t") {
#     df <- params$df
#     scale <- params$scale %||% 1
#     if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
#     return(function(x) dt(x / scale, df = df) / scale)
#   }
#   
#   if (distribution_name == "laplace") {
#     return(function(x) (1 / (2 * scale)) * exp(-abs(x) / s))
#   }
#   
#   if (distribution_name == "cauchy") {
#     scale <- params$scale %||% 1
#     return(function(x) dcauchy(x / scale) / scale)
#   }
#   
#   if (distribution_name == "logistic") {
#     scale <- params$scale %||% 1
#     return(function(x) dlogis(x / scale) /  scale)
#   }
#   
#   stop("Distribution not supported: ", distribution_name)
# }
# 
# 
# # ---- QUANTILE FUNCTION ----
# create_quantile_function <- function(distribution_name, ...) {
#   params <- list(...)
#   
#   if (distribution_name == "normal") {
#     sd <- params$sd %||% 1
#     return(function(p) sd * qnorm(p))
#   }
#   
#   if (distribution_name == "t") {
#     df <- params$df
#     scale <- params$scale %||% 1
#     if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
#     return(function(p) scale * qt(p, df = df))
#   }
#   
#   if (distribution_name == "laplace") {
#     return(function(p) {
#       ifelse(p < 0.5,
#              sd * log(2 * p),
#              -sd * log(2 * (1 - p)))
#     })
#   }
#   
#   if (distribution_name == "cauchy") {
#     scale <- params$scale %||% 1
#     return(function(p)  scale * qcauchy(p))
#   }
#   
#   if (distribution_name == "logistic") {
#     scale <- params$scale %||% 1
#     return(function(p) scale * qlogis(p))
#   }
#   
#   stop("Distribution not supported: ", distribution_name)
# }
# 
# # ---- SAMPLE NOISE ----
# sample_from_distribution <- function(n, distribution_name, ...) {
#   params <- list(...)
#   
#   if (distribution_name == "normal") {
#     sd <- params$sd %||% 1
#     return(sd * rnorm(n))
#   }
#   
#   if (distribution_name == "t") {
#     scale <- params$scale %||% 1
#     df <- params$df
#     if (is.null(df)) stop("Degrees of freedom 'df' must be provided for the t-distribution.")
#     return(scale * rt(n, df = df))
#   }
#   
#   if (distribution_name == "laplace") {
#     # via inverse CDF method
#     u <- runif(n)
#     return(scale * ifelse(u < 0.5, log(2 * u), -log(2 * (1 - u))))
#   }
#   
#   if (distribution_name == "cauchy") {
#     scale <- params$scale %||% 1
#     return(scale * rcauchy(n))
#   }
#   
#   if (distribution_name == "logistic") {
#     scale <- params$scale %||% 1
#     return(scale * rlogis(n))
#   }
#   
#   stop("Distribution not supported: ", distribution_name)
# }