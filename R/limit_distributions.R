#' @description
#' Generator for creating score, density, quantile and sampling functions for
#' various limit distributions
#' 
#' Supported distributions:
#' - "eta" (limit dist of the CUSUM estimator, Csorgo and Horvath Thm 1.6.3)
#' - "Gumbel"(TO DO)
#' - "Frechet"(TO DO)
#' - "Weibull" (TO DO)

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
  eta = list(
    params = c(),
    defaults  = list(),
    density = function(params){
      function(x) 1.5 * exp(abs(x)) * (1 - pnorm(1.5 * abs(x)^0.5))
                         - 0.5 * (1 - pnorm(0.5 * abs(x)^0.5))
    }
  )
)

# --- MAIN HANDLER ----
# points to the correct function in .distributions
.create_handler <- function(dist_name, fn_type, ...){
  if (!dist_name %in% names(.distributions)){
    stop("Distribution not supported: ", dist_name, call. = FALSE)
  }
  dist_def <- .distributions[[dist_name]]
  
  user_params <- ...
  params <- .merge_defaults(user_params, dist_def$defaults)
  
  fn_shell <- dist_def[[fn_type]]
  return(fn_shell(params))
}

# ---- PUBLIC FUNCTIONS ----
create_density_function <- function(dist_name, ...){
  .create_handler(dist_name, "density", ...)
}



  