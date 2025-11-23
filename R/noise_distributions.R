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
    score_deriv = function(params){
      function(x) - 1/params$sd^2
    },
    density = function(params){
      function(x) dnorm(x, mean = 0, sd = params$sd)
    },
    quantile = function(params){
      function(p) qnorm(p, mean = 0, sd = params$sd)
    },
    sampler = function(params){
      function(n) rnorm(n, mean = 0, params$sd)
    },
    stein_fisher_info = function(params){ 
      1/(params$sd^2)
      },
    variance = function(params){
      params$sd^2
    }
  ),
  t = list(
    params = c("df","scale"),
    defaults  = list(df = 1, scale = 1),
    score = function(params){
      function(x) -(params$df + 1) * x / (params$df*params$scale^2 + x^2)
    },
    score_deriv = function(params){
      function(x) (params$df + 1) * (x^2 - params$df*params$scale^2) / (params$df*params$scale^2 + x^2)^2
    },
    density = function(params){
      function(x) 1/params$scale*dt(x/params$scale, df = params$df)
    },
    quantile = function(params){
      function(p) params$scale*qt(p, df = params$df)
    },
    sampler = function(params){
      function(n) params$scale*rt(n, df = params$df)
    },
    stein_fisher_info = function(params){
      (params$df + 1)/(params$scale^2*(params$df + 3))
    }
  ),
  laplace = list(
    params = c("scale"),
    defaults = list(scale = 1),
    score = function(params) {
      function(x) -sign(x) / params$scale
    },
    score_deriv = function(params){
      NULL
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
    },
    stein_fisher_info = function(params){
      (1 / params$scale^2)
    }
  ),
  cauchy = list(
    params = c("scale"),
    defaults = list(scale = 1),
    score = function(params) {
      function(x) -2 * x / (params$scale^2 + x^2)
    },
    score_deriv = function(params){
      2 * (x^2 - params$scale^2) / (params$scale^2 + x^2)^2
    },
    density = function(params) {
      function(x) dcauchy(x, location = 0, scale = params$scale)
    },
    quantile = function(params) {
      function(p) qcauchy(p, location = 0, scale = params$scale)
    },
    sampler = function(params) {
      function(n) rcauchy(n, location = 0, scale = params$scale)
    },
    stein_fisher_info = function(params){
      1 / (2 * params$scale^2)
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
    },
    stein_fisher_info = function(params){
      1 / (3*params$scale^2)
    }
  ),
  bimodal = list(
    params = c("peak", "scale"),
    defaults = list(peak = 2, scale = 1),
    score = function(params) {
      # g(x) = -f'(x)/f(x)
      # f(x) = 0.5 * p1 + 0.5 * p2
      # f'(x) = 0.5 * p1 * -(x+peak)/s^2 + 0.5 * p2 * -(x-peak)/s^2
      # g(x) = [ (x+peak)*p1 + (x-peak)*p2 ] / [ s^2 * (p1 + p2) ]
      function(x) {
        s2 <- params$scale^2
        p <- params$peak
        p1 <- dnorm(x, -p, params$scale)
        p2 <- dnorm(x, p, params$scale)
        # Add a small epsilon to denominator to avoid 0/0 if p1+p2 is 0
        ((x + p) * p1 + (x - p) * p2) / (s2 * (p1 + p2) + 1e-100)
      }
    },
    density = function(params) {
      function(x) {
        0.5 * dnorm(x, mean = -params$peak, sd = params$scale) +
          0.5 * dnorm(x, mean = params$peak, sd = params$scale)
      }
    },
    quantile = function(params) {
      # No closed-form solution for mixture quantile
      function(p) stop("Quantile function for bimodal mixture is not implemented.", call. = FALSE)
    },
    sampler = function(params) {
      function(n) {
        # Sample n1 from first mode, n2 from second
        n1 <- rbinom(1, n, 0.5)
        n2 <- n - n1
        c(rnorm(n1, -params$peak, params$scale),
          rnorm(n2, params$peak, params$scale))
      }
    },
    stein_fisher_info = function(params){
      # Extremely complex integral
      NA
    }
  
),
lognormal= list(
  params = c("meanlog", "sdlog"),
  defaults = list(meanlog = 0, sdlog = 1),
  score = function(params) {
    # f(x) = dlnorm(x + mu, m, s)
    # g(x) = -f'(x)/f(x)
    # ln(f(x)) = -ln(x+mu) - ... - (ln(x+mu) - m)^2 / (2s^2)
    # d/dx ln(f(x)) = -1/(x+mu) - 2(ln(x+mu)-m)/(2s^2) * 1/(x+mu)
    # g(x) = (1 / (x+mu)) * (1 + (ln(x+mu) - m) / s^2)
    function(x) {
      m <- params$meanlog
      s <- params$sdlog
      s2 <- s^2
      mu <- exp(m + s2 / 2) # The mean, used for centering
      x_plus_mu <- x + mu
      
      # Score is undefined for x <= -mu
      ifelse(x_plus_mu <= 0, NA,
             (1 / x_plus_mu) * (1 + (log(x_plus_mu) - m) / s2)
      )
    }
  },
  density = function(params) {
    m <- params$meanlog
    s <- params$sdlog
    mu <- exp(m + s^2 / 2)
    function(x) dlnorm(x + mu, m, s)
  },
  quantile = function(params) {
    m <- params$meanlog
    s <- params$sdlog
    mu <- exp(m + s^2 / 2)
    function(p) qlnorm(p, m, s) - mu
  },
  sampler = function(params) {
    m <- params$meanlog
    s <- params$sdlog
    mu <- exp(m + s^2 / 2)
    function(n) rlnorm(n, m, s) - mu
  },
  stein_fisher_info = function(params){
    # This is not a simple location family, Fisher info is complex
    NA
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

get_stein_fisher_info <- function(dist_name, params){
  .create_handler(dist_name, "stein_fisher_info", params)
}


get_variance <- function(dist_name, params){
  .create_handler(dist_name, "variance", params)
}


# -- EXAMPLES --
# dist_name <- "lognormal"
# params <- c(meanlog = 1, sdlog=1)
# density_fn <- create_density_function(dist_name, params)