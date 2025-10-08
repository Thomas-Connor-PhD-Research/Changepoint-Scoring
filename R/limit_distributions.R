create_density_function <- function(distribution_name, ...){
  params <- list(...)
  
  if (distribution_name == "eta"){
    sd <- params$sd
    return(function(x) 1.5 * exp(abs(x)) * (1 - pnorm(1.5 * abs(x)^0.5))
                                 - 0.5 * (1 - pnorm(0.5 * abs(x)^0.5)))
  }
}

create_distribution_function <- function(distribution_name, ...){
  params <- list(...)
  
  if (distribution_name == "eta"){
    stop("Error: Not Implemented")
  }
}
  
  