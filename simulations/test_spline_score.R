plot_single_estimated_versus_true <- function(){
noise_dist_df <- 3
sample_size <- 100
noise_sample <- 2*rt(sample_size, df = 3)
true_score_fn <- create_score_function("t", scale = 2, df = 3)

estimated_score_fn <- spline_score(noise_sample, df = cv_spline_score(noise_sample)$df_min)$rho

# x_grid <- seq(min(noise_sample), max(noise_sample), length.out = 500)
x_grid <- seq(min(noise_sample), max(noise_sample), length.out = 1000)

true_scores <- true_score_fn(x_grid)
estimated_scores <- estimated_score_fn(x_grid)

# Plot the results
plot(x_grid, true_scores, type = "l", col = "blue", lwd = 2,
     main = "Comparison of True vs. Spline-Estimated Score Function",
     xlab = "x", ylab = "rho(x)",
     ylim = range(c(true_scores, estimated_scores)),
     sub = paste("Noise Source: t-distribution with df =", noise_dist_df))

lines(x_grid, estimated_scores, col = "red", lty = 2, lwd = 2)

legend("topright", 
       legend = c("True Score Function", "Spline-Estimated Score"),
       col = c("blue", "red"), 
       lwd = 2, 
       lty = c(1, 2))

grid()
}

.generate_noise <- function(n, dist, 
                           sigma, df, scale ,
                           size_perturb, prop_perturb){
  clean_noise <- sample_from_distribution(n, dist, 
                                          sd = sigma, df = df, scale = scale)
  n_perturb <- floor(n*prop_perturb)
  clean_noise[1:n_perturb] <- clean_noise[1:n_perturb] + size_perturb
  return(clean_noise)
}

compare_distance_against_n <- function(n_values, n_reps = 50,
                                       noise_dist = "normal", 
                                       sigma = 1, df = NULL, scale = NULL,
                                       p_min = 0.01, p_max = 0.99,
                                       prop_perturb = 0, size_perturb = 0){
  library(progressr)

  
  true_score_fn <- create_score_function(noise_dist, sd = sigma, 
                                         df = df, scale = scale)
  true_density_fn <- create_density_function(noise_dist, sd = sigma, 
                                             df = df, scale = scale )
  true_quantile_fn <- create_quantile_function(noise_dist, sd = sigma,
                                               df = df, scale = scale)  
  
  compute_MSE <- function(f_1, f_2, quantile_fn, pmin, pmax, n_pts = 500){
    x_grid <- seq(quantile_fn(pmin), quantile_fn(pmax), length.out = n_pts)
    return(mean((f_1(x_grid) - f_2(x_grid))^2))
  }
  
  compute_density_weighted_MSE <- function(f_1, f_2, density_fn, quantile_fn, pmin, pmax, n_pts = 500){
    p_grid <- seq(pmin, pmax, length.out = n_pts)
    x_grid <- quantile_fn(p_grid)
    squared_errors <- (f_1(x_grid) - f_2(x_grid))^2
    weights <- density_fn(x_grid)
    return(sum(squared_errors * weights)/sum(weights))
  }
  
  results <- data.frame(
    n = n_values,
    mean_MSE = NA,
    se_MSE = NA
  )
  
  with_progress({
    total_steps <- length(n_values) * n_reps
    p <- progressor(steps = total_steps)
    
  for (j in 1:length(n_values)){
    n <- n_values[j]
    MSE <- numeric(n_reps)
    
  for (i in 1:n_reps){
    p(message = sprintf("n = %d | rep = %d", n, i))
    noise_sample <- .generate_noise(n, noise_dist, sigma, df, scale, prop_perturb, size_perturb)
    
    estimated_score_fn <- spline_score(noise_sample, df = cv_spline_score(noise_sample)$df_min)$rho
    
    MSE[i] <- compute_MSE(true_score_fn, estimated_score_fn,
                         true_quantile_fn,
                         p_min, p_max)
    
     # MSE[i] <- compute_density_weighted_MSE(true_score_fn, estimated_score_fn,
     #                      true_density_fn, true_quantile_fn,
     #                      p_min, p_max)
  }
  results$mean_MSE[j] <- mean(MSE)
  results$se_MSE[j] <- sd(MSE) /sqrt(n_reps)
                  
  }
  })
  
  return(results)
}

# With results from above, plot n against mean MSE with SE bars
plot_MSE_with_error <- function(results, log_scale = FALSE){
  # Determine axis type
  log_axis <- if(log_scale) "xy" else ""
  
  # Plot mean M
  plot(
    results$n, results$mean_MSE,
    type = "b",
    pch = 19,
    ylim = c(min(results$mean_MSE - results$se_MSE), max(results$mean_MSE + results$se_MSE)),
    xlab = "Sample size (n)",
    ylab = "Mean Squared Error",
    main = "MSE vs Sample Size - N(0,1)",
    log = log_axis
  )
  
  # Add error bars
  arrows(
    x0 = results$n, y0 = results$mean_MSE - results$se_MSE,
    x1 = results$n, y1 = results$mean_MSE + results$se_MSE,
    angle = 90, code = 3, length = 0.05
  )
}

# With results from above, fit relationsip of the form MSE = a*n^-b
MSE_regression <- function(results){
  fit <- lm(log(results$mean_MSE) ~ log(results$n))
  summary(fit)
  
  b <- -coef(fit)[2]
  a <- exp(coef(fit)[1])
  cat("Estimated MSE(n) ≈", a, "* n^(-", b, ")\n")
  
  plot(results$n, results$mean_MSE, log = "xy", pch = 19, xlab = "n", ylab = "Mean MSE")
  lines(results$n, a * results$n^(-b), col = "red", lwd = 2)
}

# Perturb error dist. to mimic incorrect initial CP / jump-point to test
# scoring estimation robustness, calculate MSE

test_score_estimation_robustness <- function(){
n_values <- c(200, 500, 1000)
prop <- 0.5
size_jumps <- c(0, 0.01, 0.1, 1)

all_results <- list()  

for (i in 1:length(size_jumps)){
  cat("Running size_perturb =", size_jumps[i], "\n")
  df <- compare_distance_against_n(n_values, n_reps = 100,
                                         noise_dist = "normal", 
                                         sigma = 1, prop_perturb = prop, 
                             size_perturb = size_jumps[i]
                             )
  
  all_results[[paste0("size_", size_jumps[i])]] <- df
}

return(all_results)
  
}

plot_MSE_with_error_list <- function(results_list, log_scale = FALSE, colors = NULL) {
  log_axis <- if (log_scale) "xy" else ""
  
  # Determine global ylim
  all_mean <- unlist(lapply(results_list, function(df) df$mean_MSE))
  all_se   <- unlist(lapply(results_list, function(df) df$se_MSE))
  ylim <- c(min(all_mean - all_se), max(all_mean + all_se))
  
  # Set colors
  if (is.null(colors)) colors <- rainbow(length(results_list))
  
  # Plot the first data frame to initialize the plot
  first_df <- results_list[[1]]
  plot(
    first_df$n, first_df$mean_MSE,
    type = "b",
    pch = 19,
    col = colors[1],
    ylim = ylim,
    xlab = "Sample size (n)",
    ylab = "Mean Squared Error",
    main = "MSE vs Sample Size for Different Size Jumps",
    log = log_axis
  )
  
  # Add error bars for the first data frame
  arrows(
    x0 = first_df$n, y0 = first_df$mean_MSE - first_df$se_MSE,
    x1 = first_df$n, y1 = first_df$mean_MSE + first_df$se_MSE,
    angle = 90, code = 3, length = 0.05,
    col = colors[1]
  )
  
  # Add the rest of the data frames
  if (length(results_list) > 1) {
    for (i in 2:length(results_list)) {
      df <- results_list[[i]]
      points(df$n, df$mean_MSE, type = "b", pch = 19, col = colors[i])
      arrows(
        x0 = df$n, y0 = df$mean_MSE - df$se_MSE,
        x1 = df$n, y1 = df$mean_MSE + df$se_MSE,
        angle = 90, code = 3, length = 0.05,
        col = colors[i]
      )
    }
  }
  
  # Add a legend
  legend(
    "topright",
    legend = names(results_list),
    col = colors,
    lty = 1,
    pch = 19
  )
}

n_values <- c(200, 500, 1000)
normal_results <- compare_distance_against_n(n_values)

plot_MSE_with_error(normal_results)
