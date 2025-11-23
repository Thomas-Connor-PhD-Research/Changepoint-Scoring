library(Iso)

spline_score_estimation <- function(x, df=5, tol=1e-3, nmax=NULL){
  #' Univariate score estimation via the smoothing spline method of Cox 1985
  #'  and Ng 1994.
  #'
  #' @param x vector of datapoints
  #' @param df vector of smoothing parameters for the
  #'     non-parametric score estimator, corresponding to the
  #'     effective degrees of freedom for a smoothing spline.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return score function "rho" and derivative "drho", which take vector
  #' input and yield a vector of score estimates corresponding to each df (in a
  #' list if there are multiple df values). Also output the vector "df".
  #' @export
  #'
  #' @examples
  #' # Single bandwidth
  #' x <- stats::rlogis(100)
  #' spl <- spline_score(x, df=6)
  #' spl$rho(x)
  #' spl$drho(x)
  #'
  #' # Multiple bandwidths simultaneously
  #' x <- stats::rt(n=100, df=4)
  #' spl <- spline_score(x, df=c(2,5,10))
  #' spl$rho(x)
  
  
  
  if (!requireNamespace("graphics", quietly = TRUE)) {
    warning("Package \"graphics\" required for sort_bin. Matrices may be
            singular.")
    w <- rep(1,length(x))
    x_sort <- x
  } else {
    bin <- sort_bin(x=x, tol=tol, nmax=nmax)
    x_sort <- bin$x_sort
    w <- bin$w
    
    if (length(x_sort) < 4){
      warning("smooth.spline requires at least 4 unique x values.
                Binning would violate this, so not binning.")
      w <- rep(1,length(x))
      x_sort <- x
    }
  }
  
  pseudo_y <- ng_pseudo_response(x=x_sort, w=w)
  
  rho <- function(x){
    out <- lapply(df, function(dfval){
      sm <- stats::smooth.spline(x=x_sort, y=pseudo_y, w=w, df=dfval)
      return(as.vector(stats::predict(sm, x)$y))
    })
    # if only one df specified, unlist
    if (length(out)==1){
      out <- out[[1]]
    }
    return(out)
  }
  
  drho <- function(x){
    out <- lapply(df, function(dfval){
      sm <- stats::smooth.spline(x=x_sort, y=pseudo_y, w=w, df=dfval)
      return(as.vector(stats::predict(sm, x, deriv=1)$y))
    })
    # if only one df specified, unlist
    if (length(out)==1){
      out <- out[[1]]
    }
    return(out)
  }
  
  
  return(list("rho" = rho,
              "drho" = drho,
              "df" = df))
}

cv_spline_score <- function(x, df=2:15, nfolds=5L, tol=1e-3, nmax=NULL, foldid=NULL){
  #' K-fold cross-validation for spline_score.
  #'
  #' @param x vector of datapoints
  #' @param df vector of smoothing parameters for the
  #'     non-parametric score estimator, corresponding to the
  #'     effective degrees of freedom for a smoothing spline.
  #' @param nfolds integer number of cross-validation folds.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return list of 5 elements: df vector, cv vector of corresponding
  #'     cross-validation scores, se vector of standard error estimates,
  #'     df_min cross-validation minimiser, df_1se largest smoothing
  #'     parameter within CV score within one standard error of df_min.
  #'
  #' @export
  #' @examples
  #' set.seed(0)
  #' x <- stats::rt(100, df=4)
  #' cv_spline_score(x)
  #'
  #' x <- stats::rlogis(500)
  #' cvspl <- cv_spline_score(x)
  #' cvspl$df_min
  
  n <- length(x)
  ndf <- length(df)
  
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length.out=n))
  } else {
    nfolds <- max(foldid)
  }
  cv_folds <- matrix(NA,nrow=nfolds,ncol=ndf)
  

  
  for (fold in seq(nfolds)){
    
    which <- (foldid == fold)
    spline <- spline_score_estimation(x=x[!which], df=df, tol=tol, nmax=nmax)
    
    # Evaluate score estimate on holdout fold
    score <- spline$rho(x[which])
    dscore <- spline$drho(x[which])
    # Compute CV scores for each df
    cv_folds[fold, ] <- mapply(function(x,y){mean(x^2+2*y)},
                               x = score,
                               y = dscore)
  }
  
  cv <- colMeans(cv_folds)
  se <- apply(cv_folds,2,stats::sd)/sqrt(nfolds)
  
  dfmin_index <- which.min(cv)
  
  cv_target_1se <- cv[dfmin_index]+se[dfmin_index]
  df1se_index <- which(df==min(df[cv < cv_target_1se]))
  
  cv_target_2se <- cv[dfmin_index]+2*se[dfmin_index]
  df2se_index <- which(df==min(df[cv < cv_target_2se]))
  
  cv_target_3se <- cv[dfmin_index]+3*se[dfmin_index]
  df3se_index <- which(df==min(df[cv < cv_target_3se]))
  
  
  df_min <- df[dfmin_index]
  df_1se <- df[df1se_index]
  df_2se <- df[df2se_index]
  df_3se <- df[df3se_index]
  
  return(list('df'=df,
              'cv'=cv,
              'se'=se,
              'df_min'=df_min,
              'df_1se'=df_1se,
              'df_2se'=df_2se,
              'df_3se'=df_3se))
  
}

ng_pseudo_response <- function(x, w=rep(1,length(x))){
  #' Generate pseudo responses as in Ng 1994 to enable univariate score
  #' estimation by standard smoothing spline regression.
  #'
  #' Pseudo responses should be regarded as a computational tool, not
  #' as an estimate of the score itself.
  #'
  #' @param x vector of covariates.
  #' @param w vector of weights.
  #'
  #' @return A vector of score estimates.
  #' @export
  #'
  #' @examples
  #' x <- seq(-3,3, length.out=50)
  #' ng_pseudo_response(x)
  
  order <- order(x)
  n <- length(x)
  # Sort x and w
  x <- x[order]
  w <- w[order]
  # Compute interval widths
  h  <- x[-1] - x[-n]
  
  # Compute functions of interval widths
  wih <- c(w[1:(n-2)]/h[1:(n-2)], (w[n-1]+w[n])/h[n-1])
  wh <- c(w[1:(n-2)]*h[1:(n-2)], (w[n-1]-w[n]/2)*h[n-1])
  
  
  # Notation as in Ng [1994] and Ng [2003].
  a_vec <- c(wih,0)-c(0,wih) # -A^T P 1
  c_vec <- (wh[-(n-1)] + 2*wh[-1])/3 # C^T P 1
  
  
  ih <- 1/h
  
  R <- diag(2 * (h[-(n-1)] + h[-1]) / 3 , nrow=n-2, ncol=n-2)
  R[row(R) - col(R) == 1] <- h[-c(1,n-1)] / 3
  R[row(R) - col(R) == -1] <- h[-c(1,n-1)] / 3
  
  # Sparsify matrices
  if (requireNamespace("Matrix", quietly = TRUE)) {
    R <- Matrix::Matrix(R, sparse=TRUE)
  }
  
  Q <- diag(ih[-(n-1)], nrow=n, ncol=n-2)
  Q[(row(Q) - col(Q)) == 1] <- -(ih[-(n-1)] + ih[-1])
  Q[(row(Q) - col(Q)) == 2] <- ih[-1]
  
  
  # Sparsify matrices
  if (requireNamespace("Matrix", quietly = TRUE)) {
    Q <- Matrix::Matrix(Q, sparse=TRUE)
  }
  
  
  z <- as.vector(solve(R,c_vec))
  y <- (1/w) * as.vector(a_vec + Q%*%z)
  
  # Fix ordering of output.
  y[order] <- y
  
  return(y)
  
}

sort_bin <- function(x, tol=1e-5, nmax=NULL){
  #' Sort and bin x within a specified tolerance, using hist().
  #'
  #' @param x vector of covariates.
  #' @param tol numeric tolerance, minimum distance between neighbouring points,
  #'     to avoid singularities.
  #' @param nmax if specified, overrides tol as maximal number of unique points.
  #'
  #' @return list with three elements. x_sort is sorted and binned x,
  #'      w is a vector of weights corresponding to the frequency of each bin,
  #'      order is a vector specifying the ordering of x into the binned values
  #'      sort_x.
  
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  if (is.null(nmax)){
    br <- ceiling((max(x)-min(x))/tol) # number of bins
  } else{
    br <- nmax
  }
  
  br <- min(br, 1e+6) # hist caps br to 1e+6 anyway, so this just avoid a warning message
  
  hist <- graphics::hist(x, br, right=FALSE, plot=FALSE) # assign elements of x to bins
  counts <- hist$counts
  mids <- hist$mids
  breaks <- hist$breaks
  
  # remove empty bins
  w <- counts[counts>0] # frequencies in non-empty bins
  x_sort <- mids[counts>0] # midpoints of non-empty bins
  pos_breaks <- breaks[c(1,1+which(counts>0))] # breakpoints for non-empty bins
  
  order <- findInterval(x, pos_breaks)
  
  return(list("x_sort"=x_sort, "w"=w, "order"=order))
  
}

isotonize_score_given_density <- function(grid, densities, J_eval = NULL, k = 3000,
                                          k1=1, k2=k, ptm) {
  K <- length(grid)
  
  # Interpolate between density values on the grid
  f <- approxfun(grid, densities, method = "linear", rule = 2)
  masses <- (densities[-1] + densities[-K])/2 * diff(grid)
  Fx <- c(0, cumsum(masses))
  masses <- masses / Fx[K]
  Fx <- Fx / Fx[K]
  
  Q <- approxfun(Fx, grid, method = "linear", rule = 2, ties = max)
  
  u <- 0:k / k # customise this for light-tailed densities for plotting purposes?
  J <- f(Q(u))
  
  # Evaluate J at a given point J_eval (used to estimate the intercept standard
  # error when intercept.selection = "quantile" in asm.fit)
  if (!is.null(J_eval)) J_eval <- f(Q(J_eval))
  
  init_scores <- diff(J) * k # psihat circ Qhat(0, 1/k, 2/k, ..., 1)
  
  ## init_scores_sub = init_scores[k1:k2]
  domain = (Q(u[-1]) + Q(u[-(k + 1)])) / 2
  
  fn_vals = pava(init_scores, decreasing = TRUE,
                 long.out = FALSE, stepfun = FALSE)
  
  fn_vals = fn_vals[k1:k2]
  domain = domain[k1:k2]
  
  jumps = diff(fn_vals) != 0
  jumps_lshift = c(jumps, TRUE)
  jumps_rshift = c(TRUE, jumps)
  non_redundant_index <- jumps_lshift | jumps_rshift
  
  domain <- domain[non_redundant_index]
  fn_vals <- fn_vals[non_redundant_index]
  
  yleft <- max(fn_vals[1], 1e-7)
  yright <- min(fn_vals[length(fn_vals)], -1e-7)
  psi <- approxfun(domain, fn_vals, method = "linear",
                   yleft = yleft, yright = yright)
  
  return(psi)
}

kde_decr_score_est <- function(residuals, quantile = NULL, k = 3000, kernel_pts = 2^15, ...) {
  
  residuals <- sort(residuals)
  
  supp1 = range(residuals) + 3*bw.nrd0(residuals) * c(-1, 1)
  supp2 = median(residuals) + IQR(residuals) * c(-20, 20)
  supp = c(max(supp1[1], supp2[1]), min(supp1[2], supp2[2]))
  
  kde <- density(residuals, n = kernel_pts,
                 from=supp[1], to=supp[2], ...)
  
  n = length(residuals)
  
  k = max(k, 2 * n)
  if (n > 2) {
    skip = 2 * floor(k / n)
  } else {
    skip = 0
  }
  
  iso_res = isotonize_score_given_density(kde$x, kde$y, J_eval = quantile,
                                          k=k, k1=1+skip, k2=k-skip)
  
  return(iso_res)
}

# This can be improved to only calculate cv_spline_score(x) once if both
# spline_df_min and spline_df_1se are called
# score_estimation <- function(x, scoring_type){
#   if (scoring_type == "asm"){
#     return(kde_decr_score_est(x))
#   }
#   
#   if (scoring_type == "spline_df_min"){
#     return(spline_score_estimation(x, df = cv_spline_score(x)$df_min)$rho)
#   }
#   
#   if  (scoring_type == "spline_df_1se"){
#     return(spline_score_estimation(x, df = cv_spline_score(x)$df_1se)$rho)
#   }
# }

score_estimation <- function(x, scoring_types) {
  
  # Ensure scoring_types is a character vector
  scoring_types <- as.character(scoring_types)
  
  # Optimization: Pre-calculate spline CV results *once* if needed.
  cv_res <- NULL 
  if (any(c("spline_df_min", "spline_df_1se") %in% scoring_types)) {
    cv_res <- cv_spline_score(x)
  }
  
  # Use lapply to iterate over each scoring_type
  results_list <- lapply(scoring_types, function(st) {
    
    # switch() is a clean way to handle multiple string conditions
    switch(st,
           "asm" = {
             kde_decr_score_est(x)
           },
           "spline_df_min" = {
             spline_score_estimation(x, df = cv_res$df_min)$rho
           },
           "spline_df_1se" = {
             spline_score_estimation(x, df = cv_res$df_1se)$rho
           },
           "spline_df_2se" = {
             spline_score_estimation(x, df = cv_res$df_2se)$rho
           },
           "spline_df_3se" = {
             spline_score_estimation(x, df = cv_res$df_3se)$rho
           },
           # Default case for any unknown scoring_type
           {
             warning(paste("Unknown scoring_type:", st, "- returning NULL"))
             NULL
           }
    )
  })
  
  # Name the output list with the scoring types for easy access
  names(results_list) <- scoring_types
  
  return(results_list)
}


