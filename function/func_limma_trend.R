library(MASS)
library(ggplot2)
library(reshape2)
library(knitr)
library(REBayes)
library(latex2exp)
library(patchwork)
library(stringr)
library(matrixStats)

# BH-procedure
BH_adjust <- function(p_values, alpha = 0.05){
  sorted_indices <- order(p_values)
  sorted_p_values <- p_values[sorted_indices]
  bh_adjusted_p <- p.adjust(sorted_p_values, method = "BH")
  bh_adjusted_p_original_order <- bh_adjusted_p[order(sorted_indices)]
  significance <- bh_adjusted_p_original_order <= alpha
  return(significance)
}

# Information extractor
info_extractor <- function(data,M=NULL){
  Y <- data$Y
  design <- data$design
  contrast <- data$contrast
  
  lmfit <- lmFit(Y, design)
  
  if (!is.null(M)) {
    if (length(M) != nrow(Y)) {
      stop(sprintf("M must have length nrow(Y) = %d, but got %d.", nrow(Y), length(M)))
    }
    lmfit$Amean <- as.numeric(M)
  }
  
  fit2 <- contrasts.fit(lmfit, contrast)
  
  # Extract n,K,p (all scalar)
  n <- dim(Y)[1]
  K <- dim(Y)[2]
  df <- lmfit$df.residual[1]
  
  # Extract summary statistics
  var <- lmfit$sigma^2
  A <- lmfit$Amean
  
  # Fit mean-variance trend
  emean <- fit_trend_var_A(var,A)$emean
  V <- var/exp(emean)
  
  return(list(fit_contrasts = fit2,n=n,K=K,df=df,var=var,A=A,emean=emean,V=V))
}

# Fitted trend

fit_trend_var_A<- function(var, A) {
  y <- log(pmax(var, 1e-12))  # log S^2

  nok <- length(y)
  splinedf <- 1L + (nok >= 3L) + (nok >= 6L) + (nok >= 30L)
  splinedf <- min(splinedf, length(unique(A)))

  if (splinedf < 2L) {
    y_fit <- rep(mean(y), length(y))
    pred_y <- function(a_new) rep(mean(y), length(a_new))
  } else {
    if (!requireNamespace("splines", quietly = TRUE))
      stop("splines package required")
    basis <- splines::ns(A, df = splinedf, intercept = TRUE)
    fit <- lm.fit(basis, y)
    beta <- fit$coefficients
    y_fit <- as.numeric(basis %*% beta)

    pred_y <- function(a_new) {
      basis_new <- splines::ns(
        a_new, df = splinedf, intercept = TRUE,
        Boundary.knots = attr(basis, "Boundary.knots"),
        knots = attr(basis, "knots")
      )
      as.numeric(basis_new %*% beta)
    }
  }

  list(
    emean = y_fit,  # estimate of E[log S^2 | A] at observed A
    pred_e = pred_y  # predictor for new A
  )
}


# Prior estimates

## Parametric estimation for Untrended-Invchi and Reg-Invchi
var_invchi_prior <- function(info){
  fit_contrasts <- info$fit_contrasts
  A <- info$A
  var <- info$var
  
  fit_untrend <- eBayes(fit_contrasts, trend = FALSE)
  fit_trend <- eBayes(fit_contrasts, trend = TRUE)
  
  untrended_invchi_prior <- list(d0 = fit_untrend$df.prior, s02 = fit_untrend$s2.prior[1])
  
  emean <- fit_trend_var_A(var,A)$emean
  reg_invchi_prior <- list(d0 = fit_trend$df.prior,v02=as.numeric(fit_trend$s2.prior[1]/exp(emean)[1]))

  return(list(untrended_invchi_prior=untrended_invchi_prior,
              reg_invchi_prior=reg_invchi_prior))
}

## NPMLE estimation for Untrended-NPMLE and Reg-NPMLE
var_npmle_1d <- function(V, df, v = 30, qL = 0.01,
                          floorA = 1e-300,verbose = TRUE, ...) {
  s <- as.numeric(V)
  stopifnot(is.finite(df), df > 0)
  stopifnot(all(is.finite(s)), all(s >= 0))

  G <- length(s)
  w <- rep(1 / G, G)

  # grid for sigma^2 values
  if (length(v) == 1) {
    v <- seq(log(quantile(s, qL)),log(max(s)),length = v)
    v <- exp(v)
  }

  pv <- length(v)
  dv <- rep(1, pv)
  logv <- log(v)

  if (verbose) cat("Calculating the constraint matrix:")
  t0 <- Sys.time()

  # x_ij = df * s_i / v_j  (G x pv)
  xmat <- df * outer(s, 1 / v)

  # log kernel: log{ dchisq(df*s/v; df) * df/v }
  logAv <- matrix(dchisq(as.vector(xmat), df = df, log = TRUE),
                  nrow = G, ncol = pv, byrow = FALSE) +
           log(df) - rep(logv, each = G)

  # Row-wise scaling for optimization stability
  rowMax <- apply(logAv, 1, max, na.rm = TRUE)
  logAv_scaled <- sweep(logAv, 1, rowMax, "-")

  Av_scaled <- exp(logAv_scaled)
  Av_scaled[!is.finite(Av_scaled)] <- 0
  if (!is.null(floorA) && floorA > 0) Av_scaled <- pmax(Av_scaled, floorA)

  t1 <- Sys.time()
  if (verbose) cat(sprintf(" done. %.2f minutes\n", difftime(t1, t0, units = "mins")))

  if (verbose) cat("Solving for discretized NPMLE:")
  t2 <- Sys.time()
  f <- KWDual(Av_scaled, dv, w, ...)
  t3 <- Sys.time()
  if (verbose) cat(sprintf(" done. %.2f minutes\n", difftime(t3, t2, units = "mins")))

  z <- list(
    grid = v,
    mass =  as.numeric(f$f),
    status = f$status
  )
  class(z) <- "GIGmix"
  z
}

prior_untrended_npmle <- function(info,v = 30,qL = 0.01,floorA = 1e-300,verbose = TRUE,...){
  var <- info$var
  df <- info$df
  
  result <- var_npmle_1d(var, df, v , qL,floorA ,verbose, ...)
  return(result)
}

prior_reg_npmle <- function(info,v = 30,qL = 0.01,floorA = 1e-300,verbose = TRUE,...){
  var <- info$var
  df <- info$df
  A <- info$A
  
  emean <- fit_trend_var_A(var,A)$emean
  V <- var/exp(emean)
  
  result <- var_npmle_1d(V, df, v, qL, floorA,verbose, ...)
  return(result)
}

## NPMLE Estimation for Joint-NPMLE

logAv_func <- function(sig2_k,s,df) {
  B <- length(sig2_k)
  G <- length(s)
  out <- matrix(NA_real_, G, B)

  # if df1 constant, vectorize; else loop rows
  xmat <- outer(df * s, 1 / sig2_k, "*")     # G x B
  out  <- dchisq(xmat, df = df, log = TRUE) +
    log(df) - matrix(log(sig2_k), nrow = G, ncol = B, byrow = TRUE)
  out
}

bin_A <- function(A, nbins = 50) {
  n <- length(A)
  nbins <- min(as.integer(nbins), n)  

  ord <- order(A)
  bin_ord <- ceiling(seq_len(n) * nbins / n)  
  bin_id <- integer(n)
  bin_id[ord] <- bin_ord
  bin_id
}

get_logA_bin<- function(b,mu_list, sig2_vec,t,s,df,K) {
  mu_pts <- mu_list[[b]]
  nb <- length(mu_pts)
  J  <- length(sig2_vec)

  # log f(s^2 | sig2): G x J
  logV <- logAv_func(sig2_vec,s,df)

  # Precompute squared diffs once per bin
  DIFF <- outer(t, mu_pts, "-")   # G x nb
  D2   <- DIFF * DIFF             # G x nb

  # constants for Normal
  c0 <- dnorm(0, log = TRUE)      # -0.5*log(2*pi)
  sd2_vec <- sig2_vec / K      # sd^2 = sig2/K
  const_vec <- c0 - 0.5 * log(sd2_vec) - log(nb)  # length J

  # Accumulate into logV in-place
  for (j in seq_len(J)) {
    logU <- const_vec[j] - 0.5 * (D2 / sd2_vec[j])   # G x nb
    logV[, j] <- logV[, j] + rowLogSumExps(logU)
  }
  logV  # G x J
}

expand_to_mu_sigma <- function(fuv, atom_bin, atom_sig2, mu_list) {
  stopifnot(length(fuv) == length(atom_bin), length(fuv) == length(atom_sig2))
  B <- length(mu_list)

  out <- vector("list", length(fuv))
  for (k in seq_along(fuv)) {
    b <- atom_bin[k]
    mu_pts <- mu_list[[b]]
    nb <- length(mu_pts)

    out[[k]] <- data.frame(
      mu     = mu_pts,
      sigma2 = atom_sig2[k],
      weight = fuv[k] / nb
    )
  }

  do.call(rbind, out)
}

joint_npmle_2d <- function(info,pbin=30,pv=30,qRes_L=0.01,qRes_R=0.99,floorA = 1e-300,verbose = TRUE,...){
  A <- info$A
  var <- info$var
  df <- info$df
  K <- info$K
  
  stopifnot(length(df) == 1, is.finite(df), df > 0)
  stopifnot(length(K)  == 1, is.finite(K),  K  > 0)

  t <- A
  s <- var
  G <- length(t)
  w <- rep(1/G, G)
  
  if (verbose) cat("Calculating the constraint matrix:")
  t0 <- Sys.time()
  
  bin_id <- bin_A(t, nbins = pbin)
  B <- max(bin_id)
  mu_list <- split(t, bin_id)
  
  fitobj <- fit_trend_var_A(var,A)
  emean <-fitobj$emean
  
  r <- log(s) - emean
  rL <- as.numeric(quantile(r, qRes_L, na.rm = TRUE))
  rU <- as.numeric(quantile(r, qRes_R, na.rm = TRUE))
  r_grid <- seq(rL, rU, length.out = pv)
  
  # use the fitted value of the median A of bin
  u_bin <- tapply(t, bin_id, median)   # length B
  e_bin <- as.numeric(fitobj$pred_e(u_bin)) 
  
  sig2_list <- lapply(seq_len(B), function(b) {
    sig2_vec <- exp(e_bin[b] + r_grid)             # trend-centered grid
    sig2_vec <- pmin(pmax(sig2_vec, min(s)), max(s)) # clip to global range (as you proposed)
    sig2_vec <- pmax(sig2_vec, 1e-12)              # guard
    sort(unique(sig2_vec))
  })
    
  atom_bin  <- rep(seq_len(B), times = vapply(sig2_list, length, 1L))
  atom_sig2 <- unlist(sig2_list, use.names = FALSE)
  M <- length(atom_sig2)
  duv <- rep(1, M)
  
  # column ranges per bin
  J_by_bin  <- vapply(sig2_list, length, 1L)
  start_col <- c(1L, 1L + cumsum(J_by_bin))[seq_len(B)]
  
  # BIG STORAGE: log-likelihood matrix (G x M)
  logA_all <- matrix(NA_real_, nrow = G, ncol = M)
  
  # row-wise max (for stable exp)
  row_max <- rep(-Inf, G)
  
  for (b in seq_len(B)) {
    sig2_vec <- sig2_list[[b]]
    Jb <- length(sig2_vec)
    cols <- start_col[b]:(start_col[b] + Jb - 1L)
  
    logA_bin <- get_logA_bin(b, mu_list,sig2_vec,t,s,df,K)     # G x Jb
  
    logA_all[, cols] <- logA_bin
    row_max <- pmax(row_max, rowMaxs(logA_bin))   
  }
  A_scaled <- exp(sweep(logA_all, 1, row_max, "-"))
  A_scaled <- pmax(A_scaled, floorA)
  
  rm(logA_all)
  
  t1 <- Sys.time()
  if (verbose) cat(sprintf(" done. %.2f minutes\n", difftime(t1, t0, units = "mins")))

  if (verbose) cat("Solving for discretized NPMLE:")
  t2 <- Sys.time()
  
  f <- KWDual(A_scaled, duv, w)
  
  t3 <- Sys.time()
  if (verbose) cat(sprintf(" done. %.2f minutes\n", difftime(t3, t2, units = "mins")))

  fuv <- f$f
  
  grid_mu_sigma <- expand_to_mu_sigma(fuv, atom_bin, atom_sig2, mu_list)
  z <- list(
    grid    = as.matrix(grid_mu_sigma[, c("mu", "sigma2")]),
    mass    = grid_mu_sigma$weight,
    status  = f$status,
    bin_id  = bin_id
  )
  return(z)
}


# P-value calculation

## T-test
cal_p_value_t <- function(Z,S,d,stdev.unscaled){
  stat <- Z/(stdev.unscaled*sqrt(S))
  p_value <- 2*pt(-abs(stat),d)
  return(p_value)
}

P_value_t <- function(info,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  n <- info$n
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }

  p_set_t <- numeric(n)
  
  for(i in 1:n){
    p_value_t <- cal_p_value_t(Z[i],var[i],df,stdev.unscaled[i])
    p_set_t[i] <- p_value_t
  }
  
  return(p_set_t)
}

## Untrended-Invchi
P_value_untrended_invchi <- function(info,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  
  fit_untrend <- eBayes(fit_contrasts, trend = FALSE)
  
  if (is.null(contrast_name)){
    p_set_untrended_invchi <- fit_untrend$p.value[, 1]
  }else{
    p_set_untrended_invchi <- fit_untrend$p.value[, contrast_name]
  }
  
  return(p_set_untrended_invchi)
}

## Untrended-NPMLE
cal_p_value_untrended_npmle <- function(value,probability,Z,S,d,stdev.unscaled){
  prob <- probability * (value) ** (-d/2) * exp(-d * S/(2*value))
  k_prime <- prob / sum(prob)
  
  p_value <- 2* sum(k_prime * pnorm(-abs(Z), mean = 0, sd = sqrt(value)*stdev.unscaled, lower.tail = TRUE))
  return(p_value)
}

P_value_untrended_npmle <- function(info,grid,mass,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  n <- info$n
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  
  p_set_untrended_npmle <- numeric(n)

  for(i in 1:n){
    p_value_untrended_npmle <- cal_p_value_untrended_npmle(grid,mass,Z[i],var[i],df,stdev.unscaled[i])
    p_set_untrended_npmle[i] <- p_value_untrended_npmle
  }
  
  return(p_set_untrended_npmle)
}

## Reg-Invchi

P_value_reg_invchi <- function(info,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  
  fit_trend <- eBayes(fit_contrasts, trend = TRUE)
  
  if (is.null(contrast_name)){
    p_set_reg_invchi <- fit_trend$p.value[, 1]
  }else{
    p_set_reg_invchi <- fit_trend$p.value[, contrast_name]
  }
  
  return(p_set_reg_invchi)
}

## Reg-NPMLE

cal_p_value_reg_npmle <- function(value,probability,Z,V,df,stdev.unscaled,emean){
  prob <- probability * (value) ** (-df/2) * exp(-df * V/(2*value))
  k_prime <- prob / sum(prob)
  
  p_value <- 2* sum(k_prime * pnorm(-abs(Z), mean = 0, sd = sqrt(value * exp(emean)) * stdev.unscaled, lower.tail = TRUE))
  return(p_value)
}

P_value_reg_npmle <- function(info,grid,mass,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  A <- info$A
  n <- info$n
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  
  emean <- info$emean
  V <- var/exp(emean)
  
  p_set_reg_npmle <- numeric(n)
  for(i in 1:n){
    p_value_reg_npmle <- cal_p_value_reg_npmle(grid,mass,Z[i],V[i],df,stdev.unscaled[i],emean[i])
    p_set_reg_npmle[i] <- p_value_reg_npmle
  }
  return(p_set_reg_npmle)
}

## Joint-NPMLE

cal_p_value_joint_npmle <- function(value,probability,Z,A,var,df,K,stdev.unscaled){
  
  I <- length(probability)
  
  A_value <- value[,1]
  var_value <- value[,2]
  
  temp <- A - A_value
  
  prob <- probability * var_value ** (-(df+1)/2) * exp(-K*temp**2/(2*var_value)) * exp(-df*var/(2*var_value))
  k_prime <- prob / sum(prob)
  
  p_value <- 2* sum(k_prime * pnorm(-abs(Z), mean = 0, sd = sqrt(var_value) * stdev.unscaled, lower.tail = TRUE))
  return(p_value)
}

P_value_joint_npmle <- function(info,grid,mass,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  A <- info$A
  K <- info$K
  n <- info$n
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  
  p_set_joint_npmle <- numeric(n)
  for(i in 1:n){
    p_value_joint_npmle <- cal_p_value_joint_npmle(grid,mass,Z[i],A[i],var[i],df,K,stdev.unscaled[i])
    p_set_joint_npmle[i] <- p_value_joint_npmle
  }
  return(p_set_joint_npmle)
}



## Wrap-up function
p_value_calculator <- function(info,prior_result,contrast_name=NULL,method_list=c(1,2,3,4,5,6),verbose=FALSE){
  method_name <- c('t_test','untrended_invchi','untrended_npmle','reg_invchi','reg_npmle','joint_npmle')
  P_list_t <- NA
  P_list_untrended_invchi <- NA
  P_list_untrended_npmle <- NA
  P_list_reg_invchi <- NA
  P_list_reg_npmle <- NA
  P_list_joint_npmle <- NA
  if(verbose) print('Calculation Starts')
  result <- list()
  for (m in method_list) {
    
    if (m == 1) {
      if(verbose) print(paste('start of', method_name[m]))
      P_list_t = P_value_t(info,contrast_name)
      result$t_test = P_list_t
    }
    
    if (m == 2) {
      if(verbose) print(paste('start of', method_name[m]))
      P_list_untrended_invchi = P_value_untrended_invchi(info,contrast_name)
      result$untrended_invchi=P_list_untrended_invchi
    }
    
    if (m == 3) {
      if(verbose) print(paste('start of', method_name[m]))
      grid <- prior_result$untrended_npmle$grid
      mass <- prior_result$untrended_npmle$mass
      P_list_untrended_npmle = P_value_untrended_npmle(info,grid,mass,contrast_name)
      result$untrended_npmle =P_list_untrended_npmle
    }
    
    if (m == 4) {
      if(verbose) print(paste('start of', method_name[m]))
      P_list_reg_invchi = P_value_reg_invchi(info,contrast_name)
      result$reg_invchi = P_list_reg_invchi
    }
    
    if (m == 5) {
      if(verbose) print(paste('start of', method_name[m]))
      grid <- prior_result$reg_npmle$grid
      mass <- prior_result$reg_npmle$mass
      P_list_reg_npmle = P_value_reg_npmle(info,grid,mass,contrast_name)
      result$reg_npmle = P_list_reg_npmle
    }
    if (m == 6) {
      if(verbose) print(paste('start of', method_name[m]))
      idx <- prior_result$joint_npmle$idx
      grid <- prior_result$joint_npmle$grid[idx,]
      mass <- prior_result$joint_npmle$mass[idx]
      P_list_joint_npmle = P_value_joint_npmle(info,grid,mass,contrast_name)
      result$joint_npmle = P_list_joint_npmle
    }
  }
  if(verbose) print('Calculation Ends')
  return(result)
}





# Mimic-Joint-NPMLE
prior_mimic_joint_npmle <- function(info,pep_grp,v = 30, qL = 0.01,
                          floorA = 1e-300,verbose = TRUE,...){
  df <- info$df
  var <- info$var
  prior_byM <- vector("list", length(levels(pep_grp)))
  names(prior_byM) <- levels(pep_grp)
  
  for (g in levels(pep_grp)) {
    idx <- which(pep_grp == g)
  
    # stratum-specific prior for sigma^2
    result <- var_npmle_1d(var[idx], df, v, qL, floorA, verbose,...)
  
    value_g <- result$grid
    prob_g  <- result$mass
  
    # save prior for this stratum
    prior_byM[[g]] <- data.frame(
      pep_grp = g,
      value   = value_g,
      prob    = prob_g,
      n_in_grp = length(idx)
    )
  }
  return(prior_byM)
}

P_value_mimic_joint_npmle <- function(info,mimic_joint_npmle_prior,pep_grp,contrast_name=NULL){
  n <- info$n
  df <- info$df
  var <- info$var
  fit_contrasts <- info$fit_contrasts
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  
  p_set_mimic_joint_npmle <- numeric(n)
  for (g in levels(pep_grp)) {
    idx <- which(pep_grp == g)
    p_set_mimic_joint_npmle[idx] <- mapply(
      FUN = function(Zi, Si2, dfi, su) {
        cal_p_value_untrended_npmle(mimic_joint_npmle_prior[[g]]$value, mimic_joint_npmle_prior[[g]]$prob, Zi, Si2, df, su)
      },
      Zi  = Z[idx],
      Si2 = var[idx],
      su  = stdev.unscaled[idx]
    )
  }
  
  return(p_set_mimic_joint_npmle)
}