library(MASS)
library(ggplot2)
library(reshape2)
library(knitr)
library(latex2exp)
library(patchwork)
library(stringr)
library(viridis)
library(RColorBrewer) 

# FDR-TPR
fdr_tpr_09 <- function(significance, n){
  n0 <- n*0.9
  R  <- sum(significance)
  FP <- sum(significance[seq_len(n0)])
  TP <- sum(significance[(n0 + 1):n])
  FDR <- if (R > 0) FP / R else 0
  TPR <- if ((n-n0) > 0) TP / (n-n0) else NA_real_
  c(FDR = FDR, TPR = TPR)
}


# Oracle P-value with true prior from invchi
cal_p_value_oracle_invchi <- function(z, A, s2, df, K, stdev_unscaled,
                                      oracle_psi, d0 = 10, s0_2 = 1,
                                      mu0 = 20, mu_sd = 2,
                                      grid_n = 2000, tail_prob = 1e-6) {
  stopifnot(df > 0, K > 0, stdev_unscaled > 0, oracle_psi > 0, d0 > 0, s0_2 > 0)

  nu2 <- stdev_unscaled^2  # Z | sigma2, H0 ~ N(0, nu2*sigma2)

  # Prior: tau2 ~ scaled-inv-chi^2(d0, s0_2)
  # sigma2 = oracle_psi * tau2  =>  sigma2 ~ scaled-inv-chi^2(d0, s0_2*oracle_psi)
  s2_scale <- s0_2 * oracle_psi

  # Integration range using prior quantiles of sigma2
  # If X ~ chi^2_{d0}, then sigma2 = d0*s2_scale / X.
  q_low  <- d0 * s2_scale / qchisq(1 - tail_prob, d0)
  q_high <- d0 * s2_scale / qchisq(tail_prob, d0)
  q_low  <- max(q_low, 1e-12)
  q_high <- max(q_high, q_low * 10)

  # log-spaced grid for sigma2
  log_sig2 <- seq(log(q_low), log(q_high), length.out = grid_n)
  sig2 <- exp(log_sig2)

  # ---- log prior density of scaled-inv-chi^2(d0, s2_scale) on sigma2 ----
  # f(x) = ( (d0*s2_scale/2)^{d0/2} / Gamma(d0/2) ) * x^{-(d0/2+1)} * exp( - d0*s2_scale/(2x) )
  log_prior <- (d0/2) * log(d0 * s2_scale / 2) - lgamma(d0/2) -
    (d0/2 + 1) * log(sig2) - (d0 * s2_scale) / (2 * sig2)

  # ---- log likelihood from S^2 | sigma2 ----
  # df*s2/sigma2 ~ chi^2_df  =>  f(s2|sig2) ∝ sig2^{-df/2} * exp(-df*s2/(2*sig2))
  log_like_s2 <- -(df/2) * log(sig2) - (df * s2) / (2 * sig2)

  # ---- integrate out mu: mu ~ N(mu0, mu_sd^2), A | mu,sig2 ~ N(mu, sig2/K)
  # => A | sig2 ~ N(mu0, mu_sd^2 + sig2/K)
  log_like_A <- dnorm(A, mean = mu0, sd = sqrt(mu_sd^2 + sig2 / K), log = TRUE)

  # unnormalized log posterior for sigma2 given (A, s2) under H0
  log_post <- log_prior + log_like_s2 + log_like_A

  # stabilize weights
  m <- max(log_post)
  w <- exp(log_post - m)

  # trapezoid normalization in sigma2-space
  dsig2 <- diff(sig2)
  norm_const <- sum(0.5 * (w[-grid_n] + w[-1]) * dsig2)
  post_w <- w / norm_const  # posterior density approximated on grid

  # tail probability under H0 given sigma2: 2*Phi(-|z|/sqrt(nu2*sigma2))
  p_tail <- 2 * pnorm(-abs(z) / sqrt(nu2 * sig2))

  # integrate p_tail w.r.t posterior density
  p_val <- sum(0.5 * (p_tail[-grid_n] * post_w[-grid_n] +
                      p_tail[-1]      * post_w[-1]) * dsig2)

  p_val <- min(max(p_val, 0), 1)
  return(p_val)
}

P_value_oracle_invchi <- function(info,oracle_psi,d0 = 10, s0_2 = 1,
                                      mu0 = 20, mu_sd = 2,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  n <- info$n
  A <- info$A
  K <- info$K
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  p_set_oracle<- numeric(n)
  
  for(i in 1:n){
    p_value_oracle <- cal_p_value_oracle_invchi(Z[i],A[i],var[i],df,K,stdev.unscaled[i],oracle_psi[i],d0,s0_2,mu0,mu_sd)
    p_set_oracle[i] <- p_value_oracle 
  }
  return(p_set_oracle)
}

# Manorm2 
P_value_manorm2 <- function(data){
  YA <- data$YA
  YB <- data$YB
  
  # MAnorm2
  condA <- bioCond(norm.signal = YA, name = "A")
  condB <- bioCond(norm.signal = YB, name = "B")
  conds <- list(A = condA, B = condB)
  
  conds <- fitMeanVarCurve(
    conds,
    ratio.var = c(A = 1, B = 1),
    method = "local regression",
    occupy.only = FALSE,
    args.locfit = list(maxk = 5000),
    args.lp = list(nn = 0.9), 
    verbose = FALSE
  )
  
  conds <- estimatePriorDf(conds, occupy.only = FALSE)
  
  res <- diffTest(conds$A, conds$B)
  p_set_manorm2 <- res$pval
  return(p_set_manorm2)
}

# Oracle P-value with true prior is discrete
cal_p_value_oracle_discrete <- function(z, A, s2, df, K, stdev_unscaled,
                                        oracle_psi,
                                        value_true = c(1,10), probability_true=c(0.5,0.5),
                                        mu0 = 20, mu_sd = 2) {
  stopifnot(length(value_true) == length(probability_true))
  stopifnot(all(probability_true >= 0))
  stopifnot(abs(sum(probability_true) - 1) < 1e-8)
  stopifnot(df > 0, K > 0, stdev_unscaled > 0, oracle_psi > 0)

  nu2 <- stdev_unscaled^2
  tvals <- as.numeric(value_true)
  p0 <- as.numeric(probability_true)

  # sigma^2 takes discrete values: sigma_m^2 = psi * t_m
  sig2 <- oracle_psi * tvals
  sig2 <- pmax(sig2, 1e-12)  # guard

  # Likelihood pieces under H0 (null-calibrated)
  # 1) S^2 | sigma^2 : df*s2/sigma^2 ~ chi^2_df
  # density: f(s2|sig2) = (df/sig2) * f_chisq(df*s2/sig2)  [Jacobian form]
  # We'll compute in log form for stability:
  x <- df * s2 / sig2
  log_like_s2 <- dchisq(x, df = df, log = TRUE) + log(df) - log(sig2)

  # 2) Integrate out mu: mu ~ N(mu0, mu_sd^2), A|mu,sig2 ~ N(mu, sig2/K)
  # => A | sig2 ~ N(mu0, mu_sd^2 + sig2/K)
  log_like_A <- dnorm(A, mean = mu0, sd = sqrt(mu_sd^2 + sig2 / K), log = TRUE)

  # Posterior weights on each discrete component tau^2 = t_m
  log_w <- log(p0) + log_like_s2 + log_like_A
  m <- max(log_w)
  w <- exp(log_w - m)
  w_sum <- sum(w)

  if (!is.finite(w_sum) || w_sum <= 0) {
    return(1)  # safest fallback
  }
  w_post <- w / w_sum

  # Tail probability for Z under H0 given sigma^2
  p_tail <- 2 * pnorm(-abs(z) / sqrt(nu2 * sig2))

  # Mixture p-value
  p_val <- sum(w_post * p_tail)
  p_val <- min(max(p_val, 0), 1)
  return(p_val)
}

P_value_oracle_discrete <- function(info,oracle_psi,value_true = c(1,10), probability_true=c(0.5,0.5),mu0 = 20, mu_sd = 2,contrast_name=NULL){
  fit_contrasts <- info$fit_contrasts
  df <- info$df
  var <- info$var
  n <- info$n
  A <- info$A
  K <- info$K
  
  if (is.null(contrast_name)){
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,1]
    Z <- fit_contrasts$coefficients[,1]
  }else{
    stdev.unscaled <- fit_contrasts$stdev.unscaled[,contrast_name]
    Z <- fit_contrasts$coefficients[,contrast_name]
  }
  p_set_oracle<- numeric(n)
  
  for(i in 1:n){
    p_value_oracle <- cal_p_value_oracle_discrete(Z[i],A[i],var[i],df,K,stdev.unscaled[i],oracle_psi[i],value_true,probability_true,mu0,mu_sd)
    p_set_oracle[i] <- p_value_oracle 
  }
  return(p_set_oracle)
}

# MAP (modified)
P_value_map_modified <- function(info,contrast_name=NULL){
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
  sigma2_hat <- exp(info$emean)
  
  se_hat <- sqrt(sigma2_hat) * stdev.unscaled
  z <- Z / se_hat

  pval <- 2 * pnorm(-abs(z))
  pval
}

# -----Plot-----
# Simulation plot

make_summary_df <- function(mat, df_val) {
  rm <- rowMeans(mat, na.rm = TRUE)
  if (is.null(names(rm))) names(rm) <- rownames(mat)

  nm <- names(rm)
  data.frame(
    df = df_val,
    metric = sub("^.*\\.", "", nm),            
    method = sub("\\.(FDR|TPR)$", "", nm),      
    value = as.numeric(rm),
    stringsAsFactors = FALSE
  )
}

plot_fdr <- function(df_fdr){
  p_fdr <- ggplot(df_fdr, aes(x = df_f, y = value,
                           group = method, color = method, shape = method)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 1) +
    coord_cartesian(ylim = c(0, 0.1)) +
    labs(x = expression(bold(K-p)), y = "FDR", color = NULL, shape = NULL) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.05, 0.05)))+
    scale_color_discrete(
      breaks = method_levels,
      labels = c("Oracle","t-test",expression(bold("Untrended-Inv"*chi^2)),
                 "Untrended-NPMLE",expression(bold("Reg-Inv"*chi^2)),"Reg-NPMLE","Joint-NPMLE")
    ) +
    scale_shape_manual(
      breaks = method_levels,
      values = c(16, 17, 15, 3, 18, 7, 8),  # pick any 7 distinct shapes you like
      labels = c("Oracle","t-test",expression(bold("Untrended-Inv"*chi^2)), "Untrended-NPMLE",
                 expression(bold("Reg-Inv"*chi^2)),"Reg-NPMLE","Joint-NPMLE")
    ) +
    theme_minimal()
  
  final_fdr <- p_fdr + 
    theme_bw(base_size = 14) +  
    theme(
      plot.title = element_text(face = "bold",hjust = 0.5),
      axis.title.x = element_text(face = "bold",size = 22),  
      axis.title.y = element_text(face = "bold",size = 22),

      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.line = element_blank(),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20)
    )+
    theme(legend.position = "none")
  final_fdr
}

plot_tpr <- function(df_tpr){
  p_tpr <- ggplot(df_tpr, aes(x = df_f, y = value, group = method, color = method,shape=method)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = expression(bold(K-p)), y = "Pow", color = NULL, shape = NULL) +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0.05, 0.05)))+
    scale_color_discrete(
      breaks =method_levels,
      labels = c("Oracle","t-test",expression(bold("Untrended-Inv"*chi^2)),"Untrended-NPMLE",
                 expression(bold("Reg-Inv"*chi^2)),"Reg-NPMLE","Joint-NPMLE")
    ) +
    scale_shape_manual(
      breaks = method_levels,
      values = c(16, 17, 15, 3, 18, 7, 8),  # pick any 7 distinct shapes you like
      labels = c("Oracle","t-test",expression(bold("Untrended-Inv"*chi^2)),"Untrended-NPMLE",
                 expression(bold("Reg-Inv"*chi^2)), "Reg-NPMLE","Joint-NPMLE")
    ) +
    theme_minimal()
  
  final_tpr <- p_tpr + 
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold",hjust = 0.5),
      axis.title.x = element_text(face = "bold",size = 22),  
      axis.title.y = element_text(face = "bold",size = 22),

      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.line = element_blank(),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20)
    )+
    theme(legend.position = "none")
  final_tpr
}
