library(ggplot2)
library(reshape2)
library(progress)
library(dplyr)
library(REBayes)
library(MAnorm2)

# -----Load Function-----
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
source('../function/func_simulation.R')
# -----Simulation Function-----

simulate_data_two_group <- function(n = 10000, K_A = 2, K_B = 30, d0 = 10, s0_2 = 1, mu0 = 20, mu_sd = 0.2, 
  trend_drop = 6, trend_width = 0.15, trend_center =20, sigma2_scale = 1){

  K <- K_A + K_B
  X1 <- c(rep(1, K_A), rep(0, K_B))
  X  <- cbind(1, X1)
  group <- factor(ifelse(X1 == 1, "A", "B"), levels = c("A","B"))
  c_theta <- c(0,1)
  
  n0 <- floor(n * 0.9)
  is_alt <- rep(FALSE, n)
  is_alt[(n0 + 1):n] <- TRUE  # last 10% are non-null

 
  beta1_i <- rnorm(n, mu0, mu_sd)
  tau2 <- 1 / (rchisq(n, d0) / (d0 * s0_2))

  log_psi <- log(s0_2) - trend_drop * plogis((beta1_i - trend_center) / trend_width)
  psi <- exp(log_psi) * sigma2_scale

  sigma2 <- psi * tau2

 
  beta02_i <- rep(0, n)
  beta02_i[is_alt] <- abs(rnorm(sum(is_alt), mean = 0, sd = 4 * sqrt(sigma2[is_alt])))

  beta <- cbind(beta1_i, beta02_i)
  mu <- beta %*% t(X)

  Y <- matrix(NA, nrow = n, ncol = K)
  for (i in 1:n) {
    Y[i, ] <- rnorm(K, mean = mu[i, ], sd = sqrt(sigma2[i]))
  }

  YA <- Y[, group == "A", drop = FALSE]
  YB <- Y[, group == "B", drop = FALSE]

  list(
    design = X, group = group,
    Y = Y, YA = YA, YB = YB,
    beta1 = beta1_i, beta02 = beta02_i,
    is_alt = is_alt,
    oracle_psi = psi, sigma2 = sigma2,
    contrast = c_theta
  )
}


run_comparison <- function(simulate_times,n,K_A = 2, K_B = 30, d0 = 10, s0_2 = 1, mu0 = 20, mu_sd = 0.2,trend_drop = 6, trend_width = 0.15,
                        trend_center =20, sigma2_scale = 1,include_joint = FALSE, verbose = FALSE){
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)

  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_two_group(n,K_A,K_B,d0,s0_2,mu0,mu_sd,trend_drop,trend_width,trend_center,sigma2_scale)
    info <- info_extractor(data)
    ## Prior estimation
    param_prior <- var_invchi_prior(info)
    untrended_invchi_prior <- param_prior[[1]]
    reg_invchi_prior <- param_prior[[2]]
    
    untrended_npmle_prior <- prior_untrended_npmle(info,v=300,verbose=verbose)
    reg_npmle_prior <- prior_reg_npmle(info,v=300,verbose=verbose)
    
    prior_result <- list(
      untrended_invchi=untrended_invchi_prior,
      untrended_npmle=untrended_npmle_prior,
      reg_invchi =reg_invchi_prior,
      reg_npmle=reg_npmle_prior
    )
    if(include_joint){
      joint_npmle_prior <- joint_npmle_2d(info,pbin=50,pv=50,qRes_L = 0.1,verbose=verbose)
      threshold <- 1e-8
      idx <- which(joint_npmle_prior$mass>threshold)
      joint_npmle_prior$idx <- idx
      prior_result$joint_npmle <- joint_npmle_prior
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5,6),verbose=verbose)
    }else{
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5),verbose=verbose)
    }
    
    result$oracle <- P_value_oracle_invchi(info,data$oracle_psi,d0 = d0, s0_2 = s0_2, mu0 = mu0,mu_sd=mu_sd)
    result$manorm2 <- P_value_manorm2(data)
    result$map <- P_value_map_modified(info)
    
    significance <- lapply(result, BH_adjust, alpha = 0.05)
    fdr_tpr_result <-  lapply(significance,fdr_tpr_09,n=n)
    
    v_list <- unlist(fdr_tpr_result) 
    a <- data.frame(value = as.numeric(v_list),row.names = names(v_list))
    setTxtProgressBar(pb, i)
    a
  })
  mat <- do.call(cbind, cols)
  colnames(mat) <- paste0("sim", seq_along(sim_idx))
  mat
}

# -----Run-----
set.seed(1)
simulate_times <- 100
n <- 10000
result_compare <- run_comparison(simulate_times, n,2,10,include_joint = TRUE)
rowMeans(result_compare)

#save(result_compare,file="./data/sim_result_compare.RData")

