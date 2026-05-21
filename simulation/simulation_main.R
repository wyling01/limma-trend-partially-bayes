library(readxl)
library(quadprog)
library(ruv)
library(ggplot2)
library(reshape2)
library(progress)
library(REBayes)
library(VennDiagram)
library(data.table)

# -----Load Function-----
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
source('../function/func_simulation.R')



# -----Setting 1&2-----
## ----Simulation data for Scaled inverse prior----

simulate_data_invchi <- function(n,K, d0 = 10, s0_2 = 1, mu0 = 20, mu_sd = 2, 
                                 trend_drop = 5, trend_width = 2, trend_center =18, sigma2_scale = 5){
  
  X1 <- c(rep(1,K/2),rep(0,K/2))
  X <- cbind(rep(1,K),X1)
  group <- factor(ifelse(X1 == 1, "A", "B"), levels = c("A", "B"))
  
  c_theta <- c(0,1)
  
  beta1_i <- rnorm(n,mu0,mu_sd)
  tau2 <- 1/(rchisq(n, d0)/(d0*s0_2))
  
  log_psi <-  - trend_drop * plogis((beta1_i - trend_center) / trend_width)
  psi <- exp(log_psi) * sigma2_scale
  sigma2 <- psi * tau2
  
  effects_2 <- rnorm(n,0,sqrt(16*sigma2))
  beta02_i <- c(rep(0,n*0.9),effects_2[(n*0.9+1):n])
  beta_2  <- cbind(beta1_i,beta02_i)
  
  mu <- beta_2 %*% t(X)
  Y <- matrix(NA, nrow = n, ncol = K)
  for (g in 1:n) {
    Y[g, ] <- rnorm(K, mean = mu[g, ], sd = sqrt(sigma2[g]))
  }
  YA <- Y[, group == "A", drop = FALSE]
  YB <- Y[, group == "B", drop = FALSE]
  
  list(
    Y=Y,YA = YA,YB = YB,design=X,contrast=c_theta,oracle_psi=psi,sigma2=sigma2
  )
}

run_simulation <- function(simulate_times,simulate_data_invchi=simulate_data_invchi,n=10000,K=4,
                           d0=10,s0_2 = 1, mu0 = 20, mu_sd = 2,
                           trend_drop = 5, trend_width = 2, trend_center =18, sigma2_scale = 5,
                           include_joint = FALSE, pbin=50,pv=50,qRes_L = 0.01,verbose = FALSE,seed=1){
  set.seed(seed)
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_invchi(n,K,d0,s0_2,mu0,mu_sd,trend_drop,trend_width, trend_center, sigma2_scale)
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
      joint_npmle_prior <- joint_npmle_2d(info,pbin=pbin,pv=pv,qRes_L = qRes_L,verbose=verbose)
      threshold <- 1e-8
      idx <- which(joint_npmle_prior$mass>threshold)
      joint_npmle_prior$idx <- idx
      prior_result$joint_npmle <- joint_npmle_prior
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5,6),verbose=verbose)
    }else{
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5),verbose=verbose)
    }
    
    result$oracle <- P_value_oracle_invchi(info,data$oracle_psi,d0 = d0, s0_2 = s0_2, mu0 = mu0,mu_sd=mu_sd)
    
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
    result$manorm2 <- P_value_manorm2(data)
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
    
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

simulate_data_invchi_notrend <- function(n,K, d0 = 10, s0_2 = 1, mu0 = 20, mu_sd = 2){
  
  X1 <- c(rep(1,K/2),rep(0,K/2))
  X <- cbind(rep(1,K),X1)
  group <- factor(ifelse(X1 == 1, "A", "B"), levels = c("A", "B"))
  
  c_theta <- c(0,1)
  
  beta1_i <- rnorm(n,mu0,mu_sd)
  tau2 <- 1/(rchisq(n, d0)/(d0*s0_2))
  
  psi <- rep(1,n)
  sigma2 <- psi * tau2
  
  effects_2 <- rnorm(n,0,sqrt(16*sigma2))
  beta02_i <- c(rep(0,n*0.9),effects_2[(n*0.9+1):n])
  beta_2  <- cbind(beta1_i,beta02_i)
  
  mu <- beta_2 %*% t(X)
  Y <- matrix(NA, nrow = n, ncol = K)
  for (g in 1:n) {
    Y[g, ] <- rnorm(K, mean = mu[g, ], sd = sqrt(sigma2[g]))
  }
  
  YA <- Y[, group == "A", drop = FALSE]
  YB <- Y[, group == "B", drop = FALSE]
  
  list(
    Y=Y,YA = YA,YB = YB,design=X,contrast=c_theta,oracle_psi=psi,sigma2=sigma2
  )
}

run_simulation_notrend <- function(simulate_times,n=10000,K=4,
                                   d0=10,s0_2 = 1, mu0 = 20, mu_sd = 2,
                                   include_joint = FALSE, pbin=50,pv=50,qRes_L = 0.01,verbose = FALSE,seed=1){
  set.seed(seed)
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_invchi_notrend(n,K,d0,s0_2,mu0,mu_sd)
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
      joint_npmle_prior <- joint_npmle_2d(info,pbin=pbin,pv=pv,qRes_L = qRes_L,verbose=verbose)
      threshold <- 1e-8
      idx <- which(joint_npmle_prior$mass>threshold)
      joint_npmle_prior$idx <- idx
      prior_result$joint_npmle <- joint_npmle_prior
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5,6),verbose=verbose)
    }else{
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5),verbose=verbose)
    }
    
    result$oracle <- P_value_oracle_invchi(info,data$oracle_psi,d0 = d0, s0_2 = s0_2, mu0 = mu0,mu_sd=mu_sd)
    
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
    result$manorm2 <- P_value_manorm2(data)
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
    
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

## Simulation data for Dirac and Two-point prior
simulate_data_discrete <- function(n,K, value_true = c(1,10), probability_true=c(0.5,0.5), mu0 = 20, mu_sd = 2, 
                                   trend_drop = 5, trend_width = 2, trend_center =18, sigma2_scale = 5){
  
  X1 <- c(rep(1,K/2),rep(0,K/2))
  X <- cbind(rep(1,K),X1)
  group <- factor(ifelse(X1 == 1, "A", "B"), levels = c("A", "B"))
  
  c_theta <- c(0,1)
  
  beta1_i <- rnorm(n,mu0,mu_sd)
  tau2 <- sample(value_true, n, replace = TRUE, prob = probability_true)
  
  log_psi <-  - trend_drop * plogis((beta1_i - trend_center) / trend_width)
  psi <- exp(log_psi) * sigma2_scale
  sigma2 <- psi * tau2
  
  effects_2 <- rnorm(n,0,sqrt(16*sigma2))
  beta02_i <- c(rep(0,n*0.9),effects_2[(n*0.9+1):n])
  beta_2  <- cbind(beta1_i,beta02_i)
  
  mu <- beta_2 %*% t(X)
  Y <- matrix(NA, nrow = n, ncol = K)
  for (g in 1:n) {
    Y[g, ] <- rnorm(K, mean = mu[g, ], sd = sqrt(sigma2[g]))
  }
  
  YA <- Y[, group == "A", drop = FALSE]
  YB <- Y[, group == "B", drop = FALSE]
  
  list(
    Y=Y,YA = YA,YB = YB,design=X,contrast=c_theta,oracle_psi=psi,sigma2=sigma2
  )
}

run_simulation_discrete <- function(simulate_times,simulate_data_discrete=simulate_data_discrete,n=10000,K=4,
                                    value_true = c(1,10), probability_true=c(0.5,0.5), mu0 = 20, mu_sd = 2,
                                    trend_drop = 5, trend_width = 2, trend_center =18, sigma2_scale = 5,
                                    include_joint = FALSE, pbin=50,pv=50,qRes_L = 0.01,verbose = FALSE,seed=1){
  set.seed(seed)
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_discrete(n,K,value_true,probability_true,mu0,mu_sd,trend_drop,trend_width, trend_center, sigma2_scale)
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
      joint_npmle_prior <- joint_npmle_2d(info,pbin=pbin,pv=pv,qRes_L = qRes_L,verbose=verbose)
      threshold <- 1e-8
      idx <- which(joint_npmle_prior$mass>threshold)
      joint_npmle_prior$idx <- idx
      prior_result$joint_npmle <- joint_npmle_prior
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5,6),verbose=verbose)
    }else{
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5),verbose=verbose)
    }
    
    result$oracle <- P_value_oracle_discrete(info,data$oracle_psi,value_true = value_true, probability_true=probability_true, mu0 = mu0,mu_sd=mu_sd)
    
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
    result$manorm2 <- P_value_manorm2(data)
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
    
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

simulate_data_discrete_notrend <- function(n,K, value_true = c(1,10), probability_true=c(0.5,0.5), mu0 = 20, mu_sd = 2){
  
  X1 <- c(rep(1,K/2),rep(0,K/2))
  X <- cbind(rep(1,K),X1)
  group <- factor(ifelse(X1 == 1, "A", "B"), levels = c("A", "B"))
  
  c_theta <- c(0,1)
  
  beta1_i <- rnorm(n,mu0,mu_sd)
  tau2 <- sample(value_true, n, replace = TRUE, prob = probability_true)
  
  psi <- rep(1,n)
  sigma2 <- psi * tau2
  
  effects_2 <- rnorm(n,0,sqrt(16*sigma2))
  beta02_i <- c(rep(0,n*0.9),effects_2[(n*0.9+1):n])
  beta_2  <- cbind(beta1_i,beta02_i)
  
  mu <- beta_2 %*% t(X)
  Y <- matrix(NA, nrow = n, ncol = K)
  for (g in 1:n) {
    Y[g, ] <- rnorm(K, mean = mu[g, ], sd = sqrt(sigma2[g]))
  }
  YA <- Y[, group == "A", drop = FALSE]
  YB <- Y[, group == "B", drop = FALSE]
  
  list(
    Y=Y,YA = YA,YB = YB,design=X,contrast=c_theta,oracle_psi=psi,sigma2=sigma2
  )
}

run_simulation_discrete_notrend <- function(simulate_times,n=10000,K=4,
                                            value_true = c(1,10), probability_true=c(0.5,0.5), mu0 = 20, mu_sd = 2,
                                            include_joint = FALSE, pbin=50,pv=50,qRes_L = 0.01,verbose = FALSE,seed=1){
  set.seed(seed)
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_discrete_notrend(n,K,value_true,probability_true,mu0,mu_sd)
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
      joint_npmle_prior <- joint_npmle_2d(info,pbin=pbin,pv=pv,qRes_L = qRes_L,verbose=verbose)
      threshold <- 1e-8
      idx <- which(joint_npmle_prior$mass>threshold)
      joint_npmle_prior$idx <- idx
      prior_result$joint_npmle <- joint_npmle_prior
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5,6),verbose=verbose)
    }else{
      result <- p_value_calculator(info,prior_result,method_list=c(1,2,3,4,5),verbose=verbose)
    }
    
    result$oracle <- P_value_oracle_discrete(info,data$oracle_psi,value_true = value_true, probability_true=probability_true, mu0 = mu0,mu_sd=mu_sd)
    
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
    result$manorm2 <- P_value_manorm2(data)
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
    
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


## -----Run Simulation-----
simulate_times <- 100
n <- 10000

## Scaled inverse chi^2
mat_invchi_type1_K4 <- run_simulation_notrend(simulate_times,n =10000,K=4,d0=10,mu_sd=3,include_joint = TRUE)
mat_invchi_type1_K6 <- run_simulation_notrend(simulate_times,n=10000,K=6,d0=10,mu_sd=3,include_joint = TRUE)
mat_invchi_type1_K10 <- run_simulation_notrend(simulate_times,n=10000,K=10,d0=10,mu_sd=3,include_joint = TRUE)
mat_invchi_type1_K18 <- run_simulation_notrend(simulate_times,n=10000,K=18,d0=10,mu_sd=3,include_joint = TRUE)


mat_invchi_type2_K4 <- run_simulation(simulate_times,simulate_data_invchi,n=10000,K=4,d0=10,mu_sd=3,
                                      trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                      include_joint = TRUE)
mat_invchi_type2_K6 <- run_simulation(simulate_times,simulate_data_invchi,n=10000,K=6,d0=10,mu_sd=3,
                                      trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                      include_joint = TRUE)
mat_invchi_type2_K10 <- run_simulation(simulate_times,simulate_data_invchi,n=10000,K=10,d0=10,mu_sd=3,
                                      trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                      include_joint = TRUE)
mat_invchi_type2_K18 <- run_simulation(simulate_times,simulate_data_invchi,n=10000,K=18,d0=10,mu_sd=3,
                                      trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                      include_joint = TRUE)

## Dirac
mat_dirac_type1_K4 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=4,value_true = c(1,1), 
                                                      probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_dirac_type1_K6 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=6,value_true = c(1,1), 
                                                      probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_dirac_type1_K10 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=10,value_true = c(1,1), 
                                                      probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_dirac_type1_K18 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=18,value_true = c(1,1), 
                                                      probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)

mat_dirac_type2_K4 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=4,
                                              value_true = c(1,1), probability_true=c(0.5,0.5),mu_sd = 3,
                                              trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                              include_joint = TRUE)
mat_dirac_type2_K6 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=6,
                                              value_true = c(1,1), probability_true=c(0.5,0.5),mu_sd = 3,
                                              trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                              include_joint = TRUE)
mat_dirac_type2_K10 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=10,
                                              value_true = c(1,1), probability_true=c(0.5,0.5),mu_sd = 3,
                                              trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                              include_joint = TRUE)
mat_dirac_type2_K18 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=18,
                                              value_true = c(1,1), probability_true=c(0.5,0.5),mu_sd = 3,
                                              trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                              include_joint = TRUE)

## Two-point
mat_twopoint_type1_K4 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=4,value_true = c(1,10), 
                                                         probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_twopoint_type1_K6 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=6,value_true = c(1,10), 
                                                         probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_twopoint_type1_K10 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=10,value_true = c(1,10), 
                                                         probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)
mat_twopoint_type1_K18 <- run_simulation_discrete_notrend(simulate_times,n=10000,K=18,value_true = c(1,10), 
                                                         probability_true=c(0.5,0.5),mu_sd = 3,include_joint = TRUE)

mat_twopoint_type2_K4 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=4,
                                                 value_true = c(1,10), probability_true=c(0.5,0.5),mu_sd = 3,
                                                 trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                                 include_joint = TRUE)
mat_twopoint_type2_K6 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=6,
                                                 value_true = c(1,10), probability_true=c(0.5,0.5),mu_sd = 3,
                                                 trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                                 include_joint = TRUE)
mat_twopoint_type2_K10 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=10,
                                                 value_true = c(1,10), probability_true=c(0.5,0.5),mu_sd = 3,
                                                 trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                                 include_joint = TRUE)
mat_twopoint_type2_K18 <- run_simulation_discrete(simulate_times,simulate_data_discrete,n=10000,K=18,
                                                 value_true = c(1,10), probability_true=c(0.5,0.5),mu_sd = 3,
                                                 trend_drop = 4, trend_width = 4, trend_center =16, sigma2_scale = 12,
                                                 include_joint = TRUE)



#---- Setting 3 ----
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

## -----Run-----
set.seed(1)
simulate_times <- 100
n <- 10000
mat_unbalanced <- run_comparison(simulate_times, n,2,10,include_joint = TRUE)
rowMeans(mat_unbalanced)


# ---- Setting4 ----
simulate_data_two_group_perfect_trend <- function(n = 10000, K_A = 2, K_B = 30, mu0 = 20, mu_sd = 0.2, 
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
  tau2 <- rep(1, n)
  
  log_psi <- - trend_drop * plogis((beta1_i - trend_center) / trend_width)
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
    Y = Y, YA = YA, YB = YB,M=beta1_i,
    beta1 = beta1_i, beta02 = beta02_i,
    is_alt = is_alt,
    oracle_psi = psi, sigma2 = sigma2,
    contrast = c_theta
  )
}

run_comparison_perfect_trend <- function(simulate_times,n,K_A = 2, K_B = 30, mu0 = 20, mu_sd = 0.2,trend_drop = 6, trend_width = 0.15,
                                         trend_center =20, sigma2_scale = 1,include_joint = FALSE, verbose = FALSE){
  sim_idx <- seq_len(simulate_times)
  pb <- txtProgressBar(min = 0, max = simulate_times, style = 3)
  on.exit(close(pb), add = TRUE)
  
  cols <- lapply(seq_along(sim_idx), function(i) {
    data <- simulate_data_two_group_perfect_trend(n,K_A,K_B,mu0,mu_sd,trend_drop,trend_width,trend_center,sigma2_scale)
    info <- info_extractor(data, M = data$M)
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
    
    result$oracle <- P_value_oracle_discrete(info,data$oracle_psi,value_true =c(1,1), probability_true=c(0.5,0.5), mu0 = mu0,mu_sd=mu_sd)
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

## -----Run-----
set.seed(1)
simulate_times <- 100
n <- 10000
mat_external_side <- run_comparison_perfect_trend(simulate_times, n,3,5,trend_width = 0.15,include_joint = FALSE)
rowMeans(mat_external_side)


# Save simulation results
sim_results <- ls(envir = .GlobalEnv, pattern = "^mat")
save(list =sim_results,file="./data/sim_result_all.RData")

# ---- Plot ----
load("./data/sim_result_all.RData")

method_levels <- c("oracle","t_test","untrended_invchi","untrended_npmle",
                   "reg_invchi","reg_npmle","joint_npmle","manorm2")

## Sclaed inverse chi^2 Type 1
sum_df_invchi_type1 <- rbind(
  make_summary_df(mat_invchi_type1_K4, 2, method_levels),
  make_summary_df(mat_invchi_type1_K6, 4, method_levels),
  make_summary_df(mat_invchi_type1_K10, 8, method_levels),
  make_summary_df(mat_invchi_type1_K18, 16, method_levels)
)

sum_df_invchi_type1$method <- factor(sum_df_invchi_type1$method, levels = method_levels)

df_fdr <- subset(sum_df_invchi_type1, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_invchi_type1, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_invchi_type1 <- plot_fdr(df_fdr) +
  labs(title = expression(bold("Scaled Inv" * chi^2))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 30)) 
plot_tpr_invchi_type1 <- plot_tpr(df_tpr)

## Sclaed inverse chi^2 Type 2
sum_df_invchi_type2 <- rbind(
  make_summary_df(mat_invchi_type2_K4, 2, method_levels),
  make_summary_df(mat_invchi_type2_K6, 4, method_levels),
  make_summary_df(mat_invchi_type2_K10, 8, method_levels),
  make_summary_df(mat_invchi_type2_K18, 16, method_levels)
)

sum_df_invchi_type2$method <- factor(sum_df_invchi_type2$method, levels = method_levels)

df_fdr <- subset(sum_df_invchi_type2, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_invchi_type2, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_invchi_type2 <- plot_fdr(df_fdr) +
  labs(title = expression(bold("Scaled Inv" * chi^2))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 30)) 
plot_tpr_invchi_type2 <- plot_tpr(df_tpr)

## Dirac Type 1
sum_df_dirac_type1 <- rbind(
  make_summary_df(mat_dirac_type1_K4, 2, method_levels),
  make_summary_df(mat_dirac_type1_K6, 4, method_levels),
  make_summary_df(mat_dirac_type1_K10, 8, method_levels),
  make_summary_df(mat_dirac_type1_K18, 16, method_levels)
)

sum_df_dirac_type1$method <- factor(sum_df_dirac_type1$method, levels = method_levels)

df_fdr <- subset(sum_df_dirac_type1, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_dirac_type1, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_dirac_type1 <- plot_fdr(df_fdr) +
  labs(title = "Dirac") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size=30)) 
plot_tpr_dirac_type1 <- plot_tpr(df_tpr)

## Dirac Type 2
sum_df_dirac_type2 <- rbind(
  make_summary_df(mat_dirac_type2_K4, 2, method_levels),
  make_summary_df(mat_dirac_type2_K6, 4, method_levels),
  make_summary_df(mat_dirac_type2_K10, 8, method_levels),
  make_summary_df(mat_dirac_type2_K18, 16, method_levels)
)

sum_df_dirac_type2$method <- factor(sum_df_dirac_type2$method, levels = method_levels)

df_fdr <- subset(sum_df_dirac_type2, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_dirac_type2, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_dirac_type2 <- plot_fdr(df_fdr) +
  labs(title = "Dirac") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size=30)) 
plot_tpr_dirac_type2 <- plot_tpr(df_tpr)


## Two-point Type 1

sum_df_twopoint_type1 <- rbind(
  make_summary_df(mat_twopoint_type1_K4, 2, method_levels),
  make_summary_df(mat_twopoint_type1_K6, 4, method_levels),
  make_summary_df(mat_twopoint_type1_K10, 8, method_levels),
  make_summary_df(mat_twopoint_type1_K18, 16, method_levels)
)

sum_df_twopoint_type1$method <- factor(sum_df_twopoint_type1$method, levels = method_levels)

df_fdr <- subset(sum_df_twopoint_type1, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_twopoint_type1, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_twopoint_type1 <- plot_fdr(df_fdr) +
  labs(title = "Two-point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size=30)) 
plot_tpr_twopoint_type1 <- plot_tpr(df_tpr)


## Two-point Type 2
sum_df_twopoint_type2 <- rbind(
  make_summary_df(mat_twopoint_type2_K4, 2, method_levels),
  make_summary_df(mat_twopoint_type2_K6, 4, method_levels),
  make_summary_df(mat_twopoint_type2_K10, 8, method_levels),
  make_summary_df(mat_twopoint_type2_K18, 16, method_levels)
)

sum_df_twopoint_type2$method <- factor(sum_df_twopoint_type2$method, levels = method_levels)

df_fdr <- subset(sum_df_twopoint_type2, metric == "FDR")
df_fdr$df_f <- factor(df_fdr$df, levels = c(2,4,8,16))

df_tpr <- subset(sum_df_twopoint_type2, metric == "TPR")
df_tpr$df_f <- factor(df_tpr$df, levels = c(2,4,8,16))

plot_fdr_twopoint_type2 <- plot_fdr(df_fdr) +
  labs(title = "Two-point") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size=30)) 

plot_tpr_twopoint_type2 <- plot_tpr(df_tpr)





#ggsave(filename = "./figure/fdr_invchi_type1.pdf", plot = plot_fdr_invchi_type1, width = 8, height = 7, dpi = 300)
#ggsave(filename = "./figure/tpr_invchi_type1.pdf", plot = plot_tpr_invchi_type1, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/fdr_invchi_type2.pdf", plot = plot_fdr_invchi_type2, width = 8, height = 7, dpi = 300)
#ggsave(filename = "./figure/tpr_invchi_type2.pdf", plot = plot_tpr_invchi_type2, width = 8, height = 6, dpi = 300)

#ggsave(filename = "./figure/fdr_dirac_type1.pdf", plot = plot_fdr_dirac_type1, width = 8, height = 6.8, dpi = 300)
#ggsave(filename = "./figure/tpr_dirac_type1.pdf", plot = plot_tpr_dirac_type1, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/fdr_dirac_type2.pdf", plot = plot_fdr_dirac_type2, width = 8, height = 6.8, dpi = 300)
#ggsave(filename = "./figure/tpr_dirac_type2.pdf", plot = plot_tpr_dirac_type2, width = 8, height = 6, dpi = 300)

#ggsave(filename = "./figure/fdr_twopoint_type1.pdf", plot = plot_fdr_twopoint_type1, width = 8, height = 6.8, dpi = 300)
#ggsave(filename = "./figure/tpr_twopoint_type1.pdf", plot = plot_tpr_twopoint_type1, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/fdr_twopoint_type2.pdf", plot = plot_fdr_twopoint_type2, width = 8, height = 6.8, dpi = 300)
#ggsave(filename = "./figure/tpr_twopoint_type2.pdf", plot = plot_tpr_twopoint_type2, width = 8, height = 6, dpi = 300)
