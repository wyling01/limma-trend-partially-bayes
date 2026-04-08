library(MASS)
library(ggplot2)
library(reshape2)
library(knitr)
library(latex2exp)
library(patchwork)
library(stringr)
library(viridis)
library(RColorBrewer) 
library(tibble)
library(tidyr)
library(dplyr)


func_plot_modified <- function(plot){
  p<- plot + 
    theme_bw(base_size = 14) +  
    theme(
      plot.title = element_text(face = "bold",hjust = 0.5),
      axis.title.x = element_text(face = "bold",size = 22),  
      axis.title.y = element_text(face = "bold",size = 22),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      
      # Remove default axis lines
      axis.line = element_blank(),
      
      # Remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Legend
      legend.text = element_text(face = "bold",size = 14),
      legend.position = "inside",           # Instead of c(x, y)
      legend.position.inside = c(0.98, 0.97), 
      legend.justification = c(1, 1),
      
      # Border around the legend
      #legend.box.background = element_rect(color = "gray50", size = 0.5),
      legend.background = element_rect(fill = "NA"),
      
      # Spacing inside the legend box
      legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),   
      legend.spacing.x = unit(0, "cm"),       
      legend.spacing.y = unit(-0.3, "cm"),  
      legend.key.width = unit(2, "lines"),
      legend.key.height = unit(1.2, "lines"),
      
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20)
    )
  
  p
}


# Trend plot
posterior_weights_2d <- function(a, k, prior_tbl) {
  logw <- log(prior_tbl$prob) +
    dnorm(a,
          mean = prior_tbl$mu,
          sd   = sqrt(prior_tbl$sigma2 / k),
          log  = TRUE)
  m <- max(logw)
  w <- exp(logw - m)
  w / sum(w)
}

posterior_log_sigma2 <- function(a, k, prior_tbl) {
  w <- posterior_weights_2d(a, k, prior_tbl)
  sum(w * log(prior_tbl$sigma2))
}

posterior_logS2 <- function(a,k, d, prior_tbl) {
  posterior_log_sigma2(a, k,prior_tbl) + digamma(d/2) + log(2/d)
}

make_trend_data <- function(info,prior_result){
  idx <- prior_result$joint_npmle$idx
  prior_trend2d_df <- tibble(
    mu     = prior_result$joint_npmle$grid[idx, 1],
    sigma2 = prior_result$joint_npmle$grid[idx, 2],
    prob   = prior_result$joint_npmle$mass[idx]
  )
  E_logS2_given_A <- sapply(
    info$A,
    function(ai) posterior_logS2(ai, info$K, info$df, prior_trend2d_df)
  )
  log_var_fitted_untrended <- (log(prior_result$untrended_invchi$s02*prior_result$untrended_invchi$d0/info$df)
                               +digamma(info$df/2)-digamma(prior_result$untrended_invchi$d0/2))
  result_combined_df <- tibble(A = info$A,
                               log_sample_var = log(info$var),
                               log_var_fitted = info$emean,
                               E_log_S2 = E_logS2_given_A,
                               log_var_fitted_untrended = log_var_fitted_untrended)
  return(result_combined_df)
}

make_trend_plot <- function(trend_data,alpha,xlab){
  plot_trend<- ggplot(trend_data, aes(x = A)) +
    ## raw log-variance points (neutral, no legend)
    geom_point(aes(y = log_sample_var), color = "#AAAAAA", alpha = alpha,shape=16,stroke = 0,size=1) +
    ## fitted trend line
    geom_line(aes(y = log_var_fitted, group = 1, color = "Fitted trend"), size = 1.5,lineend = "round") +
    
    ## posterior expectation points (Joint NPMLE)
    #geom_point(aes(y = E_log_S2), color = "lightblue", size = 0.5,alpha=0.5) +
    geom_line(aes(y = E_log_S2, group = 1, color = "Joint-NPMLE"), size = 1.5,lineend = "round") +
    geom_hline(
      aes(yintercept = log_var_fitted_untrended, color = "Constant trend"),
      size = 1.5,
      linetype = "22"
    ) +
    ylim(min(trend_data$log_sample_var),max(trend_data$log_sample_var)+1)+
    labs(
      x = xlab,
      y = expression(bold(log(S^2))),
      color = NULL
    ) +
    scale_color_manual(
      values = c("Constant trend"="grey40", "Fitted trend" = "#D55E00","Joint-NPMLE" = "#377EB8"),
      breaks = c("Constant trend","Fitted trend","Joint-NPMLE"),
      name = NULL
    ) +
    theme_classic(base_size = 10)
  
  plot_trend
}

# Marginal S
dSg2_marginal <- function(x, d, d0, s0_2, log = FALSE) {
  if (log) {
    # log(f(x/s02)) - log(s02)
    return(df(x / s0_2, df1 = d, df2 = d0, log = TRUE) - log(s0_2))
  } else {
    # f(x/s02) / s02
    return(df(x / s0_2, df1 = d, df2 = d0) / s0_2)
  }
}


marginal_S_g2_npmle <- function(x, df, value, probability) {
  value <- matrix(value, ncol = 1)
  
  temp <- ( (df / 2)^(df / 2) / gamma(df / 2) * 
              value^(-df / 2) * 
              x^(df / 2 - 1) * 
              exp(-df * x / (2 * value)) )
  result <- colSums(probability * temp)
  return(result)
}

make_marginal_S_data <- function(info,prior_result,include_joint = TRUE,length.out=3000){
  x_vals <- exp(seq(log(min(info$var)), log(max(info$var)), length.out = length.out))
  y_vals_untrended_invchi <- sapply(x_vals, dSg2_marginal, d = info$df, d0 = prior_result$untrended_invchi$d0, s0_2 = prior_result$untrended_invchi$s02)
  y_vals_untrended_npmle <- sapply(x_vals, marginal_S_g2_npmle , df = info$df, value = prior_result$untrended_npmle$grid, probability = prior_result$untrended_npmle$mass)
  
  marginal_S_data <- tibble::tibble(
    x = x_vals,
    y_untrended_npmle = y_vals_untrended_npmle,
    y_untrended_invchi = y_vals_untrended_invchi
  )

  if (include_joint) {
    idx <- prior_result$joint_npmle$idx
    y_vals_joint_npmle <- sapply(
      x_vals, marginal_S_g2_npmle,
      df = info$df,
      value = prior_result$joint_npmle$grid[idx, 2],
      probability = prior_result$joint_npmle$mass[idx]
    )
    marginal_S_data$y_joint_npmle <- y_vals_joint_npmle
  }

  return(marginal_S_data)
}

make_marginal_S_plot<- function(marginal_S_data,info,xlim_R = 1,include_joint = TRUE){
  p <- ggplot() +
    geom_histogram(
      data = data.frame(S = info$var),
      aes(x = S, y = after_stat(density), fill = "Histogram"),
      color = "grey50", bins = 100, alpha = 0.3
    ) +
    geom_line(data = marginal_S_data, aes(x = x, y = y_untrended_npmle, color = "Untrended-NPMLE"), linewidth = 1) +
    geom_line(data = marginal_S_data, aes(x = x, y = y_untrended_invchi, color = "Untrended-Inv"), linewidth = 1) +
    xlim(0, xlim_R) +
    labs(x = expression(bold(S[i]^2)), y = "Density") +
    scale_fill_manual(values = c("Histogram" = "grey"), name = NULL) +
    theme_minimal()

  # Optionally add joint curve
  if (include_joint) {
    p <- p + geom_line(data = marginal_S_data, aes(x = x, y = y_joint_npmle, color = "Joint-NPMLE"), linewidth = 1)
  }

  # Build color scale conditionally so legend matches what is plotted
  color_values <- c("Untrended-Inv" = "lightsalmon1",
                    "Untrended-NPMLE" = "aquamarine3")
  color_breaks <- c("Untrended-Inv", "Untrended-NPMLE")
  color_labels <- c(expression(bold("Untrended-Inv"*chi^2)),
                    "Untrended-NPMLE")

  if (include_joint) {
    color_values <- c(color_values, "Joint-NPMLE" = "#377EB8")
    color_breaks <- c(color_breaks, "Joint-NPMLE")
    color_labels <- c(color_labels, "Joint-NPMLE")
  }
  p <- p +
    scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      labels = color_labels,
      name = NULL
    ) +
    guides(
      fill  = guide_legend(order = 1),
      color = guide_legend(order = 2)
    )

  p

}

# Prior sigma
d_scaled_invchisq <- function(t, d0, s02) {
  cst <- ( (d0*s02/2)^(d0/2) ) / gamma(d0/2)
  cst * t^(-(d0/2 + 1)) * exp(-(d0*s02)/(2*t))
}

make_prior_sigma_plot <- function(prior_result,scale_factor,length.out = 5000,xlim_R =1){
  df_npmle <- data.frame(value = prior_result$untrended_npmle$grid, probability = prior_result$untrended_npmle$mass*scale_factor)
  xgrid <- seq(min(prior_result$untrended_npmle$grid[prior_result$untrended_npmle$grid > 0]),
               max(prior_result$untrended_npmle$grid),
               length.out = length.out)
  df_invchi <- data.frame(
    value = xgrid,
    prob_curve = d_scaled_invchisq(xgrid, d0 =prior_result$untrended_invchi$d0, s02 = prior_result$untrended_invchi$s02)
  )
  
  p <-  ggplot(df_npmle, aes(x = value, y = probability)) +
    geom_segment(aes(xend = value, yend = 0, color = "Untrended-NPMLE"),linewidth = 1) +
    geom_line(data = df_invchi,aes(x = value, y = prob_curve, color = "Untrended-Inv"),linewidth = 1) +
    labs(
      x = expression(bold(sigma^2)),
      y = expression(bold(g(sigma^2))),
      #title = "Estimated prior with NPMLE_trend_1d",
      color = NULL
    ) +
    scale_color_manual(
      values = c("Untrended-Inv" = "lightsalmon1","Untrended-NPMLE" = "aquamarine3"),
      breaks = c( "Untrended-Inv","Untrended-NPMLE"),
      labels = c(expression(bold("Untrended-Inv"*chi^2)),"Untrended-NPMLE")
    ) +
    xlim(0,xlim_R)+
    theme_minimal()
  p
}

# Marginal V
make_marginal_V_data <- function(info,prior_result,length.out=3000){
  x_vals <- exp(seq(log(min(info$V)), log(max(info$V)), length.out = length.out))
  y_vals_reg_invchi <- sapply(x_vals, dSg2_marginal, d = info$df, d0 = prior_result$reg_invchi$d0, s0_2 = prior_result$reg_invchi$v02)
  y_vals_reg_npmle <- sapply(x_vals, marginal_S_g2_npmle , df = info$df, value = prior_result$reg_npmle$grid, probability = prior_result$reg_npmle$mass)
  
  marginal_V_data <- tibble(
    x = x_vals,
    y_reg_npmle = y_vals_reg_npmle,
    y_reg_invchi = y_vals_reg_invchi
  )
  return(marginal_V_data)
}

make_marginal_V_plot<- function(marginal_V_data,info,xlim_R = 1){
  p <- ggplot() +
    geom_histogram(data = data.frame(V= info$V), 
                   aes(x = V, y = after_stat(density), fill = "Histogram"),color = "grey50",
                   bins = 100, alpha = 0.3) +
    geom_line(data = marginal_V_data, aes(x = x, y = y_reg_npmle , color = "Reg-NPMLE"), linewidth = 1) +
    geom_line(data = marginal_V_data, aes(x = x, y = y_reg_invchi , color = "Reg-Inv"), linewidth = 1) +
    xlim(0, xlim_R) +
    labs( x = expression(bold(V[i]^2)), y = "Density") +
    scale_fill_manual(values = c("Histogram" = "grey"), name = NULL) +  # Legend for histogram
    scale_color_manual(
      values = c( "Reg-Inv" = "mediumpurple4","Reg-NPMLE" = "lightcoral"),
      breaks = c( "Reg-Inv","Reg-NPMLE"),
      labels = c( expression(bold("Reg-Inv"*chi^2)),"Reg-NPMLE"),
      name = NULL
    ) +
    theme_minimal() 
  p
}

# Prior tau
make_prior_tau_plot <- function(prior_result,scale_factor,length.out = 5000,xlim_R =1){
  df_npmle <- data.frame(value = prior_result$reg_npmle$grid, probability = prior_result$reg_npmle$mass*scale_factor)
  xgrid <- seq(min(prior_result$reg_npmle$grid[prior_result$reg_npmle$grid > 0]),
               max(prior_result$reg_npmle$grid),
               length.out = length.out)
  df_invchi <- data.frame(
    value = xgrid,
    prob_curve = d_scaled_invchisq(xgrid, d0 =prior_result$reg_invchi$d0, s02 = prior_result$reg_invchi$v02)
  )
  
  p <-  ggplot(df_npmle, aes(x = value, y = probability)) +
    geom_segment(aes(xend = value, yend = 0, color = "Reg-NPMLE"),linewidth = 1) +
    geom_line(data = df_invchi,aes(x = value, y = prob_curve, color = "Reg-Inv"),linewidth = 1) +
    labs(
      x = expression(bold(tau^2)),
      y = expression(bold(g(tau^2))),
      #title = "Estimated prior with NPMLE_trend_1d",
      color = NULL
    ) +
    scale_color_manual(
      values = c( "Reg-Inv" = "mediumpurple4","Reg-NPMLE" = "lightcoral"),
      breaks = c( "Reg-Inv","Reg-NPMLE"),
      labels = c( expression(bold("Reg-Inv"*chi^2)),"Reg-NPMLE")
    ) +
    xlim(0,xlim_R)+
    theme_minimal()
  p
}

# Joint Prior
make_joint_prior_plot <- function(prior_result){
  idx_plot <- which(prior_result$joint_npmle$mass>=1e-5)
  
  prob_joint <- prior_result$joint_npmle$mass[idx_plot]
  mass_joint <- prior_result$joint_npmle$grid[idx_plot,]
  
  sorted_indices <- order(prob_joint)  # Sort indices based on probability
  sorted_values <- mass_joint[sorted_indices, ]  # Sort values accordingly
  sorted_probabilities <- prob_joint[sorted_indices]  # Sorted probabilities
  
  df <- data.frame(
    mu = sorted_values[, 1],
    sigma_2 = sorted_values[, 2],
    Probability = sorted_probabilities
  )
  
  p <- ggplot(df, aes(x = mu, y = sigma_2, color = Probability)) +
    geom_point(size = 2, alpha = 1) +  
    scale_color_distiller(palette = "Blues", direction = 1) +
    scale_y_log10() +
    labs(
      x = expression(bold(mu)),
      y = expression(bold(sigma^2)),
      color = "Probability"
    ) +
    theme_minimal()
  p
}

# Trend-Mimic-Joint
E_logS2_given_Mi <- function(pep_grp, prior_byM, df1) {
  pep_grp_chr <- as.character(pep_grp)

  E_log_sigma2_bygrp <- sapply(names(prior_byM), function(g) {
    df_g <- prior_byM[[g]]
    w <- df_g$prob / sum(df_g$prob)
    sum(w * log(df_g$value))
  })

  E_log_sigma2_i <- E_log_sigma2_bygrp[pep_grp_chr]

  E_logS2_i <- E_log_sigma2_i + digamma(df1 / 2) + log(2 / df1)

  return(E_logS2_i)
}

make_trend_plot_mimic_joint <- function(info,prior_result,pep_grp,alpha,xlab){
  df<- info$df
  E_logS2_Mi <- E_logS2_given_Mi(pep_grp, prior_result$mimic_joint_npmle, df)
  log_var_fitted_untrended <- (log(prior_result$untrended_invchi$s02*prior_result$untrended_invchi$d0/info$df)
                               +digamma(info$df/2)-digamma(prior_result$untrended_invchi$d0/2))
  result_combined_df <- tibble(A = info$A,
                    log_sample_var = log(info$var),
                    log_var_fitted = info$emean,
                    E_log_S2 = E_logS2_Mi,
                    log_var_fitted_untrended = log_var_fitted_untrended)
  
  
  plot_trend <- ggplot(result_combined_df, aes(x = A)) +
    ## raw log-variance points (neutral, no legend)
    geom_point(aes(y = log_sample_var), color = "#AAAAAA", alpha = alpha,shape=16,stroke = 0,size=1) +
    ## fitted trend line
    geom_line(aes(y = log_var_fitted, group = 1, color = "Fitted trend"), size = 1.5,lineend = "round") +
  
    ## posterior expectation points (Joint NPMLE)
    #geom_point(aes(y = E_S2), color = "lightblue",size = 0.5,alpha=0.5) +
    geom_line(aes(y = E_log_S2, group = 1, color = "Joint-NPMLE"), size = 1.5,lineend = "round") +
    geom_hline(
      aes(yintercept = log_var_fitted_untrended, color = "Constant trend"),
      size = 1.5,
      linetype = "22"
    ) +
    labs(
       x = xlab,
       y = expression(bold(log(S^2))),
      color = NULL
    ) +
    scale_color_manual(
      values = c( "Constant trend"="grey40", "Fitted trend" = "#D55E00","Joint-NPMLE" = "#377EB8"),
      breaks = c( "Constant trend","Fitted trend","Joint-NPMLE"),
      name = NULL
    ) +
     theme_classic(base_size = 10)
  
  plot_trend
}

# Significance plot
summary_significance_pep_grp <- function(result,pep_grp,aloha=0.05){
  df_sig <- tibble(
    M = pep_grp,
    t          = BH_adjust(result$t_test,alpha),
    untrend    = BH_adjust(result$untrended_invchi,alpha),
    untrend_1d = BH_adjust(result$untrended_npmle,alpha),
    trend      = BH_adjust(result$reg_invchi,alpha),
    trend_1d   = BH_adjust(result$reg_npmle,alpha),
    trend_2d   = BH_adjust(result$mimic_joint_npmle,alpha)
  )
  method_cols <- setdiff(names(df_sig), "M")
  out <- df_sig %>%
    group_by(M) %>%
    summarise(
      n = dplyr::n(),
      across(all_of(method_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(M)

  df_long <- out %>%
    pivot_longer(
      cols = -c(M, n),
      names_to = "method",
      values_to = "sig"
    ) %>%
    mutate(
      M = factor(M, levels = sort(unique(M))),  # equal spacing on x
      method = factor(method,
        levels = c("t","untrend","untrend_1d","trend","trend_1d","trend_2d")
      )
    )
  df_prop_long <- out %>%
    pivot_longer(
      cols = all_of(method_cols),
      names_to = "method",
      values_to = "sig"
    ) %>%
    mutate(
      prop = sig / n,
      method = factor(method, levels = c("t","untrend","untrend_1d","trend","trend_1d","trend_2d"))
    )
  
  return(list(df_sum=df_long,df_prop=df_prop_long))
}


