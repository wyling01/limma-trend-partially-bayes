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
library(grid)
library(scales)


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

posterior_logS2 <- function(a, k, d, prior_tbl) {
  posterior_log_sigma2(a, k, prior_tbl) + digamma(d / 2) + log(2 / d)
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

# Log Marginal
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

make_logvar_marginal_data <- function(info,prior_result,length.out = 3000,
                                      xlim_log = c(-10, 10)) {
  
  eps <- 1e-12
  method_levels <- c("Untrended-Inv","Untrended-NPMLE","Reg-Inv","Reg-NPMLE")
  
  hist_df <- tibble(logS2 = log(pmax(info$var, eps)),logV2 = log(pmax(info$V, eps)))
  
  hist_df_long <- hist_df %>%
    pivot_longer(
      cols = c(logS2, logV2),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(
      type = factor(type, levels = c("logS2", "logV2"))
    )
  
  t_vals <- seq(xlim_log[1], xlim_log[2], length.out = length.out)
  x_vals <- exp(t_vals)
  
  y_untrended_invchi <- sapply(
    x_vals,
    dSg2_marginal,
    d = info$df,
    d0 = prior_result$untrended_invchi$d0,
    s0_2 = prior_result$untrended_invchi$s02
  )
  
  y_untrended_npmle <- sapply(
    x_vals,
    marginal_S_g2_npmle,
    df = info$df,
    value = prior_result$untrended_npmle$grid,
    probability = prior_result$untrended_npmle$mass
  )
  
  y_reg_invchi <- sapply(
    x_vals,
    dSg2_marginal,
    d = info$df,
    d0 = prior_result$reg_invchi$d0,
    s0_2 = prior_result$reg_invchi$v02
  )
  
  y_reg_npmle <- sapply(
    x_vals,
    marginal_S_g2_npmle,
    df = info$df,
    value = prior_result$reg_npmle$grid,
    probability = prior_result$reg_npmle$mass
  )
  
  line_df <- bind_rows(
    tibble(x = t_vals,y = x_vals * y_untrended_invchi,method = "Untrended-Inv"),
    tibble(x = t_vals,y = x_vals * y_untrended_npmle,method = "Untrended-NPMLE"),
    tibble(x = t_vals,y = x_vals * y_reg_invchi,method = "Reg-Inv"),
    tibble(x = t_vals,y = x_vals * y_reg_npmle,method = "Reg-NPMLE")
  ) %>%
    filter(is.finite(x), is.finite(y), y >= 0) %>%
    mutate(method = factor(method, levels = method_levels))
  
  list(
    hist_df_long = hist_df_long,
    line_df = line_df,
    xlim_log = xlim_log,
    method_levels = method_levels
  )
}

plot_logvar_marginal_data <- function(plot_data,bins = 60,line_size = 1.5){
  
  hist_df_long <- plot_data$hist_df_long
  line_df <- plot_data$line_df
  xlim_log <- plot_data$xlim_log
  method_levels <- plot_data$method_levels
  
  color_values <- c("Untrended-Inv" = "mediumpurple4","Untrended-NPMLE" = "lightcoral",
                    "Reg-Inv" = "aquamarine4","Reg-NPMLE" = "lightskyblue2")
  
  linetype_values <- c("Untrended-Inv" = "dashed","Untrended-NPMLE" = "solid",
                       "Reg-Inv" = "dashed","Reg-NPMLE" = "solid")
  
  legend_labels <- c(expression(bold("Untrended-Inv" * chi^2)),"Untrended-NPMLE",
                     expression(bold("Reg-Inv" * chi^2)),"Reg-NPMLE")
  
  p <- ggplot() +
    geom_histogram(data = hist_df_long,aes(x = value,y = after_stat(density),fill = type),
                   bins = bins,alpha = 0.35,position = "identity",color = "grey50",linewidth = 0.25) +
    
    geom_line(data = line_df %>%
                filter(method %in% c("Untrended-NPMLE", "Reg-NPMLE")),
              aes(x = x,y = y,color = method,linetype = method),
              linewidth = line_size,lineend = "round") +
    
    geom_line(data = line_df %>%
                filter(method %in% c("Untrended-Inv", "Reg-Inv")),
              aes(x = x,y = y,color = method,linetype = method),
              linewidth = line_size,lineend = "round") +
    
    scale_fill_manual(values = c(logS2 = "#D55E00",logV2 = "grey50"),
                      breaks = c("logS2", "logV2"),
                      labels = c(expression(bold(log(S[i]^2))),
                                 expression(bold(log(V[i]^2)))),
                      name = NULL) +
    
    scale_color_manual(values = color_values,breaks = method_levels,
                       labels = legend_labels,name = NULL) +
    scale_linetype_manual(values = linetype_values,breaks = method_levels,
                          labels = legend_labels,name = NULL) +
    coord_cartesian(xlim = xlim_log) +
    
    labs(x = expression(bold(log(S[i]^2)~"or"~log(V[i]^2))),y = "Density") +
    
    guides(
      fill = guide_legend(order = 1,ncol = 2,
                          override.aes = list(alpha = 0.35,color = "grey50",linewidth = 0.3),
                          keywidth = unit(1.8, "lines"),keyheight = unit(0.2, "lines")),
      color = guide_legend(order = 2,
                           override.aes = list(linewidth = line_size,linetype = c("dashed", "solid", "dashed", "solid")),
                           keywidth = unit(3, "lines"),keyheight = unit(0.8, "lines")),
      linetype = "none"
    )
  
  return(p)
}

# Log Prior
d_scaled_invchisq <- function(t, d0, s02) {
  cst <- ( (d0*s02/2)^(d0/2) ) / gamma(d0/2)
  cst * t^(-(d0/2 + 1)) * exp(-(d0*s02)/(2*t))
}

make_npmle_prior_log_data <- function(prior_result,threshold = 1e-8,length.out = 5000,
                                      xlim_log = c(-5, 4),npmle_scale_factor = NULL) {
  
  eps <- 1e-12
  method_levels <- c("Untrended-Inv","Untrended-NPMLE","Reg-Inv","Reg-NPMLE")
  
  untrended_npmle_df <- tibble(
    x = log(pmax(prior_result$untrended_npmle$grid, eps)),
    mass = prior_result$untrended_npmle$mass,
    method = "Untrended-NPMLE"
  )
  
  reg_npmle_df <- tibble(
    x = log(pmax(prior_result$reg_npmle$grid, eps)),
    mass = prior_result$reg_npmle$mass,
    method = "Reg-NPMLE"
  )
  
  npmle_df <- bind_rows(untrended_npmle_df, reg_npmle_df) %>%
    filter(mass > threshold,is.finite(x),x >= xlim_log[1],x <= xlim_log[2])
  
  t_vals <- seq(xlim_log[1], xlim_log[2], length.out = length.out)
  u_vals <- exp(t_vals)
  
  untrended_inv_density <- d_scaled_invchisq(u_vals,d0 = prior_result$untrended_invchi$d0,
                                             s02 = prior_result$untrended_invchi$s02)
  
  reg_inv_density <- d_scaled_invchisq(u_vals,d0 = prior_result$reg_invchi$d0,
                                       s02 = prior_result$reg_invchi$v02)
  
  inv_df <- bind_rows(
    tibble(x = t_vals,y = untrended_inv_density * u_vals,method = "Untrended-Inv"),
    tibble(x = t_vals,y = reg_inv_density * u_vals,method = "Reg-Inv")
  ) %>%
    filter(is.finite(x), is.finite(y), y >= 0)
  
  if (is.null(npmle_scale_factor)) {
    npmle_scale_factor <- max(inv_df$y, na.rm = TRUE) / 
      max(npmle_df$mass, na.rm = TRUE)
  }
  
  npmle_df <- npmle_df %>%
    mutate(y = mass * npmle_scale_factor,method = factor(method, levels = method_levels))
  
  inv_df <- inv_df %>%
    mutate(method = factor(method, levels = method_levels))
  
  list(
    npmle_df = npmle_df,
    inv_df = inv_df,
    npmle_scale_factor = npmle_scale_factor,
    xlim_log = xlim_log,
    method_levels = method_levels
  )
}

plot_npmle_prior_log_data <- function(plot_data,line_size = 1.5) {
  
  npmle_df <- plot_data$npmle_df
  inv_df <- plot_data$inv_df
  xlim_log <- plot_data$xlim_log
  method_levels <- plot_data$method_levels
  
  legend_labels <- c(expression(bold("Untrended-Inv" * chi^2)),"Untrended-NPMLE",
                     expression(bold("Reg-Inv" * chi^2)),"Reg-NPMLE"
  )
  
  color_values <- c("Untrended-Inv" = "mediumpurple4","Untrended-NPMLE" = "lightcoral", 
                    "Reg-Inv" = "aquamarine4","Reg-NPMLE" = "lightskyblue2")
  
  linetype_values <- c("Untrended-Inv" = "dashed","Untrended-NPMLE" = "solid",
                       "Reg-Inv" = "dashed","Reg-NPMLE" = "solid")
  
  # Dummy data only for legend alignment
  legend_x0 <- xlim_log[1] - diff(xlim_log) * 2
  
  legend_df <- tibble(
    x = rep(c(legend_x0, legend_x0 + 0.01), length(method_levels)),
    y = rep(0, 2 * length(method_levels)),
    method = factor(rep(method_levels, each = 2),levels = method_levels)
  )
  
  p <- ggplot() +
    geom_segment(data = npmle_df,aes(x = x,xend = x,y = 0,yend = y,color = method,linetype = method),
                 linewidth = line_size,lineend = "butt",show.legend = FALSE) +
    geom_line(data = inv_df,aes(x = x,y = y,color = method,linetype = method),
              linewidth = line_size,lineend = "round",show.legend = FALSE) +
    geom_line(data = legend_df,aes(x = x,y = y,color = method,linetype = method),
              linewidth = line_size,lineend = "round",show.legend = TRUE) +
    
    scale_color_manual(values = color_values,breaks = method_levels,
                       labels = legend_labels,name = NULL) +
    
    scale_linetype_manual(values = linetype_values,breaks = method_levels,
                          labels = legend_labels,name = NULL) +
    
    coord_cartesian(xlim = xlim_log) +
    
    labs(x = expression(bold(log(sigma^2)~"or"~log(tau^2))),
         y = "Prior")+
    
    guides(
      color = guide_legend(order = 1,override.aes = list(linewidth = line_size,
                                                         linetype = c("dashed", "solid", "dashed", "solid")),
                           keywidth = unit(3, "lines"),keyheight = unit(0.8, "lines")
      ),
      linetype = "none"
    ) +
    
    theme_minimal()
  
  return(p)
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

make_trend_plot_mimic_joint <- function(info,mimic_joint_npmle_prior,pep_grp){
  df<- info$df
  E_logS2_Mi <- E_logS2_given_Mi(pep_grp, mimic_joint_npmle_prior, df)
  
  result_combined_df <- tibble(A = info$A,
                    var = log(info$var),
                    var_fitted = info$emean,
                    E_S2 = E_logS2_Mi)
  
  
  plot_trend <- ggplot(result_combined_df, aes(x = A)) +
    ## raw log-variance points (neutral, no legend)
    geom_point(aes(y = var), color = "grey70", alpha = 0.6,shape=16,stroke = 0) +
    ## fitted trend line
    geom_line(aes(y = var_fitted, group = 1, color = "Fitted trend"), linewidth = 1.2) +
  
    ## posterior expectation points (Joint NPMLE)
    geom_point(aes(y = E_S2), color = "lightblue",size = 0.5,alpha=0.5) +
    geom_line(aes(y = E_S2, group = 1, color = "Joint-NPMLE"), linewidth = 1.2) +
    labs(
       x = expression(bold(log[2](count))),
       y = expression(bold(log(S^2))),
      color = NULL
    ) +
    scale_color_manual(
      values = c( "Fitted trend" = "#FB8072","Joint-NPMLE" = "#377EB8"),
      breaks = c( "Fitted trend","Joint-NPMLE"),
      name = NULL
    ) +
    theme_minimal()
  
  plot_trend
}

# Significance plot
summary_significance_pep_grp <- function(result,pep_grp,alpha=0.05){
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


