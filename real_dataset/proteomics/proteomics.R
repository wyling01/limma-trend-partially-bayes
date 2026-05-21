library(oligo)
library(edgeR)
library(readxl)
library(patchwork)
library(gtable)
library(grid)
library(scales)

# -----Load Function-----
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
# -----Data Preprocessing-----

## Download data
url <- "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
download.file(url, destfile = "./miR_Proteintable.txt",method = "auto")

df.prot = read.table("miR_Proteintable.txt",stringsAsFactors = FALSE,
                     header = TRUE, quote = "", comment.char = "",sep = "\t")

## Extract quant data columns for DEqMS

# filter at 1% protein FDR and extract TMT quantifications
TMT_columns = seq(15,33,2)
dat = df.prot[df.prot$miR.FASP_q.value<0.01,TMT_columns]
df.prot_new <- df.prot[df.prot$miR.FASP_q.value < 0.01, ]
rownames(df.prot_new) <- df.prot_new$Protein.accession
rownames(dat) = df.prot[df.prot$miR.FASP_q.value<0.01,]$Protein.accession
dat.log = log2(dat)
#remove rows with NAs
dat.log = na.omit(dat.log)
df.prot_new <- df.prot_new[rownames(dat.log), ]
boxplot(dat.log,las=2,main="TMT10plex data PXD004163")

# Make design table
cond = as.factor(c("ctrl","miR191","miR372","miR519","ctrl",
                   "miR372","miR519","ctrl","miR191","miR372"))
design = model.matrix(~0+cond) 
colnames(design) = gsub("cond","",colnames(design))

# Make contrasts
x <- c("miR372-ctrl","miR519-ctrl","miR191-ctrl",
       "miR372-miR519","miR372-miR191","miR519-miR191")
contrast =  makeContrasts(contrasts=x,levels=design)

# Use number of peptides as $M_i$
M <- log2(df.prot_new[, 8] + 1)

### Check if there is NA in the expression matrix. If yes, remove rows that contain NAs.
any_na <- anyNA(dat.log)

if (any_na) {
  keep <- complete.cases(dat.log)  
  dat.log <- dat.log[keep, , drop = FALSE]
}

pep_num <- df.prot_new[, 8]
## Save data
data_proteomics <- list(Y=dat.log,design = design,contrast=contrast,M=M,pep_num = pep_num)

save(data_proteomics, file = "./data/filtered_data_proteomics.RData")


# -----DE analysis-----
load("./data/filtered_data_proteomics.RData")



info <- info_extractor(data_proteomics,M=data_proteomics$M)
info$pep_num <- data_proteomics$pep_num


## Prior estimation
param_prior <- var_invchi_prior(info)
untrended_invchi_prior <- param_prior[[1]]
reg_invchi_prior <- param_prior[[2]]

untrended_npmle_prior <- prior_untrended_npmle(info,v=300)
reg_npmle_prior <- prior_reg_npmle(info,v=300)

prior_result <- list(
  untrended_invchi=untrended_invchi_prior,
  untrended_npmle=untrended_npmle_prior,
  reg_invchi =reg_invchi_prior,
  reg_npmle=reg_npmle_prior
)

## P value calculation (change contreast_name from c(1,2,3,4,5,6))
result <- p_value_calculator(info,prior_result,contrast_name=1,method_list=c(1,2,3,4,5),verbose=TRUE)

alpha = 0.05
sum(BH_adjust(result$t_test,alpha))
sum(BH_adjust(result$untrended_invchi,alpha))
sum(BH_adjust(result$untrended_npmle,alpha))
sum(BH_adjust(result$reg_invchi,alpha))
sum(BH_adjust(result$reg_npmle,alpha))

# -----Plot-----

## Log Marginal
logvar_data <- make_logvar_marginal_data(
  info = info,
  prior_result = prior_result,
  length.out = 3000,
  xlim_log = c(-10, 9)
)

plot_marginal_combined <-  plot_logvar_marginal_data(plot_data = logvar_data,bins = 100)
plot_marginal_combined <- func_plot_modified(plot_marginal_combined ) +ylim(0,0.55)

## Log Prior
prior_log_data <- make_npmle_prior_log_data(
  prior_result = prior_result,
  threshold = 1e-8,
  length.out = 5000,
  xlim_log = c(-8, 8)
)
plot_prior_combined  <- plot_npmle_prior_log_data(plot_data = prior_log_data)
plot_prior_combined <- func_plot_modified(plot_prior_combined) 



#ggsave(filename = "./figure/log_marginal_combined.pdf", plot = plot_marginal_combined , width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/log_prior_combined.pdf", plot = plot_prior_combined, width = 8, height = 6, dpi = 300)


# -----Mimic-joint-npmle-----
# self-defined group
pep <- data_proteomics$pep_num  # raw peptide count M_i
n   <- length(pep)

## ---- define bins ----
pep_grp <- as.character(pep)

# keep exact for 1...11
pep_grp[pep >= 12 & pep <= 13]  <- "12-13"
pep_grp[pep >= 14 & pep <= 16]  <- "14-16"
pep_grp[pep >= 17 & pep <= 20]  <- "17-20"
pep_grp[pep >= 21 & pep <= 25]  <- "21-25"
pep_grp[pep >= 26 & pep <= 30]  <- "26-30"
pep_grp[pep >= 31 & pep <= 45]  <- "31-45"
pep_grp[pep >= 46]  <- "46+"


pep_grp <- factor(pep_grp, levels = c(
  as.character(1:11),
  "12-13","14-16","17-20","21-25","26-30","31-45","46+"
))

table(pep_grp)

# Estimating prior
mimic_joint_npmle_prior <- prior_mimic_joint_npmle(info,pep_grp,v=80,verbose = FALSE)

# P value calculation
P_list_mimic_joint_npmle <- P_value_mimic_joint_npmle(info,mimic_joint_npmle_prior,pep_grp,contrast_name=1)
result$mimic_joint_npmle <- P_list_mimic_joint_npmle
sum(BH_adjust(P_list_mimic_joint_npmle,alpha))

# -----Mimic-joint-npmle plot-----
## trend plot

plot_trend_mimic_joint <- make_trend_plot_mimic_joint(info,mimic_joint_npmle_prior,pep_grp)
plot_trend_mimic_joint <- func_plot_modified(plot_trend_mimic_joint)

#ggsave(filename = "./figure/combined_trend.pdf", plot = plot_trend_mimic_joint, width = 20, height = 6, dpi = 300)


## Prop-significance plot
out <- summary_significance_pep_grp(result,pep_grp)
method_breaks <- c("t", "untrend", "untrend_1d", "trend", "trend_1d", "trend_2d")

method_labels <- c("t-test",expression(bold("Untrended-Inv" * chi^2)),"Untrended-NPMLE",
                   expression(bold("Reg-Inv" * chi^2)),"Reg-NPMLE","Joint-NPMLE")

method_colors <- c("t" = "gray40","untrend" = "mediumpurple4","untrend_1d" = "lightcoral",
                   "trend" = "aquamarine4","trend_1d" = "lightskyblue2","trend_2d" = "lightsalmon1")

method_shapes <- c("t" = 17,"untrend" = 15,"untrend_1d" = 3,"trend" = 18,"trend_1d" = 7,"trend_2d" = 8)

p_prop <- ggplot(out$df_prop,aes(x = M, y = prop, color = method, group = method, shape = method)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1),expand = expansion(mult = c(0, 0.3))) +
  labs(x = expression(bold(M[i])),y = "Proportion of significance",color = NULL,shape = NULL) +
  scale_color_manual(values = method_colors,breaks = method_breaks,labels = method_labels,name = NULL) +
  scale_shape_manual(values = method_shapes,breaks = method_breaks,labels = method_labels,name = NULL) +
  theme_minimal(base_size = 14)

p_prop_final<- func_plot_modified(p_prop)+
  guides(
    color = guide_legend(
      nrow = 2,byrow = TRUE,
      override.aes = list(shape = unname(method_shapes[method_breaks]),linewidth = 1.5,size = 3),
      keywidth = unit(1.6, "lines"),keyheight = unit(0.9, "lines")
    ),
    shape = "none"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.06, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.25, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.key.width = unit(1.6, "lines"),
    legend.key.height = unit(0.9, "lines"),
    axis.text.x = element_text(size = 15, angle = 40, hjust = 1, vjust = 1)
  )
#ggsave(filename = "./figure/prop_discoveries.pdf", plot = p_prop_final, width = 8, height = 6, dpi = 300)

## Mimic-joint-prior plot
make_horizontal_joint_prior_plot <- function(
  prior_result,
  ylim_sigma = c(5e-4, 1e-1),
  max_width = 0.75,
  grid_width = 0.85,
  line_size = 1,
  color = "aquamarine3"
) {
  
  df_spike <- bind_rows(lapply(names(prior_result$mimic_joint_npmle), function(g) {
    d <- prior_result$mimic_joint_npmle[[g]]
    data.frame(pep_grp = g,value = d$value,prob = d$prob)
  }))
  
  bins <- unique(df_spike$pep_grp)
  
  df_plot <- df_spike %>%
    mutate(pep_grp = factor(pep_grp, levels = bins),bin_id = as.numeric(pep_grp))
  
  max_prob <- max(df_plot$prob, na.rm = TRUE)
  
  df_plot <- df_plot %>%
    mutate(x_start = bin_id - max_width / 2,x_end = x_start + max_width * prob / max_prob)
  
  major_breaks <- c(1e-3, 1e-2, 1e-1)
  minor_breaks <- c(3e-3, 3e-2)
  
  grid_major_df <- expand.grid(
    bin_id = seq_along(bins),
    y = major_breaks
  ) %>%
    mutate(x_start = bin_id - grid_width / 2,x_end = bin_id + grid_width / 2)
  
  grid_minor_df <- expand.grid(
    bin_id = seq_along(bins),
    y = minor_breaks
  ) %>%
    mutate(x_start = bin_id - grid_width / 2,x_end = bin_id + grid_width / 2)
  
  p <- ggplot(df_plot) +
    
    geom_segment(
      data = grid_minor_df,
      aes(x = x_start,xend = x_end,y = y,yend = y),
      inherit.aes = FALSE,color = "grey93",linewidth = 0.35
    ) +
    
    geom_segment(
      data = grid_major_df,
      aes(x = x_start,xend = x_end,y = y,yend = y),
      inherit.aes = FALSE,
      color = "grey88",
      linewidth = 0.55
    ) +
    
    geom_segment(
      aes(x = x_start,xend = x_end,y = value,yend = value),
      linewidth = line_size,
      color = color,
      lineend = "butt"
    ) +
    
    scale_x_continuous(breaks = seq_along(bins),labels = bins,expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_log10(breaks = major_breaks,minor_breaks = NULL,labels = c("0.001", "0.010", "0.100")) +
    coord_cartesian(ylim = ylim_sigma) +
    labs(x = expression(bold(M[i])),y = expression(bold(sigma^2))
    ) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text.x  = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = element_text(size = 12),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  p
}
p_horizontal <- make_horizontal_joint_prior_plot(
  prior_result = prior_result,
  ylim_sigma = c(5e-4, 1e-1),
  max_width = 0.75,
  grid_width = 0.85,
  line_size = 0.7,
  color = "aquamarine3"
)
#ggsave(filename = "./figure/mimic-joint-prior-flip.pdf", plot = p_horizontal, width = 8, height = 6, dpi = 300)

## Log marginal with different values of M
marginal_S2_reg_npmle_given_A <- function(x, A, df, value, probability, mhat_fun) {
  xi2 <- exp(mhat_fun(A))
  (1 / xi2) * marginal_S_g2_npmle(
    x / xi2,
    df = df,
    value = value,
    probability = probability
  )
}

make_marginal_S2_data_by_A <- function(info_sub, prior_result, group_id, mhat_fun,
                                  length.out = 3000) {
  x_vals <- exp(seq(log(min(info_sub$var)), log(max(info_sub$var)), length.out = length.out))
  prior_group <- prior_result$mimic_joint_npmle[[group_id]]

  y_vals_joint_npmle_discrete <- sapply(x_vals, marginal_S_g2_npmle , df = info_sub$df, 
                                        value = prior_group$value, probability = prior_group$prob)

  y_vals_reg_npmle <- sapply(
    x_vals,
    marginal_S2_reg_npmle_given_A,
    A = info_sub$A[1],
    df = info_sub$df,
    value = prior_result$reg_npmle$grid,
    probability = prior_result$reg_npmle$mass,
    mhat_fun = mhat_fun
  )

  tibble::tibble(
    x = x_vals,
    y_joint_npmle_discrete = y_vals_joint_npmle_discrete,
    y_reg_npmle = y_vals_reg_npmle
  )
}

make_marginal_S_by_A_plot<- function(marginal_S_data_by_A,info,xlim_R = 1,line_size=1.5){
  p <- ggplot() +
    geom_histogram(
      data = data.frame(S = info$var),
      aes(x = S, y = after_stat(density), fill = "Histogram"),
      color = "grey50", bins = 50, alpha = 0.3
    ) +
    geom_line(data = marginal_S_data_by_A, aes(x = x, y = y_reg_npmle, color = "Reg-NPMLE"), size = line_size, lineend = "round") +
    geom_line(data = marginal_S_data_by_A, aes(x = x, y = y_joint_npmle_discrete, color = "Joint-NPMLE"), size = line_size, lineend = "round") +
    xlim(0, xlim_R) +
    labs(x = expression(bold(S[i]^2)), y = "Density") +
    scale_fill_manual(values = c("Histogram" = "grey"), name = NULL) +
    theme_minimal()


  # Build color scale conditionally so legend matches what is plotted
  color_values <- c( "Reg-NPMLE" = "lightskyblue2",
                    "Joint-NPMLE" = "deepskyblue4")
  color_breaks <- c("Reg-NPMLE", "Joint-NPMLE")
  color_labels <- c("Reg-NPMLE","Joint-NPMLE")

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

idx <- which(pep_grp == "1")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "1",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_1 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.2)
plot_marginal_S_by_A_1 <- func_plot_modified(plot_marginal_S_by_A_1) + 
  labs(title = expression(M[i]==1)) + theme(axis.title.x = element_blank()) 

idx <- which(pep_grp == "2")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "2",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_2 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.15)
plot_marginal_S_by_A_2 <- func_plot_modified(plot_marginal_S_by_A_2) + 
  labs(title = expression(M[i]==2)) + theme(legend.position = "none",
                                            axis.title.y = element_blank(),
                                            axis.title.x = element_blank()) 

idx <- which(pep_grp == "3")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "3",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_3 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.1)
plot_marginal_S_by_A_3 <- func_plot_modified(plot_marginal_S_by_A_3) + 
  labs(title = expression(M[i]==3)) + theme(legend.position = "none",
                                            axis.title.y = element_blank(),
                                            axis.title.x = element_blank()) 

idx <- which(pep_grp == "4")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "4",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_4 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.1)
plot_marginal_S_by_A_4 <- func_plot_modified(plot_marginal_S_by_A_4)+ 
  labs(title = expression(M[i]==4)) + theme(legend.position = "none")

idx <- which(pep_grp == "5")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "5",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_5 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.1)
plot_marginal_S_by_A_5 <- func_plot_modified(plot_marginal_S_by_A_5)+ 
  labs(title = expression(M[i]==5)) + theme(legend.position = "none") + theme(axis.title.y = element_blank())


idx <- which(pep_grp == "6")
info_sub <- list(
  var = info$var[idx],
  A   = info$A[idx],
  df  = info$df
)
marginal_S2_data_by_A <- make_marginal_S2_data_by_A(
  info_sub = info_sub,
  prior_result = prior_result,
  group_id = "6",
  mhat_fun = fit_trend_var_A(info$var,info$A)$pred_e
)
plot_marginal_S_by_A_6 <- make_marginal_S_by_A_plot(marginal_S2_data_by_A,info_sub,xlim_R =0.1)
plot_marginal_S_by_A_6 <- func_plot_modified(plot_marginal_S_by_A_6) + 
  labs(title = expression(M[i]==6)) + theme(legend.position = "none") + theme(axis.title.y = element_blank())


panel_plot <- (plot_marginal_S_by_A_1 | plot_marginal_S_by_A_2 | plot_marginal_S_by_A_3) /
              (plot_marginal_S_by_A_4 | plot_marginal_S_by_A_5 | plot_marginal_S_by_A_6)
#ggsave("./figure/marginal_logS2_panel.pdf", panel_plot, width = 16, height = 8)
