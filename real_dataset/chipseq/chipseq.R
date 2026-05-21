library(oligo)
library(edgeR)
library(readxl)


# -----Load Function-----
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
# -----Data Preprocessing-----

## Read data.
profile_bins <- read.table("./data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
n1 <- 3
n <- 6
normalized <- read.table("./data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)
head(normalized)
rownames(normalized) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                              sep = "_")


## Construct the design matrix.
design <- cbind(cond1 = 1, cond2_minus_cond1 = c(rep(0, n1), rep(1, n - n1)))
rownames(design) <- names(normalized)
contrast <- makeContrasts(Treatment = cond2_minus_cond1, levels = design)

## Save data
data_chipseq <- list(Y=as.matrix(normalized),design = design,contrast=contrast)

#save(data_chipseq, file = "./data/filtered_data_chipseq.RData")


# DE analysis
#load("./data/filtered_data_chipseq.RData")



info <- info_extractor(data_chipseq)



## Prior estimation
param_prior <- var_invchi_prior(info)
untrended_invchi_prior <- param_prior[[1]]
reg_invchi_prior <- param_prior[[2]]

untrended_npmle_prior <- prior_untrended_npmle(info,v=300)
reg_npmle_prior <- prior_reg_npmle(info,v=300)

joint_npmle_prior <- joint_npmle_2d(info,pbin=100,pv=30,qRes_L = 0.1)


threshold <- 1e-8
idx <- which(joint_npmle_prior$mass>threshold)
joint_npmle_prior$idx <- idx
prior_result <- list(
  untrended_invchi=untrended_invchi_prior,
  untrended_npmle=untrended_npmle_prior,
  reg_invchi =reg_invchi_prior,
  reg_npmle=reg_npmle_prior,
  joint_npmle=joint_npmle_prior
)


## P value calculation
result <- p_value_calculator(info,prior_result,verbose=TRUE)

alpha = 0.0001
sum(BH_adjust(result$t_test,alpha)) # 0
sum(BH_adjust(result$untrended_invchi,alpha)) # 197
sum(BH_adjust(result$untrended_npmle,alpha)) # 2525
sum(BH_adjust(result$reg_invchi,alpha)) # 4638
sum(BH_adjust(result$reg_npmle,alpha)) # 4620
sum(BH_adjust(result$joint_npmle,alpha)) # 5674

# Plot

## Trend plot
trend_data <- make_trend_data(info,prior_result)
plot_trend <- make_trend_plot(trend_data)
plot_trend <- func_plot_modified(plot_trend)

## Log Marginal
logvar_data <- make_logvar_marginal_data(
  info = info,
  prior_result = prior_result,
  length.out = 3000,
  xlim_log = c(-10, 10)
)
plot_marginal_combined <-  plot_logvar_marginal_data(plot_data = logvar_data,bins = 100)
plot_marginal_combined <- func_plot_modified(plot_marginal_combined ) +ylim(0,0.5)

## Log Prior
prior_log_data <- make_npmle_prior_log_data(
  prior_result = prior_result,
  threshold = 1e-8,
  length.out = 5000,
  xlim_log = c(-8, 8)
)
plot_prior_combined  <- plot_npmle_prior_log_data(plot_data = prior_log_data)
plot_prior_combined <- func_plot_modified(plot_prior_combined) 

## Joint prior
plot_joint_prior <- make_joint_prior_plot(prior_result)
plot_joint_prior <- func_plot_modified(plot_joint_prior)+theme(legend.position = "none")



#ggsave(filename = "./figure/log_marginal_combined.pdf", plot = plot_marginal_combined , width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/log_prior_combined.pdf", plot = plot_prior_combined, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_2d.pdf", plot = plot_joint_prior, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/combined_trend.png", plot = plot_trend, width = 8, height = 6, dpi = 300)
