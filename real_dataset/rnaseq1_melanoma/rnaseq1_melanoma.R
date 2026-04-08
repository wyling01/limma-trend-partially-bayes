library(oligo)
library(GEOquery)
library(edgeR)
library(readxl)


# Load Function
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
# Data Preprocessing

## Get raw data
raw_df <- read_excel("./data/GSE114716_raw.counts.hs.xlsx")
gset <- getGEO("GSE114716", GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gset[[1]])

## Extract design from the metadata

### Extract Patient ID

metadata$Patient <- str_extract(metadata$title, "patient \\d+")
metadata$Patient <- gsub("patient ", "P", metadata$Patient)
metadata$Patient <- factor(metadata$Patient)

### Extract Treatment Group
metadata$Treatment <- ifelse(grepl("baseline", metadata$title, ignore.case=TRUE), 
                             "Control", 
                             "Treated")
metadata$Treatment <- factor(metadata$Treatment, levels = c("Control", "Treated"))
design <- model.matrix(~ Treatment + Patient, data = metadata)

### Extract contrast
contrast <- makeContrasts(Treatment = TreatmentTreated, levels = design)

## Filtering steps

### Create the DGEList object
dge <- DGEList(counts = raw_df, samples = metadata)


### Run and apply the Adaptive Filter
### This calculates the minimum CPM needed to have ~10 reads in the smallest group
keep <- filterByExpr(dge, design)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

### edgeR normalization and dispersion estimation
dge_filtered<- calcNormFactors(dge_filtered)
dge_filtered <- estimateCommonDisp(dge_filtered)
dge_filtered <- estimateGLMTrendedDisp(dge_filtered)
dge_filtered <- estimateTagwiseDisp(dge_filtered)
plotBCV(dge_filtered)

### Compute log2 counts-per-million
logCPM <- cpm(dge_filtered, log = TRUE, prior.count = 3)

### Check if there is NA in the expression matrix. If yes, remove rows that contain NAs.
any_na <- anyNA(logCPM)

if (any_na) {
  keep <- complete.cases(logCPM)  
  logCPM <- logCPM[keep, , drop = FALSE]
}

## Save data
data_rnaseq1_melanoma <- list(Y=logCPM,design = design,contrast=contrast)

save(data_rnaseq1_melanoma, file = "./data/filtered_data_rnaseq1.RData")


# DE analysis
load("./data/filtered_data_rnaseq1.RData")



info <- info_extractor(data_rnaseq1_melanoma)



## Prior estimation
param_prior <- var_invchi_prior(info)
untrended_invchi_prior <- param_prior[[1]]
reg_invchi_prior <- param_prior[[2]]

untrended_npmle_prior <- prior_untrended_npmle(info,v=300)
reg_npmle_prior <- prior_reg_npmle(info,v=300)

joint_npmle_prior <- joint_npmle_2d(info,pbin=50,pv=50,qRes_L = 0.1)


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

alpha = 0.05
sum(BH_adjust(result$t_test,alpha)) # 2
sum(BH_adjust(result$untrended_invchi,alpha)) # 0
sum(BH_adjust(result$untrended_npmle,alpha)) # 0
sum(BH_adjust(result$reg_invchi,alpha)) # 133
sum(BH_adjust(result$reg_npmle,alpha)) # 76
sum(BH_adjust(result$joint_npmle,alpha)) # 106

#save(info,prior_result, file = "data/plot_ready_data_rnaseq1.RData")

# Plot

#load("data/plot_ready_data_rnaseq1.RData")

## Trend plot
x_rnaseq <- expression(bold(A[i])~"(avg. intensity)")
plot_trend<- make_trend_plot(trend_data,0.2,x_rnaseq)
plot_trend <- make_trend_plot(trend_data)
plot_trend <- func_plot_modified(plot_trend)

## Marginal S
marginal_S_data <- make_marginal_S_data(info,prior_result)
plot_marginal_S <- make_marginal_S_plot(marginal_S_data,info,xlim_R =0.5)
plot_marginal_S <- func_plot_modified(plot_marginal_S)

## Prior sigma
plot_prior_sigma <- make_prior_sigma_plot(prior_result,scale_factor=50,xlim_R =1)
plot_prior_sigma <- func_plot_modified(plot_prior_sigma)

## Marginal V
marginal_V_data <- make_marginal_V_data(info,prior_result)
plot_marginal_V <- make_marginal_V_plot(marginal_V_data,info,xlim_R =10)
plot_marginal_V <- func_plot_modified(plot_marginal_V)

## Prior tau
plot_prior_tau <- make_prior_tau_plot(prior_result,scale_factor=2,xlim_R =7.5)
plot_prior_tau <- func_plot_modified(plot_prior_tau)

## Joint prior
plot_joint_prior <- make_joint_prior_plot(prior_result)
plot_joint_prior <- func_plot_modified(plot_joint_prior)+theme(legend.position = "none")



#ggsave(filename = "./figure/marginal_S.pdf", plot = plot_marginal_S, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_sigma.pdf", plot = plot_prior_sigma, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/marginal_V.pdf", plot = plot_marginal_V, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_tau.pdf", plot = plot_prior_tau, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_2d.pdf", plot = plot_joint_prior, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/combined_trend.png", plot = plot_trend, width = 8, height = 6, dpi = 300)
