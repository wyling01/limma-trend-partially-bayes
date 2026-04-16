library(oligo)
library(edgeR)
library(readxl)
library(GEOquery)
library(quadprog)
library(ruv)
library(EDASeq)
library(RColorBrewer)
library(VennDiagram)
library(data.table)
library(pheatmap)

# -----Load Function-----
source('../function/func_limma_trend.R')
source('../function/func_plot.R')
# -----Data Preprocessing-----

## ---Load gene ID mappings
geneID_mappings <- fread("./data/geneID_mappings.txt", data.table = FALSE, stringsAsFactors = FALSE)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
s <- str_split(as.character(geneID_mappings$`[Previous ID(s)]`), ',')
geneID_mappings <- data.frame(current=rep(geneID_mappings$`[Gene ID]`, sapply(s, FUN=length)), old=unlist(s))
geneID_mappings$old <- str_trim(geneID_mappings$old)

## ---Load read counts
load("./data/fc_genes_with_all_samples_aligned_human_Pf_WithOutVivax_subread_uH.RData")
fc2 <- fc_genes_aligned_w_Pfhuman_withoutVivax
t <- fc2$counts

x <- DGEList(counts=fc2$counts, genes=fc2$annotation)

## ---Initial filter
# filter out genes that do not have sufficient coverage across our samples for reliable estimates to be made. 
before <- nrow(x$counts)
keep <- rowSums(cpm(x)>2) >= 10
x <- x[keep,,keep.lib.sizes=FALSE]
bplot <- reshape2::melt(colSums(x$counts))
colnames(bplot) <- c("Read.Counts")
bplot$Sample <- rownames(bplot)

gg <- ggplot(bplot, aes(x=factor(Sample), y=Read.Counts)) + geom_bar(stat = "identity")
gg <- gg + theme_bw() + scale_y_sqrt(breaks=c(0,10000,1000000,10000000,100000000))
gg <- gg + theme(axis.text.x = element_text(size=12, angle = 90)
                 , axis.text.y = element_text(size=12, angle = 0)
                 , axis.title=element_text(size=14,face="bold"))
gg <- gg + labs(x='Sample', y='Read Count') 
gg <- gg + geom_hline(aes(yintercept=1000000), col="red")

x$counts <- x$counts[,!grepl("SFD1|SFM9|SFC\\.025", colnames(x$counts))]
filt.annnotation <- fc2$annotation[fc2$annotation$GeneID %in% rownames(x$counts),]
dge <- DGEList(counts=x$counts, genes=filt.annnotation)

gg
# Samples SFD1, SFM9 and SFC.025 were not sequenced to a sufficient depth to be reliable and thus are removed.

## ---Remove drug treated samples
x$counts <- x$counts[,!grepl("SFU2|SFU\\.3", colnames(x$counts))]
x$counts <- x$counts[,!grepl("SFC\\.023|SFM\\.7|IFM12|IFM21", colnames(x$counts))]
filt.annnotation <- fc2$annotation[fc2$annotation$GeneID %in% rownames(x$counts),]
dge <- DGEList(counts=x$counts, genes=filt.annnotation)
# Now we take the log transform of our sample counts (RPKM) to estimate the proportion of different stages present.
our_log_rpkm <- rpkm(dge)
our_log_rpkm <- log2(1+our_log_rpkm)


## ---Estimate Contributions from Different Life-cycle Stages
su_rpkm <- read.csv("./data/_rpkm_7SexualAndAsexualLifeStages_suetal.csv", header=TRUE, sep=",")

#ignore the Ookinete profile as this should be absent from our data.
su_rpkm$Ookinete <- NULL
rownames(su_rpkm) <- su_rpkm$ID
su_rpkm$ID <- NULL

#take the log transform of the data which helps to manage outlier genes.
su_log_rpkm <- log2(1+su_rpkm)

#use mixture model
findMix <- function(Y, X){  
  X[is.na(X)] <- t(replicate(ncol(X), apply(X,1,mean, na.rm=T)))[is.na(X)]
  Rinv <- solve(chol(t(X) %*% X))
  C <- cbind(rep(1,ncol(X)), diag(ncol(X)))
  b <- c(1,rep(0,ncol(X)))
  d <- t(Y) %*% X  
  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq=1)
  sol <- qp$solution
  sol[sol<1e-10] <- 0
  return(sol)
}
#Calculate Sample Life-Cycle Profiles
inter <- intersect(rownames(our_log_rpkm), rownames(su_log_rpkm))

O <- our_log_rpkm[rownames(our_log_rpkm) %in% inter, ]
O <- O[order(rownames(O)),]
S <- su_log_rpkm[rownames(su_log_rpkm) %in% inter, ]
S <- S[order(rownames(S)), ]

#Now lets fit some samples!
ourPlotData <- data.frame()
for (i in 1:ncol(O)){
  mix <- findMix(O[,i], as.matrix(S))
  ourPlotData <- rbind(ourPlotData, data.frame(sample=rep(colnames(O)[i], ncol(S)), stage=colnames(S), proportion=mix))
}

#Organise the results.
ourPlotData$stage <- gsub("Gametocyte.*","Gametocyte",ourPlotData$stage)
ourPlotData <- aggregate(proportion~sample+stage,data=ourPlotData,FUN=sum)
ourPlotData <- within(ourPlotData, stage <- factor(stage, levels=c("Ring", "Early.Trophozoite", "Late.Trophozoite", "Schizont","Gametocyte")))
ourPlotData$phenotype <- ifelse(substring(ourPlotData$sample,1,1)=="I" ,"non-severe", "severe")

#Make a plot.
gg <- ggplot(ourPlotData, aes(x=factor(sample), y=proportion, fill=factor(phenotype))) + geom_bar(stat="identity")
gg <- gg + scale_fill_manual(values = c("non-severe" = "#2c7bb6", "severe" = "#d7191c"))
gg <- gg + facet_wrap(~ stage, ncol = 1)
gg <- gg + theme_bw()
gg <- gg + theme(axis.text.x = element_text(size=12, angle = 90)
                 , axis.text.y = element_text(size=12, angle = 0)
                 , axis.title=element_text(size=14,face="bold")
                 , strip.text.x = element_text(size = 16, face="bold"))
gg <- gg + labs(x='Sample', y='Proportion', fill='Stage')
gg <- gg + theme(legend.text=element_text(size=14))
gg <- gg + theme(legend.key.size =  unit(0.25, "in"))
gg <- gg + theme(legend.title = element_text(size=16, face="bold"))
gg <- gg + guides(fill=guide_legend(title="Phenotype"))
gg

#Investigate Differential Expression
options(width=150)
categories<- as.factor(substring(rownames(dge$samples),1,1))
dge <- calcNormFactors(dge, method="TMM")
dge$samples$group <- categories
design <- model.matrix(~group, data=dge$samples)
v <- voom(dge, design=design, plot=TRUE)


# PCA Normalising for Library Size, Ring, Gametocyte and Schizont effects.

covs <- data.frame(v$design[,2])
ourPlotData$phenotype <- NULL
c <-  reshape2::dcast(ourPlotData, sample ~ ..., value.var="proportion")
covs <- merge(covs, c, by.x=0, by.y="sample")

colnames(covs) <- c("sample", "disease", colnames(c)[2:ncol(c)])
rownames(covs) <- covs$sample
covs$sample <- NULL
covs <- covs[match(colnames(v$E), rownames(covs)),]
stopifnot(colnames(v$E)==rownames(covs))
# head(covs)

covs <- covs[,c(1,2,5,6)]
# head(covs)

mod = model.matrix(
  as.formula(paste("~", paste(colnames(covs), collapse=" + "), sep="")), data=covs)

norm_counts_ring <- removeBatchEffect(v$E, covariates=mod[,3:ncol(mod)], design=mod[,1:2])
colors <- c("#0571b0","#ca0020")
plotPCA(norm_counts_ring, col=colors[categories], cex=0.8, isLog=TRUE)

covs <- covs[,c(1,2)]

mod = model.matrix(
  as.formula(paste("~", paste(colnames(covs), collapse=" + "), sep="")), data=covs)

norm_counts_ring <- removeBatchEffect(v$E, covariates=mod[,3:ncol(mod)], design=mod[,1:2])
plotPCA(norm_counts_ring, col=colors[categories], cex=0.8, isLog=TRUE)

contrast <- makeContrasts(Treatment = disease, levels = mod)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

## Save data
data_rnaseq2_malaria <- list(Y=logCPM,design = mod,contrast=contrast)

save(data_rnaseq2_malaria, file = "./data/filtered_data_rnaseq2.RData")




# -----DE analysis-----
load("./data/filtered_data_rnaseq2.RData")

info <- info_extractor(data_rnaseq2_malaria)

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
sum(BH_adjust(result$t_test,alpha)) #13
sum(BH_adjust(result$untrended_invchi,alpha)) #77
sum(BH_adjust(result$untrended_npmle,alpha)) #11
sum(BH_adjust(result$reg_invchi,alpha)) # 6
sum(BH_adjust(result$reg_npmle,alpha)) # 6
sum(BH_adjust(result$joint_npmle,alpha)) # 10

#save(info,prior_result, file = "data/plot_ready_data_rnaseq2.RData")

# Plot

#load("data/plot_ready_data_rnaseq2.RData")

## Trend plot
trend_data <- make_trend_data(info,prior_result)
x_rnaseq <- expression(bold(A[i])~"(avg. intensity)")
plot_trend<- make_trend_plot(trend_data,0.7,x_rnaseq)
plot_trend <- func_plot_modified(plot_trend)

## Marginal S
marginal_S_data <- make_marginal_S_data(info,prior_result)
plot_marginal_S <- make_marginal_S_plot(marginal_S_data,info,xlim_R =5)
plot_marginal_S <- func_plot_modified(plot_marginal_S)

## Prior sigma
plot_prior_sigma <- make_prior_sigma_plot(prior_result,scale_factor=5,xlim_R =6)
plot_prior_sigma <- func_plot_modified(plot_prior_sigma)

## Marginal V
marginal_V_data <- make_marginal_V_data(info,prior_result)
plot_marginal_V <- make_marginal_V_plot(marginal_V_data,info,xlim_R =6)
plot_marginal_V <- func_plot_modified(plot_marginal_V)

## Prior tau
plot_prior_tau <- make_prior_tau_plot(prior_result,scale_factor=3,xlim_R =6)
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
