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

## Marginal S
marginal_S_data <- make_marginal_S_data(info,prior_result,include_joint = FALSE)
plot_marginal_S <- make_marginal_S_plot(marginal_S_data,info,xlim_R =0.1,include_joint = FALSE)
plot_marginal_S <- func_plot_modified(plot_marginal_S)

## Prior sigma
plot_prior_sigma <- make_prior_sigma_plot(prior_result,scale_factor=200,xlim_R =0.1)
plot_prior_sigma <- func_plot_modified(plot_prior_sigma)

## Marginal V
marginal_V_data <- make_marginal_V_data(info,prior_result)
plot_marginal_V <- make_marginal_V_plot(marginal_V_data,info,xlim_R =7.5)
plot_marginal_V <- func_plot_modified(plot_marginal_V)

## Prior tau
plot_prior_tau <- make_prior_tau_plot(prior_result,scale_factor=1.5,xlim_R =8)
plot_prior_tau <- func_plot_modified(plot_prior_tau)


#ggsave(filename = "./figure/marginal_S.pdf", plot = plot_marginal_S, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_sigma.pdf", plot = plot_prior_sigma, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/marginal_V.pdf", plot = plot_marginal_V, width = 8, height = 6, dpi = 300)
#ggsave(filename = "./figure/prior_tau.pdf", plot = plot_prior_tau, width = 8, height = 6, dpi = 300)


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
prior_result$mimic_joint_npmle<-mimic_joint_npmle_prior

# P value calculation
P_list_mimic_joint_npmle <- P_value_mimic_joint_npmle(info,mimic_joint_npmle_prior,pep_grp,contrast_name=1)
result$mimic_joint_npmle <- P_list_mimic_joint_npmle
sum(BH_adjust(P_list_mimic_joint_npmle,alpha))

#save(info,prior_result,result, file = "data/plot_ready_data_proteomics.RData")

# -----Mimic-joint-npmle plot-----

#load("data/plot_ready_data_proteomics.RData")

## trend plot
x_proteome <- expression(bold(M[i])~"(" * "peptide counts)")
plot_trend_mimic_joint<- make_trend_plot_mimic_joint(info,prior_result,pep_grp,0.6,x_proteome)
plot_trend_mimic_joint <- func_plot_modified(plot_trend_mimic_joint)

#ggsave(filename = "./figure/combined_trend.pdf", plot = plot_trend_mimic_joint, width = 20, height = 6, dpi = 300)

## significance plot
out <- summary_significance_pep_grp(result,pep_grp)
p_bar <- ggplot(out$df_sum, aes(x = M, y = sig, fill = method)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.6,alpha = 0.8,color = "black",linewidth = 0.4) +
  labs(x = expression(bold(M[i])),y = expression(bold("# significant discoveries")),fill = NULL) +
  scale_fill_manual(
    values = c("gray40","mediumpurple4","lightcoral","lightsalmon1","aquamarine3","#377EB8"),
    breaks = c("t","untrend","untrend_1d","trend","trend_1d","trend_2d"),
    labels = c("t-test",expression(bold("Untrended-Inv"*chi^2)), "Untrended-NPMLE", expression(bold("Reg-Inv"*chi^2)),
               "Reg-NPMLE$","Joint-NPMLE")
  ) +
  theme_minimal(base_size = 14)
p_bar_final <- func_plot_modified(p_bar)

#ggsave(filename = "./figure/discoveries.pdf", plot = p_bar_final, width = 20, height = 6, dpi = 300)

## Prop-significance plot

p_prop <- ggplot(out$df_prop, aes(x = M, y = prop, color = method, group = method,shape=method)) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1),expand = expansion(mult = c(0, 0.3))) +
  labs(x = expression(bold(M[i])), y = "Proportion of significance",color = NULL) +
  scale_color_manual(
    values = c("gray40","mediumpurple4","lightcoral","lightsalmon1","aquamarine3","#377EB8"),
    breaks = c("t","untrend","untrend_1d","trend","trend_1d","trend_2d"),
    labels = c("t-test",expression(bold("Untrended-Inv"*chi^2)),"Untrended-NPMLE", expression(bold("Reg-Inv"*chi^2)),
               "Reg-NPMLE","Joint-NPMLE") ) +
  scale_shape_manual(
    breaks = c("t","untrend","untrend_1d","trend","trend_1d","trend_2d"),
    values = c( 17, 15, 3, 18, 7, 8),  # pick any 7 distinct shapes you like
    labels = c("t-test",expression(bold("Untrended-Inv"*chi^2)),"Untrended-NPMLE",expression(bold("Reg-Inv"*chi^2)),    
               "Reg-NPMLE","Joint-NPMLE"),
    guide = "none")+
  guides(color = guide_legend(override.aes = list(shape = unname(c( 17, 15, 3, 18, 7, 8))))) +
  theme_minimal(base_size = 14)
p_prop_final<- func_plot_modified(p_prop) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x  = element_text(size = 10,angle = 30, hjust = 1, vjust = 1))

#ggsave(filename = "./figure/prop_discoveries.pdf", plot = p_prop_final, width = 8, height = 6, dpi = 300)

## Mimic-joint-prior plot
df_spike <- bind_rows(lapply(names(mimic_joint_npmle_prior), function(g) {
  d <- mimic_joint_npmle_prior[[g]]
  data.frame(
    pep_grp = g,
    value   = d$value,
    prob    = d$prob
  )
}))
bins <- unique(df_spike$pep_grp)

make_one_bin_plot <- function(g, show_x = FALSE) {
  d <- df_spike %>% filter(pep_grp == g)

  p <- ggplot(d, aes(x = value)) +
    geom_segment(aes(xend = value, y = 0, yend = prob),
                 linewidth = 1, color = "aquamarine3") +
    coord_cartesian(xlim = c(5e-4, 1e-1),
                    ylim = c(0, max(df_spike$prob, na.rm = TRUE))) +
    scale_x_log10() +
    scale_y_continuous(breaks = NULL) +   # <-- removes 0.0/0.25/... everywhere
    labs(y = g, x = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin = margin(4, 10, 4, 10),
      axis.title.y = element_text(size = 9, face = "plain", angle = 0, vjust = 0.5),

      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),

      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  } else {
    p <- p + labs(x = expression(bold(sigma^2))) +
      theme(axis.title.x = element_text(face = "bold", size = 16),
            axis.text.x  = element_text(size = 12))
  }

  p
}

plist <- lapply(seq_along(bins), \(j) make_one_bin_plot(bins[j], show_x = (j == length(bins))))
p_stack <- wrap_plots(plist, ncol = 1) + plot_layout(heights = rep(1, length(plist)))
ylab_col <- wrap_elements(full = textGrob(expression(M[i]), rot = 90,
                                         gp = gpar(fontface = "bold", fontsize = 14)))
p_with_ylab <- ylab_col + p_stack + plot_layout(widths = c(0.08, 1))

add_panel_border <- function(p, col = "black", lwd = 2) {
  g <- patchworkGrob(p)

  panels <- grepl("^panel", g$layout$name)
  t <- min(g$layout$t[panels]); b <- max(g$layout$b[panels])
  l <- min(g$layout$l[panels]); r <- max(g$layout$r[panels])

  g <- gtable_add_grob(
    g,
    rectGrob(gp = gpar(fill = NA, col = col, lwd = lwd)),
    t = t, b = b, l = l, r = r, z = Inf
  )
  g
}

g_stack_bordered <- add_panel_border(p_stack, lwd = 2)
p_stack_bordered  <- wrap_elements(full = g_stack_bordered)
mimic_prior_final <- ylab_col + p_stack_bordered + plot_layout(widths = c(0.01, 1))
#ggsave(filename = "./figure/mimic-joint-prior.pdf", plot = mimic_prior_final, width = 4.5, height = 11, dpi = 300)
