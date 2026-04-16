library(ggplot2)
library(dplyr)
library(patchwork)


out_dir <- "camera_ready_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----Load data-----
load("./trend_data.RData")

# -----Plot Function-----
make_trend_panel <- function(trend_data, alpha, xlab, ylab = expression(log(S[i]^2))) {
  td <- trend_data %>% arrange(M)
  ggplot(td, aes(x = M)) +
    geom_point(
      aes(y = log_sample_var),
      shape = 16,
      color = "#AAAAAA",
      stroke = 0, 
      size = 0.5,
      alpha = alpha,
    ) +
    geom_line(
      aes(y = log_var_fitted, color = "Fitted trend"),
      size = 0.9,
      lineend = "round"
    ) +
    geom_hline(
      aes(yintercept = log_var_fitted_untrended, color = "Constant trend"),
      size = 0.9,
      linetype = "22"
    ) +
    scale_color_manual(
      values = c("Fitted trend" = "#D55E00", "Constant trend" = "#0072B2"),
      ,
      breaks = c("Fitted trend", "Constant trend"),
      name = NULL
    ) +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 10) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 9),
      legend.position = c(0.7, 0.2),
      legend.text = element_text(size = 9),
      legend.key.width = unit(2, "lines"),
      legend.key.height = unit(1.2, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(1, 1, 1, 1, unit = "pt")
    )
}
save_panel <- function(p, stem, w = 3.15, h = 2.75) {
  ggsave(
    filename = file.path(out_dir, paste0(stem, ".pdf")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE
  )
  ggsave(
    filename = file.path(out_dir, paste0(stem, ".png")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# -----Plot-----

out_dir <- "camera_ready_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

x_rnaseq   <- expression(bold(A[i])~"(avg. intensity)")
x_chipseq  <- expression(bold(A[i])~"(avg. intensity)")
x_proteome <- expression(bold(M[i])~"(" * log[2] * " peptide counts)")

q <- quantile(rnaseq1_trend_data$M,0.999)
rnaseq1_trend_data_filter <- filter(rnaseq1_trend_data, M <= q)
p_rnaseq1 <- make_trend_panel(rnaseq1_trend_data_filter, 0.2, xlab = x_rnaseq)

p_rnaseq2 <- make_trend_panel(rnaseq2_trend_data, 0.8, xlab = x_rnaseq, ylab=NULL) +  theme(legend.position = "none")

q_chipseq_low <- quantile(chipseq_trend_data$M, 0)
q_chipseq_hi <- quantile(chipseq_trend_data$M, 0.99)
chipseq_trend_data_filt <- filter(chipseq_trend_data, M <= q_chipseq_hi, M >= q_chipseq_low)
p_chipseq   <- make_trend_panel(chipseq_trend_data_filt, 0.17, xlab = x_chipseq, ylab=NULL) + theme(legend.position = "none")

q_proteome <- quantile(proteomic_trend_data$M, 0.999)
proteomic_trend_data_filt <- filter(proteomic_trend_data, M <= q_proteome)

p_proteome  <- make_trend_panel(proteomic_trend_data_filt,0.6, xlab = x_proteome, ylab=NULL) + theme(legend.position = "none") + scale_x_continuous(trans = scales::log2_trans())

#save_panel(p_rnaseq1,   "panel1rnaseq1", w = 4.2, h = 3)
#save_panel(p_rnaseq2,   "panel2rnaseq2", w = 4, h = 3)
#save_panel(p_chipseq,   "panel3chipseq", w = 4, h = 3)
#save_panel(p_proteome,  "panel4proteomic", w = 4, h = 3)
