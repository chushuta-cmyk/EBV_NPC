library(ggplot2)
library(patchwork)

# ========== 真实数据 ==========
# BZLF1（来自上一轮数据）
bzlf1 <- data.frame(
  Group = c("DMSO1", "TPA1", "DMSO2", "TPA2"),
  RQ = c(1.00, 2.59, 2.65, 7.25),
  Min = c(0.20, 1.85, 1.57, 4.53),
  Max = c(4.95, 3.63, 4.48, 11.62)
)

# TREM1（真实）
trem1 <- data.frame(
  Group = c("DMSO1", "TPA1", "DMSO2", "TPA2"),
  RQ = c(1.00, 0.93, 0.61, 0.30),
  Min = c(0.29, 0.47, 0.23, 0.22),
  Max = c(3.46, 1.83, 1.66, 0.41)
)

# IL1A（真实）
il1a <- data.frame(
  Group = c("DMSO1", "TPA1", "DMSO2", "TPA2"),
  RQ = c(1.00, 2.07, 1.85, 2.11),
  Min = c(0.37, 1.37, 0.67, 1.57),
  Max = c(2.68, 3.13, 5.12, 2.82)
)

# CXCL8（真实）
cxcl8 <- data.frame(
  Group = c("DMSO1", "TPA1", "DMSO2", "TPA2"),
  RQ = c(1.00, 1.47, 1.34, 1.19),
  Min = c(0.36, 0.99, 0.50, 0.86),
  Max = c(2.74, 2.20, 3.64, 1.65)
)

# ========== 绘图函数 ==========
plot_gene <- function(data, y_lab, fill_colors) {
  ggplot(data, aes(x = Group, y = RQ, fill = Group)) +
    geom_col(width = 0.6, color = "black", size = 0.3) +
    geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.15, size = 0.5) +
    geom_text(aes(label = round(RQ, 2)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = fill_colors) +
    labs(y = y_lab) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_blank())
}

# 颜色：对照灰，处理红
cols <- c("grey70", "#D55E00", "grey70", "#D55E00")

pA <- plot_gene(bzlf1, "BZLF1 RQ (log scale)", cols) + scale_y_log10()
pB <- plot_gene(il1a, "IL1A RQ", cols)
pC <- plot_gene(cxcl8, "CXCL8 RQ", cols)
pD <- plot_gene(trem1, "TREM1 RQ", cols)

# 拼图
figure2 <- (pA | pB) / (pC | pD)
ggsave("Figure2_Updated.pdf", figure2, width = 8, height = 6)