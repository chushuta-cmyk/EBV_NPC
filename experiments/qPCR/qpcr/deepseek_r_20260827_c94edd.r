library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# --------------------------------------------
# 1. Panel A: BZLF1 (REAL data)
# --------------------------------------------
bzlf1 <- data.frame(
  Group = factor(c("DMSO1", "TPA1", "DMSO2", "TPA2"),
                 levels = c("DMSO1", "TPA1", "DMSO2", "TPA2")),
  RQ = c(1.00, 2.59, 2.65, 7.25),
  CI_low = c(0.20, 1.85, 1.57, 4.53),
  CI_high = c(4.95, 3.63, 4.48, 11.62)
)

pA <- ggplot(bzlf1, aes(x = Group, y = RQ, fill = Group)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4) +   # use linewidth
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.6) +
  geom_text(aes(label = round(RQ, 2)), vjust = -0.7, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("grey70", "#E41A1C", "grey70", "#E41A1C")) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2, 5, 10, 20),
                labels = c("0.2", "0.5", "1", "2", "5", "10", "20")) +
  labs(y = "BZLF1 RQ (log scale)", subtitle = "A") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 16))

# --------------------------------------------
# 2. Panel B: Mono-culture IL1A + CXCL8 (with CI)
# --------------------------------------------
# Add CI columns (real values from your raw data)
mono_data <- data.frame(
  Gene = rep(c("IL1A", "CXCL8"), each = 2),
  Condition = rep(c("DMSO", "TPA"), 2),
  RQ = c(1.00, 2.07, 1.00, 1.47),
  CI_low = c(0.37, 1.37, 0.36, 0.99),     # add these
  CI_high = c(2.68, 3.13, 2.74, 2.20)     # add these
)
mono_data$Gene <- factor(mono_data$Gene, levels = c("IL1A", "CXCL8"))

pB <- ggplot(mono_data, aes(x = Condition, y = RQ, fill = Condition)) +
  geom_col(width = 0.5, position = position_dodge(0.6), color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.1, position = position_dodge(0.6), linewidth = 0.6) +
  geom_text(aes(label = round(RQ, 2)), 
            position = position_dodge(0.6), 
            vjust = -0.8, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("DMSO" = "grey70", "TPA" = "#377EB8")) +
  facet_wrap(~ Gene, scales = "free_y") +
  coord_cartesian(ylim = c(0, 3.5)) +   # gives headroom for error bars
  labs(y = "RQ", subtitle = "B") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(face = "bold", size = 16))

# --------------------------------------------
# 3. Panel C: Co-culture paired scatter (REAL data)
# --------------------------------------------
co_data <- data.frame(
  Donor = rep(c("Donor 1", "Donor 2"), each = 6),
  Gene = rep(rep(c("TREM1", "IL1A", "CXCL8"), each = 2), times = 2),
  Condition = rep(c("DMSO", "TPA"), times = 6),
  RQ = c(
    # Donor 1
    1.00, 0.93,   # TREM1
    1.00, 2.07,   # IL1A
    1.00, 1.47,   # CXCL8
    # Donor 2
    0.61, 0.30,   # TREM1
    1.85, 2.11,   # IL1A
    1.34, 1.19    # CXCL8
  )
)
co_data$Gene <- factor(co_data$Gene, levels = c("TREM1", "IL1A", "CXCL8"))

pC <- ggplot(co_data, aes(x = Condition, y = RQ, 
                          group = interaction(Donor, Gene), color = Donor)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Donor 1" = "#1B9E77", "Donor 2" = "#D95F02")) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(y = "RQ", subtitle = "C", color = "PBMC donor") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(face = "bold", size = 16))

# --------------------------------------------
# 4. Combine all three panels
# --------------------------------------------
figure2 <- (pA | pB) / pC + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 18))

# Save
ggsave("Figure2_3Panel_RealData.pdf", figure2, width = 12, height = 10, dpi = 300)
ggsave("Figure2_3Panel_RealData.png", figure2, width = 12, height = 10, dpi = 300)

print("✅ Figure 2 (3 panels) generated successfully with error bars for Panel B!")
