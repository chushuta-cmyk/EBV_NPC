# -------------------------------------------------------------------------
# Nature / Lancet Style Cytokine Barplot
# Minimal publication-ready version
# -------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. Read data
# -------------------------------------------------------------------------

# 读取原始数据
data_raw <- read.csv("cytokine_replicate.csv", stringsAsFactors = FALSE)

# ==========================================
# 新增修改：精选 SCI 核心双极化细胞因子
# ==========================================
# 定义你需要展示的目标因子（需严格与你的 CSV 列名前缀一致）
target_cytokines <- c("CXCL10", "CXCL9", "CCL5", "IFN-g", "IL-2", "CXCL8", "CXCL1", "IL-6", "VEGF")

# 提取所有的列名
all_cols <- colnames(data_raw)

# 第一列是 Condition，必须保留。对于后续列，通过正则去掉 _R1, _R2 后缀，匹配目标因子
selected_cols <- all_cols[1] # 保留 "Condition" 列
for (col in all_cols[-1]) {
  base_name <- sub("_R[0-9]+$", "", col) # 去掉下划线和数字
  if (base_name %in% target_cytokines) {
    selected_cols <- c(selected_cols, col)
  }
}

# 覆盖原数据集，只保留精选列
data_raw <- data_raw[, selected_cols, drop = FALSE]


# -------------------------------------------------------------------------
# 2. Extract cytokine names
# -------------------------------------------------------------------------

col_names <- colnames(data_raw)

cytokine_bases <- unique(
  sub("_R[0-9]+$", "", col_names)
)

# -------------------------------------------------------------------------
# 3. Build long-format replicate dataframe
# -------------------------------------------------------------------------

replicate_list <- list()

for (cyt in cytokine_bases) {
  
  rep_cols <- grep(
    paste0("^", cyt, "_R[0-9]+$"),
    col_names,
    value = TRUE
  )
  
  if (length(rep_cols) > 0) {
    
    for (cond in rownames(data_raw)) {
      
      for (i in seq_along(rep_cols)) {
        
        replicate_list[[length(replicate_list) + 1]] <- data.frame(
          Condition = cond,
          Cytokine = cyt,
          Replicate = i,
          Value = data_raw[cond, rep_cols[i]]
        )
      }
    }
  }
}

df_replicates <- bind_rows(replicate_list)

# -------------------------------------------------------------------------
# 4. Summary statistics
# -------------------------------------------------------------------------

df_summary <- df_replicates %>%
  group_by(Condition, Cytokine) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD   = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# 5. Select Top 8 cytokines by TPA / UT fold-change
# -------------------------------------------------------------------------

fc_table <- df_summary %>%
  pivot_wider(
    id_cols = Cytokine,
    names_from = Condition,
    values_from = Mean
  ) %>%
  mutate(
    FC_TPA_vs_UT = TPA / UT
  ) %>%
  arrange(desc(FC_TPA_vs_UT))

top_cytokines <- fc_table$Cytokine[1:min(8, nrow(fc_table))]

# -------------------------------------------------------------------------
# 6. Filter plotting data
# -------------------------------------------------------------------------

df_plot <- df_summary %>%
  filter(Cytokine %in% top_cytokines) %>%
  mutate(
    Cytokine = factor(Cytokine, levels = top_cytokines)
  )

df_replicates_plot <- df_replicates %>%
  filter(Cytokine %in% top_cytokines) %>%
  mutate(
    Cytokine = factor(Cytokine, levels = top_cytokines)
  )

# -------------------------------------------------------------------------
# Publication-style cytokine plot
# -------------------------------------------------------------------------

p <- ggplot(
  df_plot,
  aes(x = Condition, y = Mean, fill = Condition)
) +
  
  # bars
  geom_col(
    width = 0.66,
    colour = "black",
    linewidth = 0.4,
    alpha = 0.95
  ) +
  
  # error bars
  geom_errorbar(
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.22,
    linewidth = 0.7,
    colour = "black"
  ) +
  
  # facet panels
  facet_wrap(
    ~ Cytokine,
    scales = "free_y",
    ncol = 4
  ) +
  
  # softer publication colors
  scale_fill_manual(
    values = c(
      "UT"   = "#6AAED6",
      "DMSO" = "#B3B3B3",
      "TPA"  = "#F28E00",
      "MOCK" = "#2C7FB8",
      "Zta"  = "#18A87D"
    )
  ) +
  
  # slightly expand y-axis for better SD visibility
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.06))
  ) +
  
  # labels
  labs(
    title = "Top Cytokines Responding to TPA Treatment",
    subtitle = "Mean ± SD (n = 3). TPA shown in orange",
    x = NULL,
    y = "Concentration (pg/mL)"
  ) +
  
  # base theme
  theme_bw(base_size = 11) +
  
  theme(
    
    # light grey panel background
    panel.background = element_rect(
      fill = "#F7F7F7",
      colour = NA
    ),
    
    # grids
    panel.grid.major = element_line(
      colour = "#D9D9D9",
      linewidth = 0.45
    ),
    
    panel.grid.minor = element_line(
      colour = "#EEEEEE",
      linewidth = 0.28
    ),
    
    # panel border
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    ),
    
    # facet strips
    strip.background = element_rect(
      fill = "#E6E6E6",
      colour = "black",
      linewidth = 0.45
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10,
      colour = "black"
    ),
    
    # title
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    # subtitle
    plot.subtitle = element_text(
      size = 11,
      face = "italic",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    # x-axis text
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      colour = "black"
    ),
    
    # y-axis text
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    
    # axis title
    axis.title.y = element_text(
      size = 13,
      face = "bold"
    ),
    
    # axis lines
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    # ticks
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    
    # spacing
    panel.spacing = unit(0.6, "lines"),
    
    # remove legend
    legend.position = "none",
    
    # margins
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

# show plot
print(p)
# -------------------------------------------------------------------------
# 10. Save
# -------------------------------------------------------------------------

ggsave(
  "Nature_style_cytokine_barplot.pdf",
  p,
  width = 8.5,
  height = 5.8,
  useDingbats = FALSE
)

ggsave(
  "Nature_style_cytokine_barplot.tiff",
  p,
  width = 8.5,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# -------------------------------------------------------------------------
# 11. Output selected cytokines
# -------------------------------------------------------------------------

cat("\nTop cytokines by TPA/UT fold-change:\n")

print(
  fc_table[1:8, c("Cytokine", "FC_TPA_vs_UT")]
)

