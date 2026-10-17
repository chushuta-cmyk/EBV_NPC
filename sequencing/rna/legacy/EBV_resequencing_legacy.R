#Based on GSE224450(https://pmc.ncbi.nlm.nih.gov/articles/PMC10184620/)
######数据下载与加载
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(pathview)
library(CIBERSORT)
library(tibble)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(pROC)
library(ROCR)
library(xgboost)
library(glmnet)
library(randomForest)
library(kernlab)
library(gbm)
library(e1071)
library(ggpubr)
library(vip)
library(DALEX)
library(reshape2)
library(survival)
library(survminer)
library(venn)
library(VennDiagram)
library(corrplot)
library(viridis)
library(RColorBrewer)
library(tibble)

#设置工作目录和保存路径
working_dir <- getwd()  # 获取当前工作目录
tables_dir <- file.path(working_dir, "tables")  # 定义 tables 目录路径
graphs_dir <- file.path(working_dir, "graphs")  # 定义 graphs 目录路径
#创建目录（如果不存在）
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)  # 创建 tables 目录
}
if (!dir.exists(graphs_dir)) {
  dir.create(graphs_dir)  # 创建 graphs 目录
}
########数据预处理
# 下载 GSE224450 数据
gse224450 <- getGEO("GSE224450", GSEMatrix = TRUE, destdir = getwd())
data224450 <- exprs(gse224450[[1]])
pheno_data <- pData(gse224450[[1]])
#分组为 High_CP 和 Low_CP
Low_CP_group <- data224450[, colnames(data224450) %in% c("GSM7024418", "GSM7024402", "GSM7024424", "GSM7024411", "GSM7024416")]
High_CP_group <- data224450[, colnames(data224450) %in% c("GSM7024387", "GSM7024389", "GSM7024391", "GSM7024392", "GSM7024394", "GSM7024396", "GSM7024398", "GSM7024404", "GSM7024407", "GSM7024413", "GSM7024414", "GSM7024422", "GSM7024385", "GSM7024400", "GSM7024420", "GSM7024409", "GSM7024405")]
#合并数据
combined_data <- cbind(Low_CP_group, High_CP_group)
#########DESeq2 标准化和差异表达分析
#创建分组信息
sample_info <- data.frame(
  row.names = colnames(combined_data),
  Group = rep(c("Low_CP", "High_CP"), c(ncol(Low_CP_group), ncol(High_CP_group)))
)
#确保 Group 是因子类型
sample_info$Group <- factor(sample_info$Group, levels = c("Low_CP", "High_CP"))
#过滤低表达基因
keep <- rowSums(combined_data >= 10) >= 2
combined_data <- combined_data[keep, ]
#创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = combined_data, colData = sample_info, design = ~ Group)
#运行 DESeq2 分析
dds <- DESeq(dds)
#获取差异分析结果
res <- results(dds)
summary(res)
#提取显著差异表达基因
significant_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
###########数据标准化与可视化
#进行 rlog 变换
rld <- rlogTransformation(dds)
#提取标准化后的表达矩阵
exprSet_new <- assay(rld)
#可视化标准化前后的数据分布
cols <- c("blue", "red")
par(mfrow = c(1, 2))
boxplot(combined_data, col = cols, main = "Expression Value (Before Normalization)", las = 2)
boxplot(exprSet_new, col = cols, main = "Expression Value (After Normalization)", las = 2)
#直方图可视化
par(mfrow = c(1, 2))
hist(combined_data, main = "Histogram Before Normalization", col = "grey", border = "white")
hist(exprSet_new, main = "Histogram After Normalization", col = "grey", border = "white")
# ============================================================================
# 数据准备与规范化（已有部分）
# ============================================================================
# rlog 变换
rld <- rlogTransformation(dds)
exprSet_new <- assay(rld)

# 设置颜色向量（LowCP vs HighCP）
cols <- c("blue", "red")

# ============================================================================
# 高分辨率四图保存方案（修复版）- 1200×1200 px, 600 DPI
# ============================================================================

if (!dir.exists("Outcome")) {
  dir.create("Outcome", showWarnings = FALSE)
}

# 设置颜色向量
cols <- c("blue", "red")

cat("Generating publication-quality figures (1200×1200 px, 600 DPI)...\n\n")

# ============================================================================
# 方案：使用更大的物理尺寸（inch），DPI 指定为 600
# ============================================================================

# FigA: Boxplot Before Normalization
png(filename = "Outcome/FigA.png", 
    width = 2400, height = 2400, res = 300)  # 2400 px @ 300 DPI = 8 inch
par(mar = c(5, 5, 3, 2), mgp = c(3, 0.8, 0))
boxplot(combined_data, col = cols, 
        main = "Expression Value (Before Normalization)", 
        las = 2, cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, lwd = 1.5)
dev.off()
cat("✓ FigA.png saved successfully\n")

# FigB: Boxplot After Normalization
png(filename = "Outcome/FigB.png", 
    width = 2400, height = 2400, res = 300)
par(mar = c(5, 5, 3, 2), mgp = c(3, 0.8, 0))
boxplot(exprSet_new, col = cols, 
        main = "Expression Value (After Normalization)", 
        las = 2, cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, lwd = 1.5)
dev.off()
cat("✓ FigB.png saved successfully\n")

# FigC: Histogram Before Normalization
png(filename = "Outcome/FigC.png", 
    width = 2400, height = 2400, res = 300)
par(mar = c(5, 5, 3, 2), mgp = c(3, 0.8, 0))
hist(combined_data, 
     main = "Histogram Before Normalization", 
     col = "grey60", border = "white", 
     xlab = "Expression Value", ylab = "Frequency",
     cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
dev.off()
cat("✓ FigC.png saved successfully\n")

# FigD: Histogram After Normalization
png(filename = "Outcome/FigD.png", 
    width = 2400, height = 2400, res = 300)
par(mar = c(5, 5, 3, 2), mgp = c(3, 0.8, 0))
hist(exprSet_new, 
     main = "Histogram After Normalization", 
     col = "grey60", border = "white", 
     xlab = "Expression Value", ylab = "Frequency",
     cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
dev.off()
cat("✓ FigD.png saved successfully\n")

# ============================================================================
# 验证文件是否成功保存
# ============================================================================

fig_files <- c("Outcome/FigA.png", "Outcome/FigB.png", 
               "Outcome/FigC.png", "Outcome/FigD.png")

cat("\n--- File Verification ---\n")
for (fig in fig_files) {
  if (file.exists(fig)) {
    file_size <- file.info(fig)$size / (1024^2)
    cat(sprintf("✓ %s (%.2f MB)\n", fig, file_size))
  } else {
    cat(sprintf("✗ %s NOT FOUND\n", fig))
  }
}

# ============================================================================
# macOS ARM64 专用方案：纯 R 方案（无需 ImageMagick）
# ============================================================================

cat("Setting up packages...\n")

# 安装必需的包（如果还没有）
required_pkgs <- c("png", "grid")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("✓ All packages ready\n\n")

library(magick)
library(png)

cat("Creating compact four-panel composite figure...\n")

# 读取四张图片
imgA <- image_read("Outcome/FigA.png")
imgB <- image_read("Outcome/FigB.png")
imgC <- image_read("Outcome/FigC.png")
imgD <- image_read("Outcome/FigD.png")

cat("✓ PNG files loaded\n")

# 统一缩放到 1200×1200
imgA <- image_scale(imgA, "1200x1200!")
imgB <- image_scale(imgB, "1200x1200!")
imgC <- image_scale(imgC, "1200x1200!")
imgD <- image_scale(imgD, "1200x1200!")

cat("✓ Images scaled to 1200×1200\n")

# 第一行：A 和 B 并排
row1 <- image_append(c(imgA, imgB), stack = FALSE)

# 第二行：C 和 D 并排
row2 <- image_append(c(imgC, imgD), stack = FALSE)

# 垂直叠放两行（无间距）
composite <- image_append(c(row1, row2), stack = TRUE)

cat("✓ Composite created: ", 
    paste(image_info(composite)$width, "x", image_info(composite)$height), "\n")

# 保存拼接图
image_write(composite, 
            path = "Outcome/Fig_Composite_unlabeled.png", 
            format = "png", density = "600x600")

# ================================================================
# 使用 R 原生方法添加小标注（不改变间距）
# ================================================================

png("Outcome/Fig_Composite_PanelABCD.png", width = 2400, height = 2400, res = 300)

# 绘制拼接好的图
plot.new()
rasterImage(readPNG("Outcome/Fig_Composite_unlabeled.png"), 0, 0, 1, 1)

# 添加小标注（仅在四个角上）
# 总宽度为 2400，总高度为 2400（2×1200）
# A：左上角（第一象限左上）
text(0.025, 0.975, "A", cex = 0.8, font = 2, 
     adj = c(0, 1), col = "black", xpd = NA)

# B：右上角（第二象限右上）
text(0.975, 0.975, "B", cex = 0.8, font = 2, 
     adj = c(1, 1), col = "black", xpd = NA)

# C：左下角（第三象限左下）
text(0.025, 0.025, "C", cex = 0.8, font = 2, 
     adj = c(0, 0), col = "black", xpd = NA)

# D：右下角（第四象限右下）
text(0.975, 0.025, "D", cex = 0.8, font = 2, 
     adj = c(1, 0), col = "black", xpd = NA)

dev.off()

cat("\n✓✓✓ SUCCESS ✓✓✓\n")
cat("File saved: Outcome/Fig_Composite_PanelABCD.png\n")
cat("Layout: Compact (no spacing) with small labels\n")

if (file.exists("Outcome/Fig_Composite_PanelABCD.png")) {
  file_size <- file.info("Outcome/Fig_Composite_PanelABCD.png")$size / (1024^2)
  cat(sprintf("File size: %.2f MB\n", file_size))
  cat("✓ Ready for publication!\n")
}
##########检查异常样本
#################################################################################################
# Compute PCA
pca_res <- prcomp(t(exprSet_new), scale. = TRUE)

# Fix rownames to match expression matrix
rownames(sample_info) <- colnames(combined_data)

# Build PCA dataframe
df_pca <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = sample_info$Group,
  Sample = colnames(combined_data)
)

# Detect outliers
df_pca$Outlier <- ifelse(df_pca$PC1 > 0.5 & df_pca$PC2 > 0.5, "Outlier", "Normal")

# Set group colors (extendable)
group_colors <- c("Low_CP" = "#1f78b4", "High_CP" = "#e31a1c")  # prettier palette
group_shapes <- c("Normal" = 16, "Outlier" = 17)  # shapes for Normal vs Outlier

# Base PCA plot
p <- ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group, shape = Outlier), size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Sample), max.overlaps = 20, size = 3.2)

# Loop to draw group-specific ellipses (with fill)
unique_groups <- unique(df_pca$Group)
for (grp in unique_groups) {
  group_data <- df_pca %>% filter(Group == grp)
  if (nrow(group_data) >= 3) {  # Ellipse needs at least 3 points
    p <- p +
      stat_ellipse(
        data = group_data,
        aes(x = PC1, y = PC2, fill = Group),
        geom = "polygon",
        alpha = 0.2,
        level = 0.95,
        show.legend = FALSE
      ) +
      stat_ellipse(
        data = group_data,
        aes(x = PC1, y = PC2, color = Group),
        geom = "path",
        level = 0.95,
        size = 1.2,
        linetype = "solid"
      )
  }
}

# Final polish
p <- p +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(
    title = "PCA Plot with Group Ellipses and Outlier Highlighting",
    x = "PC1",
    y = "PC2",
    color = "Group",
    shape = "Outlier"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90")
  )

# Save plot
ggsave(file.path(graphs_dir, "PCA_Plot_Enhanced.png"), p, width = 9, height = 7, dpi = 300)

# Optional: print to view
print(p)
##################################################################################################

############
### Differential gene filtering
deg_genes <- rownames(subset(res, padj < 0.01 & abs(log2FoldChange) > 2))
write.csv(deg_genes, file.path(tables_dir, "DEGS_Results_Significant_genes.csv"))

# Upregulated & Downregulated
upregulated_genes <- rownames(subset(res, padj < 0.01 & log2FoldChange > 2))
downregulated_genes <- rownames(subset(res, padj < 0.01 & log2FoldChange < -2))

cat("上调基因数量:", length(upregulated_genes), "\n")
cat("下调基因数量:", length(downregulated_genes), "\n")

# If not already done
volcano_data <- as.data.frame(res) %>%
  mutate(
    regulation = case_when(
      padj < 0.01 & log2FoldChange > 2 ~ "Upregulated",
      padj < 0.01 & log2FoldChange < -2 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    gene_name = rownames(res)
  )

# Select top labeled genes for volcano plot
up_labeled <- volcano_data %>%
  filter(regulation == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

down_labeled <- volcano_data %>%
  filter(regulation == "Downregulated") %>%
  arrange(log2FoldChange) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgray") +
  geom_text_repel(data = up_labeled, aes(label = gene_name),
                  size = 3, color = "darkred", fontface = "italic") +
  geom_text_repel(data = down_labeled, aes(label = gene_name),
                  size = 3, color = "darkblue", fontface = "italic") +
  scale_color_manual(values = c("Upregulated" = "#E64B35FF",
                                "Downregulated" = "#4DBBD5FF",
                                "Not Significant" = "gray70")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~p~value),
    color = "Regulation"
  )

# Save it
ggsave("Volcano_Plot_DEGs.pdf", volcano_plot, width = 8, height = 6, device = cairo_pdf)

### Heatmap of DEGs
# Get expression of significant genes
expr_deg <- exprSet_new[deg_genes, ]

# Z-score normalization per gene (row)
expr_scaled <- t(scale(t(expr_deg)))

# Sample annotation
sample_info <- as.data.frame(colData(dds))  # Use the same `dds` object
annotation_col <- sample_info["Group"]
ann_colors <- list(Group = c("High_CP" = "#E64B35FF", "Low_CP" = "#4DBBD5FF"))

# Optional: load library if not already loaded
library(pheatmap)

png("Heatmap_DEGs.png", width = 1000, height = 800, res = 150)

pheatmap(expr_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize = 10,
         main = "Heatmap of Top Differentially Expressed Genes",
         legend = TRUE,
         border_color = NA)

dev.off()


# ============================================================================
# 第一步：生成火山图和热图（高分辨率）
# ============================================================================

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(png)

cat("Generating volcano plot and heatmap...\n\n")

# ========== 火山图 ==========
cat("Creating volcano plot...\n")

volcano_data <- as.data.frame(res) %>%
  mutate(
    regulation = case_when(
      padj < 0.01 & log2FoldChange > 2 ~ "Upregulated",
      padj < 0.01 & log2FoldChange < -2 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    gene_name = rownames(res)
  )

# Select top labeled genes for volcano plot
up_labeled <- volcano_data %>%
  filter(regulation == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

down_labeled <- volcano_data %>%
  filter(regulation == "Downregulated") %>%
  arrange(log2FoldChange) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgray") +
  geom_text_repel(data = up_labeled, aes(label = gene_name),
                  size = 3, color = "darkred", fontface = "italic") +
  geom_text_repel(data = down_labeled, aes(label = gene_name),
                  size = 3, color = "darkblue", fontface = "italic") +
  scale_color_manual(values = c("Upregulated" = "#E64B35FF",
                                "Downregulated" = "#4DBBD5FF",
                                "Not Significant" = "gray70")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Volcano Plot",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~p~value),
    color = "Regulation"
  )

# 保存火山图为高分辨率 PNG
ggsave("Outcome/Fig_VolcanoPlot.png", volcano_plot, 
       width = 10, height = 8, dpi = 300, device = "png")

cat("✓ Volcano plot saved: Outcome/Fig_VolcanoPlot.png\n")

# ========== 热图 ==========
cat("Creating heatmap...\n")

# Get expression of significant genes
expr_deg <- exprSet_new[deg_genes, ]

# Z-score normalization per gene (row)
expr_scaled <- t(scale(t(expr_deg)))

# Sample annotation
sample_info <- as.data.frame(colData(dds))
annotation_col <- sample_info["Group"]
ann_colors <- list(Group = c("High_CP" = "#E64B35FF", "Low_CP" = "#4DBBD5FF"))

# 保存热图为高分辨率 PNG
png("Outcome/Fig_Heatmap.png", width = 2400, height = 2400, res = 300)

pheatmap(expr_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize = 10,
         main = "Heatmap of Top DEGs",
         legend = TRUE,
         border_color = NA)

dev.off()

cat("✓ Heatmap saved: Outcome/Fig_Heatmap.png\n")

# ============================================================================
# 第二步：将火山图和热图拼接成 1×2 面板图
# ============================================================================

cat("\nCompositing volcano and heatmap into two-panel figure...\n")

library(magick)

tryCatch({
  # 读取两张图片
  imgVolcano <- image_read("Outcome/Fig_VolcanoPlot.png")
  imgHeatmap <- image_read("Outcome/Fig_Heatmap.png")
  
  cat("✓ Images loaded\n")
  
  # 统一缩放到 1200×1200
  imgVolcano <- image_scale(imgVolcano, "1200x1200!")
  imgHeatmap <- image_scale(imgHeatmap, "1200x1200!")
  
  cat("✓ Images scaled to 1200×1200\n")
  
  # 并排拼接（左：火山图，右：热图）
  composite <- image_append(c(imgVolcano, imgHeatmap), stack = FALSE)
  
  cat("✓ Composite created: ", 
      paste(image_info(composite)$width, "x", image_info(composite)$height), "\n")
  
  # 保存拼接图
  image_write(composite, 
              path = "Outcome/Fig_Composite_unlabeled_VolHeat.png", 
              format = "png", density = "600x600")
  
  # ================================================================
  # 使用 R 原生方法添加小标注
  # ================================================================
  
  png("Outcome/Fig_Composite_VolcanoHeatmap.png", 
      width = 2400, height = 1200, res = 300)
  
  # 绘制拼接好的图
  plot.new()
  rasterImage(readPNG("Outcome/Fig_Composite_unlabeled_VolHeat.png"), 0, 0, 1, 1)
  
  # 添加小标注
  # A：左侧（火山图）
  text(0.025, 0.975, "A", cex = 0.8, font = 2, 
       adj = c(0, 1), col = "black", xpd = NA)
  
  # B：右侧（热图）
  text(0.975, 0.975, "B", cex = 0.8, font = 2, 
       adj = c(1, 1), col = "black", xpd = NA)
  
  dev.off()
  
  cat("\n✓✓✓ SUCCESS ✓✓✓\n")
  cat("Two-panel figure saved: Outcome/Fig_Composite_VolcanoHeatmap.png\n")
  
  if (file.exists("Outcome/Fig_Composite_VolcanoHeatmap.png")) {
    file_size <- file.info("Outcome/Fig_Composite_VolcanoHeatmap.png")$size / (1024^2)
    cat(sprintf("File size: %.2f MB\n", file_size))
    cat("Layout: Volcano Plot (A) + Heatmap (B), side by side\n")
    cat("Resolution: 2400×1200 px @ 300 DPI\n")
    cat("✓ Ready for publication!\n")
  }
  
}, error = function(e) {
  cat("\n✗✗✗ ERROR ✗✗✗\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\nNote: This error may occur if magick annotation fails.\n")
  cat("Alternative: Use the unlabeled composite at:\n")
  cat("Outcome/Fig_Composite_unlabeled_VolHeat.png\n")
})



########################################################################################
# 提取基因列表并转换为 ENTREZID
######## GO 和 KEGG 富集分析（分别对上调和下调基因）
# 提取上调和下调基因的 ENTREZID
upregulated_entrez <- bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
downregulated_entrez <- bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 上调和下调基因的 ENTREZID 向量
upregulated_entrez_vector <- upregulated_entrez$ENTREZID
downregulated_entrez_vector <- downregulated_entrez$ENTREZID

# 运行 enrichGO 并合并结果 --------------------------------------------------
run_enrichGO <- function(ontology) {
  enrichGO(
    gene = upregulated_entrez_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}

# 获取 BP/CC/MF 富集结果
GO_BP_up <- run_enrichGO("BP")
GO_CC_up <- run_enrichGO("CC")
GO_MF_up <- run_enrichGO("MF")

# 添加 ONTOLOGY 列并转换 GeneRatio 为数值型
process_go <- function(go_obj, ontology) {
  result <- go_obj@result
  result$ONTOLOGY <- ontology
  result$GeneRatio <- apply(result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))  # 将 "10/100" 转换为 0.1
  })
  result
}

# 合并结果（每个类别取前10）
GO_combined <- rbind(
  process_go(GO_BP_up, "BP") %>% arrange(desc(Count)) %>% head(10),
  process_go(GO_CC_up, "CC") %>% arrange(desc(Count)) %>% head(10),
  process_go(GO_MF_up, "MF") %>% arrange(desc(Count)) %>% head(10)
)

# 确保 Description 顺序
GO_combined$Description <- factor(
  GO_combined$Description,
  levels = rev(unique(GO_combined$Description))
)

#保存GO_up的数据到tables文件夹
write.csv(GO_BP_up, file.path(tables_dir, "GO_BP_Upregulated_Results.csv"))
write.csv(GO_CC_up, file.path(tables_dir, "GO_CC_Upregulated_Results.csv"))
write.csv(GO_MF_up, file.path(tables_dir, "GO_MF_Upregulated_Results.csv"))


# 绘图 ------------------------------------------------------------------------
ontology_labels <- c(
  BP = "Biological Process",
  CC = "Cellular Component",
  MF = "Molecular Function"
)

p <- ggplot(GO_combined, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  # 分区块背景色（禁止继承全局 aes）
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ONTOLOGY),
    data = GO_combined %>% distinct(ONTOLOGY),  # 独立数据框
    alpha = 0.1,
    inherit.aes = FALSE,  # 关键修正：禁止继承主图 aes
    show.legend = FALSE
  ) +
  # 绘制点图
  geom_point(alpha = 0.8) +
  # 分面与主题设置
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free", 
             labeller = labeller(ONTOLOGY = ontology_labels)) +
  scale_fill_manual(values = c(BP = "#FFE4B5", CC = "#B0E0E6", MF = "#98FB98")) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(title = "GO Enrichment (BP, CC, MF) - Upregulated Genes",
       x = "Gene Ratio", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 14, face = "bold", hjust = 0),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.major.y = element_line(color = "grey90")
  )

# 保存图形
ggsave(file.path(graphs_dir, "GO_Combined_Dotplot_Fixed.png"), 
       p, width = 14, height = 18, dpi = 300)


run_enrichGO_down <- function(ontology) {
  enrichGO(
    gene = downregulated_entrez_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}

# 获取 BP/CC/MF 富集结果
GO_BP_down <- run_enrichGO_down("BP")
GO_CC_down <- run_enrichGO_down("CC")
GO_MF_down <- run_enrichGO_down("MF")

# 处理结果（添加 ONTOLOGY 列并转换 GeneRatio）
process_go_down <- function(go_obj, ontology) {
  result <- go_obj@result
  result$ONTOLOGY <- ontology
  result$GeneRatio <- apply(result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))  # 将 "10/100" 转换为 0.1
  })
  result
}

# 合并结果（每个类别取前10）
GO_combined_down <- rbind(
  process_go_down(GO_BP_down, "BP") %>% arrange(desc(Count)) %>% head(10),
  process_go_down(GO_CC_down, "CC") %>% arrange(desc(Count)) %>% head(10),
  process_go_down(GO_MF_down, "MF") %>% arrange(desc(Count)) %>% head(10)
)

# 确保 Description 顺序
GO_combined_down$Description <- factor(
  GO_combined_down$Description,
  levels = rev(unique(GO_combined_down$Description))
)

# 绘图参数设置
ontology_labels <- c(
  BP = "Biological Process",
  CC = "Cellular Component",
  MF = "Molecular Function"
)

# 绘制分区块点图
p_down <- ggplot(GO_combined_down, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ONTOLOGY),
    data = GO_combined_down %>% distinct(ONTOLOGY),
    alpha = 0.1,
    inherit.aes = FALSE,  # 禁止继承主图 aes
    show.legend = FALSE
  ) +
  geom_point(alpha = 0.8) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free", 
             labeller = labeller(ONTOLOGY = ontology_labels)) +
  scale_fill_manual(values = c(BP = "#FFE4B5", CC = "#B0E0E6", MF = "#98FB98")) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(title = "GO Enrichment (BP, CC, MF) - Downregulated Genes",
       x = "Gene Ratio", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 14, face = "bold", hjust = 0),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.major.y = element_line(color = "grey90")
  )

# 保存图形（文件名区分上调/下调）
ggsave(file.path(graphs_dir, "GO_Combined_Dotplot_Downregulated.png"), 
       p_down, width = 14, height = 18, dpi = 300)


#保存GO_down的数据到tables文件夹
write.csv(GO_BP_down, file.path(tables_dir, "GO_BP_Downregulated_Results.csv"))
write.csv(GO_CC_down, file.path(tables_dir, "GO_CC_Downregulated_Results.csv"))
write.csv(GO_MF_down, file.path(tables_dir, "GO_MF_Downregulated_Results.csv"))
#########################################################################################################
# 运行 KEGG 富集分析
KEGG_up <- enrichKEGG(
  gene = upregulated_entrez_vector,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# 保存结果
write.csv(KEGG_up, file.path(tables_dir, "KEGG_Upregulated_Results.csv"))

# 转换 GeneRatio 为数值型
KEGG_up@result$GeneRatio <- apply(KEGG_up@result, 1, function(x) {
  eval(parse(text = x["GeneRatio"]))  # 将 "10/100" 转换为 0.1
})

# 提取前10个通路（按 GeneRatio 降序）
kegg_data <- KEGG_up@result %>%
  arrange(desc(GeneRatio)) %>%
  head(10)

# 固定通路描述的顺序（按 GeneRatio 降序）
kegg_data$Description <- factor(
  kegg_data$Description,
  levels = rev(unique(kegg_data$Description))  # 反转顺序以匹配降序排列
)

# 绘制长形点图（横向布局）
p <- ggplot(kegg_data, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count", range = c(3, 8)) +  # 调整点的大小范围
  labs(
    title = "KEGG Pathway Enrichment - Upregulated Genes",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90"),  # 添加横向网格线
    legend.position = "right"
  )

# 保存高清图
ggsave(
  file.path(graphs_dir, "KEGG_Upregulated_Dotplot_Enhanced.png"),
  p,
  width = 12,
  height = 8,
  dpi = 300
)

# 可选：横向条形图（按 GeneRatio 排序）
p_bar <- ggplot(kegg_data, aes(x = GeneRatio, y = Description, fill = p.adjust)) +
  geom_col(width = 0.8) +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(
    title = "KEGG Pathway Enrichment - Upregulated Genes",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  file.path(graphs_dir, "KEGG_Upregulated_Barplot_Enhanced.png"),
  p_bar,
  width = 12,
  height = 8,
  dpi = 300
)

# KEGG 通路可视化（原代码保持不变）
if (!is.null(KEGG_up)) {
  pathways <- KEGG_up$ID[1:min(5, nrow(KEGG_up))]
  for (pathway_id in pathways) {
    pathview(
      gene.data = upregulated_entrez_vector,
      pathway.id = pathway_id,
      species = "hsa",
      gene.idtype = "entrez",
      limit = list(gene = max(abs(res$log2FoldChange))),
      kegg.native = TRUE,
      same.layer = TRUE,
      out.suffix = "upregulated"
    )
    # 自动保存为PNG，无需手动保存
  }
}
##########################################################################

# 运行 KEGG 富集分析（下调基因）
KEGG_down <- enrichKEGG(
  gene = downregulated_entrez_vector,  # 需确保此变量已定义
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# 保存结果
write.csv(KEGG_down, file.path(tables_dir, "KEGG_Downregulated_Results.csv"))

# 转换 GeneRatio 为数值型
if (!is.null(KEGG_down)) {
  KEGG_down@result$GeneRatio <- apply(KEGG_down@result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))  # 将 "10/100" 转换为 0.1
  })
}

# 提取前10个通路（按 GeneRatio 降序）
kegg_data_down <- KEGG_down@result %>%
  arrange(desc(GeneRatio)) %>%
  head(10)

# 固定通路描述的顺序（按 GeneRatio 降序）
kegg_data_down$Description <- factor(
  kegg_data_down$Description,
  levels = rev(unique(kegg_data_down$Description))  # 反转顺序以匹配降序排列
)

# 绘制长形点图（横向布局）
p_down <- ggplot(kegg_data_down, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
  labs(
    title = "KEGG Pathway Enrichment - Downregulated Genes",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "right"
  )

# 保存高清图
ggsave(
  file.path(graphs_dir, "KEGG_Downregulated_Dotplot_Enhanced.png"),
  p_down,
  width = 12,
  height = 8,
  dpi = 300
)

# 横向条形图（按 GeneRatio 排序）
p_bar_down <- ggplot(kegg_data_down, aes(x = GeneRatio, y = Description, fill = p.adjust)) +
  geom_col(width = 0.8) +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(
    title = "KEGG Pathway Enrichment - Downregulated Genes",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  file.path(graphs_dir, "KEGG_Downregulated_Barplot_Enhanced.png"),
  p_bar_down,
  width = 12,
  height = 8,
  dpi = 300
)

# KEGG 通路可视化（映射下调基因 log2FoldChange）
if (!is.null(KEGG_down)) {
  pathways_down <- KEGG_down$ID[1:min(5, nrow(KEGG_down))]
  for (pathway_id in pathways_down) {
    pathview(
      gene.data = res[rownames(res) %in% downregulated_entrez_vector, "log2FoldChange"],  # 提取下调基因的 log2FC
      pathway.id = pathway_id,
      species = "hsa",
      gene.idtype = "entrez",
      limit = list(gene = max(abs(res$log2FoldChange))),
      kegg.native = TRUE,
      same.layer = TRUE,
      out.suffix = "downregulated",  # 文件名后缀
    )
  }
}


# ============================================================================
# GO 和 KEGG 富集分析合并可视化方案
# ============================================================================

library(magick)
library(png)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

cat("Creating GO and KEGG enrichment figures...\n\n")

# ========== 步骤1：生成 GO 富集图（上调基因）==========

cat("1. Generating GO enrichment plot (Upregulated genes)...\n")

# 提取上调基因的 ENTREZID
upregulated_entrez <- bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
upregulated_entrez_vector <- upregulated_entrez$ENTREZID

# 运行 enrichGO
run_enrichGO <- function(ontology) {
  enrichGO(
    gene = upregulated_entrez_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}

GO_BP_up <- run_enrichGO("BP")
GO_CC_up <- run_enrichGO("CC")
GO_MF_up <- run_enrichGO("MF")

# 处理结果
process_go <- function(go_obj, ontology) {
  result <- go_obj@result
  result$ONTOLOGY <- ontology
  result$GeneRatio <- apply(result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))
  })
  result
}

GO_combined <- rbind(
  process_go(GO_BP_up, "BP") %>% arrange(desc(Count)) %>% head(10),
  process_go(GO_CC_up, "CC") %>% arrange(desc(Count)) %>% head(10),
  process_go(GO_MF_up, "MF") %>% arrange(desc(Count)) %>% head(10)
)

GO_combined$Description <- factor(
  GO_combined$Description,
  levels = rev(unique(GO_combined$Description))
)

ontology_labels <- c(
  BP = "Biological Process",
  CC = "Cellular Component",
  MF = "Molecular Function"
)

p_GO_up <- ggplot(GO_combined, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ONTOLOGY),
    data = GO_combined %>% distinct(ONTOLOGY),
    alpha = 0.1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(alpha = 0.8) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free", 
             labeller = labeller(ONTOLOGY = ontology_labels)) +
  scale_fill_manual(values = c(BP = "#FFE4B5", CC = "#B0E0E6", MF = "#98FB98")) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
  labs(title = "GO Enrichment - Upregulated",
       x = "Gene Ratio", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold", hjust = 0),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave("Outcome/Fig_GO_Upregulated.png", p_GO_up, 
       width = 12, height = 10, dpi = 300, device = "png")

cat("   ✓ GO Upregulated saved\n")

# ========== 步骤2：生成 GO 富集图（下调基因）==========

cat("2. Generating GO enrichment plot (Downregulated genes)...\n")

downregulated_entrez <- bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
downregulated_entrez_vector <- downregulated_entrez$ENTREZID

run_enrichGO_down <- function(ontology) {
  enrichGO(
    gene = downregulated_entrez_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
}

GO_BP_down <- run_enrichGO_down("BP")
GO_CC_down <- run_enrichGO_down("CC")
GO_MF_down <- run_enrichGO_down("MF")

process_go_down <- function(go_obj, ontology) {
  result <- go_obj@result
  result$ONTOLOGY <- ontology
  result$GeneRatio <- apply(result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))
  })
  result
}

GO_combined_down <- rbind(
  process_go_down(GO_BP_down, "BP") %>% arrange(desc(Count)) %>% head(10),
  process_go_down(GO_CC_down, "CC") %>% arrange(desc(Count)) %>% head(10),
  process_go_down(GO_MF_down, "MF") %>% arrange(desc(Count)) %>% head(10)
)

GO_combined_down$Description <- factor(
  GO_combined_down$Description,
  levels = rev(unique(GO_combined_down$Description))
)

p_GO_down <- ggplot(GO_combined_down, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ONTOLOGY),
    data = GO_combined_down %>% distinct(ONTOLOGY),
    alpha = 0.1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(alpha = 0.8) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free", 
             labeller = labeller(ONTOLOGY = ontology_labels)) +
  scale_fill_manual(values = c(BP = "#FFE4B5", CC = "#B0E0E6", MF = "#98FB98")) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
  labs(title = "GO Enrichment - Downregulated",
       x = "Gene Ratio", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold", hjust = 0),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave("Outcome/Fig_GO_Downregulated.png", p_GO_down, 
       width = 12, height = 10, dpi = 300, device = "png")

cat("   ✓ GO Downregulated saved\n")

# ========== 步骤3：生成 KEGG 富集图（上调基因）==========

cat("3. Generating KEGG enrichment plot (Upregulated genes)...\n")

KEGG_up <- enrichKEGG(
  gene = upregulated_entrez_vector,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

KEGG_up@result$GeneRatio <- apply(KEGG_up@result, 1, function(x) {
  eval(parse(text = x["GeneRatio"]))
})

kegg_data <- KEGG_up@result %>%
  arrange(desc(GeneRatio)) %>%
  head(10)

kegg_data$Description <- factor(
  kegg_data$Description,
  levels = rev(unique(kegg_data$Description))
)

p_KEGG_up <- ggplot(kegg_data, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
  scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
  labs(
    title = "KEGG Pathways - Upregulated",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave("Outcome/Fig_KEGG_Upregulated.png", p_KEGG_up, 
       width = 12, height = 8, dpi = 300, device = "png")

cat("   ✓ KEGG Upregulated saved\n")

# ========== 步骤4：生成 KEGG 富集图（下调基因）==========

cat("4. Generating KEGG enrichment plot (Downregulated genes)...\n")

KEGG_down <- enrichKEGG(
  gene = downregulated_entrez_vector,
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

KEGG_down@result$GeneRatio <- apply(KEGG_down@result, 1, function(x) {
  eval(parse(text = x["GeneRatio"]))
})

kegg_data_down <- KEGG_down@result %>%
  arrange(desc(GeneRatio)) %>%
  head(10)

kegg_data_down$Description <- factor(
  kegg_data_down$Description,
  levels = rev(unique(kegg_data_down$Description))
)

p_KEGG_down <- ggplot(kegg_data_down, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
  scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
  labs(
    title = "KEGG Pathways - Downregulated",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom"
  )

ggsave("Outcome/Fig_KEGG_Downregulated.png", p_KEGG_down, 
       width = 12, height = 8, dpi = 300, device = "png")

cat("   ✓ KEGG Downregulated saved\n")

# ============================================================================
# 步骤5：拼接成 2×2 四面板图（GO Upregulated, GO Downregulated, KEGG Upregulated, KEGG Downregulated）
# ============================================================================

cat("\n5. Compositing into 2×2 four-panel figure...\n")

tryCatch({
  # 读取四张图片
  imgGO_up <- image_read("Outcome/Fig_GO_Upregulated.png")
  imgGO_down <- image_read("Outcome/Fig_GO_Downregulated.png")
  imgKEGG_up <- image_read("Outcome/Fig_KEGG_Upregulated.png")
  imgKEGG_down <- image_read("Outcome/Fig_KEGG_Downregulated.png")
  
  cat("   ✓ All images loaded\n")
  
  # 统一缩放到 1200×1200
  imgGO_up <- image_scale(imgGO_up, "1200x1200!")
  imgGO_down <- image_scale(imgGO_down, "1200x1200!")
  imgKEGG_up <- image_scale(imgKEGG_up, "1200x1200!")
  imgKEGG_down <- image_scale(imgKEGG_down, "1200x1200!")
  
  cat("   ✓ Images scaled to 1200×1200\n")
  
  # 第一行：GO上调 和 GO下调
  row1 <- image_append(c(imgGO_up, imgGO_down), stack = FALSE)
  
  # 第二行：KEGG上调 和 KEGG下调
  row2 <- image_append(c(imgKEGG_up, imgKEGG_down), stack = FALSE)
  
  # 垂直叠放两行
  composite <- image_append(c(row1, row2), stack = TRUE)
  
  cat("   ✓ Composite created: ", 
      paste(image_info(composite)$width, "x", image_info(composite)$height), "\n")
  
  # 保存拼接图
  image_write(composite, 
              path = "Outcome/Fig_Composite_unlabeled_GO_KEGG.png", 
              format = "png", density = "600x600")
  
  # 使用 R 原生方法添加小标注
  png("Outcome/Fig_Composite_GO_KEGG.png", width = 2400, height = 2400, res = 300)
  
  plot.new()
  rasterImage(readPNG("Outcome/Fig_Composite_unlabeled_GO_KEGG.png"), 0, 0, 1, 1)
  
  # 添加小标注在四个角
  text(0.025, 0.975, "A", cex = 0.8, font = 2, 
       adj = c(0, 1), col = "black", xpd = NA)
  text(0.975, 0.975, "B", cex = 0.8, font = 2, 
       adj = c(1, 1), col = "black", xpd = NA)
  text(0.025, 0.025, "C", cex = 0.8, font = 2, 
       adj = c(0, 0), col = "black", xpd = NA)
  text(0.975, 0.025, "D", cex = 0.8, font = 2, 
       adj = c(1, 0), col = "black", xpd = NA)
  
  dev.off()
  
  cat("\n✓✓✓ SUCCESS ✓✓✓\n")
  cat("Four-panel figure saved: Outcome/Fig_Composite_GO_KEGG.png\n")
  cat("Layout:\n")
  cat("  A (Top-Left): GO Enrichment - Upregulated Genes\n")
  cat("  B (Top-Right): GO Enrichment - Downregulated Genes\n")
  cat("  C (Bottom-Left): KEGG Pathways - Upregulated Genes\n")
  cat("  D (Bottom-Right): KEGG Pathways - Downregulated Genes\n")
  
  if (file.exists("Outcome/Fig_Composite_GO_KEGG.png")) {
    file_size <- file.info("Outcome/Fig_Composite_GO_KEGG.png")$size / (1024^2)
    cat(sprintf("\nFile size: %.2f MB\n", file_size))
    cat("Resolution: 2400×2400 px @ 300 DPI\n")
    cat("✓ Ready for publication!\n")
  }
  
}, error = function(e) {
  cat("\n✗ Error during composite creation:\n")
  cat(conditionMessage(e), "\n")
  cat("\nUsing unlabeled composite at:\n")
  cat("Outcome/Fig_Composite_unlabeled_GO_KEGG.png\n")
})


#############免疫浸润分析########################################################

# 读取数据与运行CIBERSORT
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, exprSet_new, perm = 0, QN = TRUE)
results <- as.data.frame(results)

# 筛选低可信样本
low_confidence_samples <- results[results$Correlation < 0.3 | results$RMSE > 1.0, ]
results_filtered <- results[!rownames(results) %in% rownames(low_confidence_samples), ]

# 保存结果
write.csv(results, "CIBERSORT_Results.csv")
write.csv(results_filtered, "CIBERSORT_Results_Filtered.csv")

immune_data <- results[, 1:22]

# 把行名变成新的一列 sample
immune_data <- immune_data %>%
  rownames_to_column(var = "sample")

# 转成长格式
immune_long <- immune_data %>%
  pivot_longer(
    cols = -sample,
    names_to = "cell_type",
    values_to = "proportion"
  )

head(immune_long)



# 从 Low_CP_group 和 High_CP_group 提取样本名
low_cp_samples <- colnames(Low_CP_group)
high_cp_samples <- colnames(High_CP_group)

# 将immune profiling转为long format（假设你已有immune_wide数据框）
immune_long <- immune_data %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "proportion") %>%
  mutate(group = case_when(
    sample %in% low_cp_samples ~ "Low_CP",
    sample %in% high_cp_samples ~ "High_CP",
    TRUE ~ "Unknown"
  ))

### 2. 发表级别箱线图（带显著性星号） -----

p_box <- ggplot(immune_long, aes(x = cell_type, y = proportion, fill = group)) +
  geom_boxplot(
    outlier.shape = 21, outlier.size = 1.8, lwd = 0.5, alpha = 0.8, width=0.6
  ) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Low_CP" = "#4E79A7", "High_CP" = "#E15759")) +
  labs(
    x = "Immune Cell Type", y = "Proportion",
    title = "Immune Profiling Comparison", fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
    axis.text.y = element_text(face="bold"),
    axis.title = element_text(face="bold", size=14),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face="bold", size=16)
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif", hide.ns = TRUE, size = 5
  )

# 保存箱线图
ggsave("Immune_Boxplot_Publication.png", p_box, width = 10, height = 6, dpi = 300)

### 3. 热图部分 -----

# 确保数据格式匹配
valid_samples <- intersect(rownames(results_filtered), c(low_cp_samples, high_cp_samples))
group <- data.frame(
  Group = ifelse(rownames(results_filtered) %in% low_cp_samples, "Low_CP", "High_CP"),
  row.names = rownames(results_filtered)
)

annotation_colors <- list(Group = c(Low_CP = "#4E79A7", High_CP = "#E15759"))

# 热图绘制
pheatmap(
  t(results_filtered[, 1:22]),
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = group,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 12,
  main = "Immune Cell Composition Heatmap"
)

# 保存高清热图
png("Immune_Heatmap_Publication.png", width = 3000, height = 2000, res = 300)
pheatmap(
  t(results_filtered[, 1:22]),
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = group,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  fontsize_row = 12,
  main = "Immune Cell Composition Heatmap"
)
dev.off()

### 4. 堆叠柱状图（Composition Plot） -----

# 数据准备
plot_data <- results_filtered[, 1:22] %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Group = ifelse(Sample %in% low_cp_samples, "Low_CP", "High_CP"))

# 颜色可以用viridis更美观
n_colors <- length(unique(plot_data$CellType))
color_palette <- viridis(n_colors)

p_stack <- ggplot(plot_data, aes(Sample, Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Immune Cell Composition (Stacked Barplot)", x=NULL, y="Proportion", fill="Cell Type") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
    plot.title = element_text(hjust=0.5, face="bold", size=16),
    legend.position = "right"
  )

# 保存堆叠图
ggsave("Immune_StackedBar_Publication.png", p_stack, width = 12, height = 8, dpi = 300)

#####################################################################################
##################################################
#draw venn plot based on
#upregulated pathway：
#GO BP: KLRK1, SH2D1A, BCL2, CCL19, CCL21, VCAM1, CD5，CD3E, CD247, LCK, IL2RG, HLA-DPB1, MS4A1
#KEGG: CD3D, CD3G, CD19, CD34, HLA-DQB1，CD3E, CD247, LCK, IL2RG, HLA-DPB1, MS4A1
#Common: CD3E, CD247, LCK, IL2RG, HLA-DPB1, MS4A1
#downregulated pathway：
#GO BP:CXCL8, IL1B, C5AR1, TREM1, CXCR2, CSF3R, S100A9, S100A8, NLRP3, IL1A, THBD, ITGAX
#KEGG:CXCL8, IL1B, C5AR1, TREM1, CXCR2, CSF3R, S100A9, S100A8, NLRP3, IL1A, THBD, ITGAX
#Common: CXCL8, IL1B, C5AR1, TREM1, CXCR2, CSF3R, S100A9, S100A8, NLRP3, IL1A, THBD, ITGAX

# Required Libraries
library(VennDiagram)
library(grid)
library(gridExtra)
library(ggplot2)
library(gtable)

# Output Directory Setup
venn_dir <- file.path(graphs_dir, "Venn_Diagrams")
if (!dir.exists(venn_dir)) dir.create(venn_dir)

# Gene Sets Definition
upregulated_GO_BP <- c("KLRK1", "SH2D1A", "BCL2", "CCL19", "CCL21", "VCAM1", "CD5", "CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
upregulated_KEGG <- c("CD3D", "CD3G", "CD19", "CD34", "HLA-DQB1", "CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
downregulated_GO_BP <- c("CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX")
downregulated_KEGG <- c("CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX")

# Intersections
up_intersect <- intersect(upregulated_GO_BP, upregulated_KEGG)
down_intersect <- intersect(downregulated_GO_BP, downregulated_KEGG)

# Venn Diagram Generator
create_venn_plot <- function(set1, set2, cat_names, title, fill_colors) {
  venn <- venn.diagram(
    x = list(set1, set2),
    category.names = cat_names,
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    height = 1000,
    width = 1000,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    col = "black",
    fill = fill_colors,
    alpha = 0.7,
    cex = 1.5,
    fontfamily = "sans",
    cat.cex = 1.5,
    cat.fontfamily = "sans",
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05),
    margin = 0.1
  )
  return(gTree(children = venn))
}

# Generate Plots
venn_up <- create_venn_plot(upregulated_GO_BP, upregulated_KEGG, c("GO_BP", "KEGG"), NULL, c("#4E79A7", "#F28E2B"))
venn_down <- create_venn_plot(downregulated_GO_BP, downregulated_KEGG, c("GO_BP", "KEGG"), NULL, c("#59A14F", "#E15759"))

# Combine with Title and Footer
footer_text <- paste(
  "Upregulated Intersection:", paste(up_intersect, collapse = ", "), "\n",
  "Downregulated Intersection:", paste(down_intersect, collapse = ", ")
)
footer_grob <- textGrob(footer_text, gp = gpar(fontsize = 12), just = "centre")
main_title <- textGrob("Pathway Gene Overlap Analysis", gp = gpar(fontsize = 20, fontface = "bold"))

# Layout and Save Composite Figure
layout <- grid.arrange(
  main_title,
  arrangeGrob(venn_up, venn_down, ncol = 2),
  footer_grob,
  nrow = 3,
  heights = c(0.1, 0.8, 0.1)
)

png(file.path(venn_dir, "Combined_Pathway_Venn.png"), width = 2000, height = 1200, res = 300)
grid.draw(layout)
dev.off()

# Save Individual Venn Diagrams
png(file.path(venn_dir, "Upregulated_Pathway_Venn.png"), width = 1000, height = 1000, res = 300)
grid.draw(venn_up)
dev.off()

png(file.path(venn_dir, "Downregulated_Pathway_Venn.png"), width = 1000, height = 1000, res = 300)
grid.draw(venn_down)
dev.off()

# Save Intersection Gene Lists
write.csv(data.frame(Gene = up_intersect), file.path(tables_dir, "Upregulated_Intersection_Genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = down_intersect), file.path(tables_dir, "Downregulated_Intersection_Genes.csv"), row.names = FALSE)

# Console Output
cat("Venn diagrams saved to:", venn_dir, "\n")
cat("Upregulated Intersection Genes:", paste(up_intersect, collapse = ", "), "\n")
cat("Downregulated Intersection Genes:", paste(down_intersect, collapse = ", "), "\n")

########################################
###从新筛选的hub genes从新进行分析
# 加载所需 R 包
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(randomForest)
library(dplyr)
library(tidyr)
library(corrplot)

# 定义 hub 基因
upregulated_genes <- c("CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
downregulated_genes <- c("CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX")

# 读取免疫细胞比例 (CIBERSORT 结果)
immune_data <- results[, 1:22]  # 假设 results 已经包含免疫细胞比例

# 读取基因表达数据
gene_expr <- exprSet_new[c(upregulated_genes, downregulated_genes), ]  # 筛选 hub 基因

# 计算相关性（Pearson）
cor_matrix <- cor(t(gene_expr), immune_data, method = "pearson")

# 绘制相关性热图
png("correlation_heatmap.png", width = 1200, height = 1000)
pheatmap(cor_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "Gene-Immune Cell Correlation Heatmap")
dev.off()

# 转换相关性数据用于柱状图
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Gene", "Immune_Cell", "Correlation")
write.csv(cor_df, file.path(tables_dir, "Hub_Genes_Immune_Correlation.csv"), row.names = FALSE)
# 绘制柱状图
png("correlation_barplot.png", width = 1200, height = 800)
ggplot(cor_df, aes(x = reorder(Immune_Cell, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Correlation Between Hub Genes and Immune Cells",
       x = "Immune Cell Type",
       y = "Correlation") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
dev.off()

# ---- 随机森林分析 ----
# 构建随机森林模型预测基因表达
X <- immune_data  # 免疫细胞比例
importance_list <- list()

for (gene in c(upregulated_genes, downregulated_genes)) {
  y <- gene_expr[gene, ]  # 目标基因表达
  
  rf_model <- randomForest(x = X, y = y, importance = TRUE)
  
  # 获取特征重要性
  importance_scores <- importance(rf_model)
  importance_df <- data.frame(Immune_Cell = rownames(importance_scores), Importance = importance_scores[, "IncNodePurity"])
  importance_df$Gene <- gene
  
  importance_list[[gene]] <- importance_df
}

# 合并所有基因的特征重要性数据
importance_all <- do.call(rbind, importance_list)
write.csv(importance_all, file.path(tables_dir, "Hub_Genes_RF_Importance.csv"), row.names = FALSE)
# 绘制随机森林特征重要性柱状图
png("random_forest_importance.png", width = 1200, height = 800)
ggplot(importance_all, aes(x = reorder(Immune_Cell, Importance), y = Importance, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Feature Importance of Immune Cells for Hub Genes",
       x = "Immune Cell Type",
       y = "Importance Score") +
  scale_fill_brewer(palette = "Set1")
dev.off()


# Load required packages
library(ggplot2)
library(viridis)  # For better color scaling
library(RColorBrewer)

# Use a larger color palette for more distinct colors
num_genes <- length(unique(importance_all$Gene))
color_palette <- colorRampPalette(brewer.pal(9, "Set1"))(num_genes)  # Expanding Set1

# Plot with the updated color scale
ggplot(importance_all, aes(x = reorder(Immune_Cell, Importance), y = Importance, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Feature Importance of Immune Cells for Hub Genes",
       x = "Immune Cell Type",
       y = "Importance Score") +
  scale_fill_manual(values = color_palette)  # Use expanded color palette

# Save the plot
ggsave("random_forest_importance_fixed.png", width = 12, height = 8)


# ---- 基因表达热图 ----
# Z-score 标准化
gene_expr_scaled <- t(scale(t(gene_expr)))

# 定义分组信息
annotation_col <- data.frame(Group = sample_info$Group)
rownames(annotation_col) <- colnames(gene_expr_scaled)

# 颜色
ann_colors <- list(Group = c(High_CP = "red", Low_CP = "blue"))

# 绘制热图
png("gene_expression_heatmap.png", width = 1200, height = 800)
pheatmap(gene_expr_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Hub Gene Expression Heatmap")
dev.off()

# ---- 基因表达箱线图 ----
# 转换数据格式
gene_expr_df <- as.data.frame(t(gene_expr))
gene_expr_df$Group <- sample_info$Group
gene_expr_long <- gather(gene_expr_df, key = "Gene", value = "Expression", -Group)

# 绘制箱线图
png("gene_expression_boxplot.png", width = 1200, height = 800)
ggplot(gene_expr_long, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  stat_compare_means(aes(group = Group), label = "p.format", method = "wilcox.test") +
  theme_bw() +
  scale_fill_manual(values = c("red", "blue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Hub Gene Expression in High_CP and Low_CP Groups",
       x = "Gene",
       y = "Expression Level")
dev.off()





###################
# [模块1] 免疫激活通路活性计算 - 确保结果正确
immune_genes <- c("CD3D","CD3E","CD8A","CD4","HLA-DRA","HLA-DRB1","GZMB","PRF1","IFNG","IL2RA")
immune_genes <- immune_genes[immune_genes %in% rownames(exprSet_new)]

# 检查并转换数据类型
if(!is.numeric(exprSet_new)) exprSet_new <- apply(exprSet_new, 2, as.numeric)

pathway_activity <- colMeans(exprSet_new[immune_genes, ], na.rm = TRUE)
results$Immune_Activation <- as.numeric(pathway_activity)  # 确保为数值型

# [模块2] 相关性分析 - 使用正确列名
# 获取实际的22个免疫细胞列名
cell_types <- colnames(results)[1:22]  # 直接使用数据中的列名

# 计算相关系数和P值
cor_data <- do.call(rbind, lapply(cell_types, function(col) {
  # 确保两列都是数值型
  x <- as.numeric(results[[col]])
  y <- as.numeric(results$Immune_Activation)
  
  ct <- cor.test(x, y, use = "complete.obs")
  data.frame(
    Immune_Cell = col,
    Correlation = ct$estimate,
    P_value = ct$p.value
  )
}))

# 添加显著性标记
cor_data$Significance <- cut(cor_data$P_value, 
                             breaks = c(0, 0.001, 0.01, 0.05, 1),
                             labels = c("***", "**", "*", ""))

# 排序数据
cor_data <- cor_data[order(cor_data$Correlation), ]

# [模块3] 增强可视化

ggplot(cor_data, aes(x = reorder(Immune_Cell, Correlation), 
                     y = Correlation, 
                     fill = Correlation)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.2f%s", Correlation, Significance)),
            hjust = ifelse(cor_data$Correlation > 0, -0.1, 1.1),
            size = 3.5, color = "black") +
  coord_flip() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0) +
  scale_y_continuous(limits = c(min(cor_data$Correlation)*1.1, 
                                max(cor_data$Correlation)*1.1)) +
  labs(title = "Correlation Between Immune Activation Pathway and Immune Cells",
       x = "Immune Cell Type", 
       y = "Pearson Correlation Coefficient") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(color = "black", size = 11),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

# 保存结果
write.csv(cor_data, "Immune_Activation_Correlation.csv")
ggsave("Immune_Activation_Correlation.png", width = 10, height = 8, dpi = 300)











##############################




# [模块3] 免疫浸润与临床特征关联分析

save.image(file = "pre_model_training.RData", 
           compress = "xz", 
           safe = TRUE)
### ---------------------------
### 崩溃恢复方案
### ---------------------------
# 崩溃后重启R，运行：
load("pre_model_training.RData")
library(caret)


##############################already changed to use version two, this is previous draft
### ---------------------------
### 数据预处理与特征工程
### ---------------------------
#library(tidyverse)
library(caret)
library(pROC)
library(xgboost)
library(glmnet)
library(randomForest)
library(kernlab)
library(gbm)
library(e1071)
library(ggpubr)
library(vip)
set.seed(123)  # Ensure reproducibility

### Data Preparation & Preprocessing ---------------------------
X <- results[, 1:22]  # Immune cell proportions
y <- factor(ifelse(sample_info$Group == "High_CP", 1, 0))
y <- factor(y, levels = c(0, 1), labels = c("Low_CP", "High_CP"))

# Identify highly correlated features
cor_matrix <- cor(X)
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.8)
X_filtered <- X[, -highly_correlated]  # Remove correlated features
removed_features <- colnames(X)[highly_correlated]  # Store removed features

# Create stratified folds using createDataPartition
cv_folds <- createDataPartition(y, p = 0.8, list = TRUE, times = 5)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  index = cv_folds,
  preProcOptions = c("center", "scale")  # Adding scaling for KNN
)


### Model Definitions & Training ---------------------------
rf_model <- train(X_filtered, y, method = "rf", metric = "ROC", trControl = ctrl, tuneGrid = expand.grid(mtry = c(3,5,7)))
xgb_model <- train(X_filtered, y, method = "xgbTree", metric = "ROC", trControl = ctrl, 
                   tuneGrid = expand.grid(nrounds = 100, max_depth = c(3,5), eta = c(0.01, 0.1), 
                                          gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8))
en_model <- train(X_filtered, y, method = "glmnet", metric = "ROC", trControl = ctrl, 
                  tuneGrid = expand.grid(alpha = seq(0, 1, 0.2), lambda = 10^seq(-3, 0, length=20)))
lasso_model <- train(X_filtered, y, method = "glmnet", metric = "ROC", trControl = ctrl, 
                     tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length=20)))
logit_model <- train(X_filtered, y, method = "glmnet", metric = "ROC", trControl = ctrl, family = "binomial",
                     tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length = 10)))


# SVM with radial kernel
svm_model <- train(X_filtered, y, method = "svmRadial", metric = "ROC", trControl = ctrl, tuneLength = 10)


# K-Nearest Neighbors (with scaling)
knn_model <- train(X_filtered, y, method = "knn", metric = "ROC", trControl = ctrl, tuneLength = 10, preProcess = "scale")

# Extract ROC values from each model
models <- list(lasso_model, knn_model, svm_model, logit_model, rf_model, xgb_model, en_model)
model_names <- c("lasso","KNN", "SVM", "Logistic", "RF", "XGBoost", "ElasticNet")
results <- data.frame(Model = model_names, ROC = sapply(models, function(model) max(model$results$ROC, na.rm = TRUE)))
print(results)

# Bar plot for ROC comparison
library(ggplot2)
ggplot(results, aes(x = Model, y = ROC, fill = Model)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Model ROC Comparison", x = "Model", y = "ROC")
ggsave(file.path(graphs_dir, "Model_ROC_Comparison.png"), width = 10, height = 6, dpi = 300)

#保存rdata
save.image(file = "model_training_results.RData", compress = "xz", safe = TRUE)

###################################################################################
##################modified version#################################################
# Load necessary libraries (确保这些库已安装)
library(caret)      # 用于机器学习模型训练和评估
library(ggplot2)    # 用于数据可视化
library(dplyr)      # 用于数据处理和管道操作
library(corrplot)   # 用于 findCorrelation 函数

set.seed(123) # 确保结果可重现性

### 数据准备与预处理 ------------------------------------------------------
# 假设 'results' (CIBERSORT输出的免疫细胞比例) 和 'sample_info' (包含Group信息)
# 已在之前的步骤中加载并可用。

# X: 免疫细胞比例特征矩阵 (从 CIBERSORT 结果中提取前22列)
X <- results[, 1:22]

# y: 目标变量，将 'High_CP' 和 'Low_CP' 转换为因子类型，并指定水平和标签
# 确保目标变量是 caret 包期望的二分类因子格式
y <- factor(ifelse(sample_info$Group == "High_CP", 1, 0))
y <- factor(y, levels = c(0, 1), labels = c("Low_CP", "High_CP"))

# 识别并移除高度相关的特征，以减少多重共线性，提高模型稳定性和解释性
# findCorrelation 函数 (来自 'caret' 包) 能够识别并返回需要移除的列名，
# 以确保任何一对特征之间的绝对相关性低于设定的阈值 (此处为 0.8)。
cor_matrix <- cor(X) # 计算特征之间的相关性矩阵
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE) # 获取要移除的特征名称
X_filtered <- X[,!(colnames(X) %in% highly_correlated)] # 从原始特征矩阵中移除这些特征
removed_features <- highly_correlated # 存储被移除的特征名称，以便后续报告

# 创建分层交叉验证折叠，用于模型训练和评估
# 分层抽样确保每个折叠中目标变量 (Low_CP/High_CP) 的比例与原始数据集保持一致，
# 这对于不平衡数据集尤其重要，能提供更稳健的性能评估。
cv_folds <- createDataPartition(y, p = 0.8, list = TRUE, times = 5) # 5折交叉验证，80%训练，20%测试

# 定义 caret::train 函数的训练控制参数
ctrl <- trainControl(
  method = "cv",             # 使用交叉验证 (Cross-Validation)
  number = 5,                # 交叉验证的折叠数设置为 5
  classProbs = TRUE,         # 计算类别概率，这对于 ROC 曲线的计算是必需的
  summaryFunction = twoClassSummary, # 使用 twoClassSummary 函数来计算性能指标，
  # 其中包括 ROC (Area Under the Receiver Operating Characteristic Curve)
  savePredictions = "final", # 保存最终预测结果，以便后续分析
  index = cv_folds,          # 使用预定义的分层交叉验证折叠
  preProcOptions = c("center", "scale") # 对所有模型在训练前进行特征预处理：
  # center (减去均值) 和 scale (除以标准差)，
  # 这有助于许多机器学习算法的收敛和性能。
)


### 模型定义与训练 ------------------------------------------------------
# 使用 caret 框架训练多种机器学习模型，并进行超参数调优，以最大化 ROC 值。

# 1. 随机森林 (Random Forest, RF)
# 'mtry' 是每次分裂时随机采样的变量数量。
rf_model <- train(X_filtered, y, method = "rf", metric = "ROC", trControl = ctrl, tuneGrid = expand.grid(mtry = c(3, 5, 7)))

# 2. eXtreme Gradient Boosting (XGBoost)
# 'nrounds': 提升迭代次数
# 'max_depth': 树的最大深度
# 'eta': 学习率
# 'gamma': 节点分裂所需的最小损失减少
# 'colsample_bytree': 构建每棵树时列的子采样比率
# 'min_child_weight': 子节点所需的最小实例权重和
# 'subsample': 训练实例的子采样比率
xgb_model <- train(X_filtered, y, method = "xgbTree", metric = "ROC", trControl = ctrl,
                   tuneGrid = expand.grid(nrounds = 100, max_depth = c(3, 5), eta = c(0.01, 0.1),
                                          gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8))

# 3. 弹性网络 (Elastic Net, EN) - 一种结合了 Lasso 和 Ridge 惩罚的正则化方法
# 'alpha': 混合参数 (0 = Ridge, 1 = Lasso)。此处尝试从纯 Ridge 到纯 Lasso 的多种组合。
# 'lambda': 正则化强度。
en_model <- train(X_filtered, y, method = "glmnet", metric = "ROC", trControl = ctrl,
                  tuneGrid = expand.grid(alpha = seq(0, 1, 0.2), lambda = 10^seq(-3, 0, length = 20)))

# 4. Lasso 回归 - 弹性网络的一种特殊情况，当 alpha = 1 时 (L1 惩罚)
# 注意：对于二分类结果，glmnet 默认使用 binomial (逻辑回归) 家族。
# 因此，此模型在功能上与下面的 'logit_model' 非常相似，甚至可能完全相同。
lasso_model <- train(X_filtered, y, method = "glmnet", metric = "ROC", trControl = ctrl,
                     tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length = 20)))

# 5. 逻辑回归 (Logistic Regression) - 通过 glmnet 实现，并带有 Lasso 惩罚 (alpha = 1)
# 'family = "binomial"' 明确指定了用于二分类结果的逻辑回归。
# 如上所述，此模型与 'lasso_model' 在此设置下是冗余的。
logit_model <- train(X_filtered, y, method = "glmnet", family = "binomial", metric = "ROC", trControl = ctrl,
                     tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length = 10)))

# 6. 支持向量机 (Support Vector Machine, SVM) - 使用径向基函数 (RBF) 核
# 'tuneLength' 指定要尝试的超参数组合的数量。
svm_model <- train(X_filtered, y, method = "svmRadial", metric = "ROC", trControl = ctrl, tuneLength = 10)

# 7. K-近邻 (K-Nearest Neighbors, KNN)
# 注意：原始代码中的 'preProcess = "scale"' 是冗余的，因为 'center' 和 'scale'
# 已经在 'ctrl$preProcOptions' 中为所有模型全局应用。已将其移除。
knn_model <- train(X_filtered, y, method = "knn", metric = "ROC", trControl = ctrl, tuneLength = 10)


### 结果提取与可视化 ------------------------------------------------------

# 提取每个训练模型的最佳 ROC 值
models_list <- list(lasso_model, knn_model, svm_model, logit_model, rf_model, xgb_model, en_model)
# 为模型提供更具描述性的名称，以便在图表中清晰显示
model_names_display <- c("Lasso","KNN", "SVM", "Logistic Regression", "Random Forest", "XGBoost", "Elastic Net")

# 创建一个数据框来存储模型性能结果
results_df <- data.frame(Model = model_names_display,
                         ROC = sapply(models_list, function(model) max(model$results$ROC, na.rm = TRUE)))

# 按照 ROC 值降序排列结果，以便在条形图中更好地比较模型性能
results_df <- results_df %>%
  arrange(desc(ROC)) %>%
  mutate(Model = factor(Model, levels = Model)) # 将 Model 列转换为因子，以保持 ggplot 中的顺序

print("模型性能 (ROC 值):")
print(results_df)

# 绘制 ROC 比较的条形图 (增强可视化效果)
roc_plot <- ggplot(results_df, aes(x = Model, y = ROC, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) + # 绘制条形图，并调整条形宽度
  geom_text(aes(label = sprintf("%.3f", ROC)), # 在条形图上方添加 ROC 值标签，格式化为三位小数
            vjust = -0.5, # 调整标签的垂直位置，使其略高于条形
            size = 4,     # 调整标签字体大小
            color = "black") + # 设置标签颜色
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # 设置 Y 轴范围从 0 到 1，并定义刻度
  labs(
    title = "Machine Learning models comparison (ROC)", # 图表主标题
    x = "Model",                           # X 轴标签
    y = "ROC", # Y 轴标签
    fill = "Model"                         # 填充颜色图例标题 (虽然在此图中将隐藏)
  ) +
  theme_minimal(base_size = 12) + # 使用简洁主题，并设置基础字体大小
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # 居中、加粗并增大标题字体
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), # 旋转 X 轴标签，右对齐，加粗并调整大小
    axis.text.y = element_text(size = 10), # 调整 Y 轴标签大小
    axis.title = element_text(face = "bold", size = 12), # 加粗轴标题
    legend.position = "none", # 隐藏图例，因为模型名称已在 X 轴上显示，图例会显得冗余
    panel.grid.major.x = element_blank(), # 移除垂直主网格线
    panel.grid.minor.y = element_blank()  # 移除水平次网格线
  ) +
  scale_fill_brewer(palette = "Set2") # 使用一个颜色友好的调色板 (例如 Set2，对色盲友好)

# 保存增强后的图表 (请确保 'graphs_dir' 变量已定义，或直接保存到当前工作目录)
ggsave(
   file.path(graphs_dir, "Model_ROC_Comparison_modified.png"),
   roc_plot,
   width = 10,
   height = 7,
   dpi = 300, # 设置高分辨率 (每英寸点数)
   bg = "white" # 确保背景为纯白色
 )


### Model Performance Evaluation ---------------------------
# Compute predictions and evaluate performance
library(ROCR)
library(caret)
library(pROC)

model_performance <- function(model, X_test, y_test) {
  preds <- predict(model, X_test, type = "prob")[,2]
  pred_classes <- predict(model, X_test)
  
  roc_obj <- roc(y_test, preds, levels = rev(levels(y_test)))
  auc_value <- auc(roc_obj)
  cm <- confusionMatrix(pred_classes, y_test)
  
  return(list(auc = auc_value, cm = cm))
}

# Split data into train and test
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X_filtered[trainIndex, ]
X_test <- X_filtered[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Evaluate each model
performance_list <- lapply(models, model_performance, X_test = X_test, y_test = y_test)
names(performance_list) <- model_names

# Print AUC and Confusion Matrix
for (name in model_names) {
  cat("Model:", name, "\n")
  cat("AUC:", performance_list[[name]]$auc, "\n")
  print(performance_list[[name]]$cm)
  cat("-----------------------------------\n")
}

### ---------------------------
### 1️⃣ 特征重要性分析（修复与优化）
### ---------------------------

get_feature_importance <- function(model, model_name) {
  # 处理不同模型类型
  if (model_name == "rf") {
    # 确保加载caret包
    if (!require(caret, quietly = TRUE)) install.packages("caret")
    
    # 处理caret训练的随机森林
    if (inherits(model, "train")) {
      imp <- caret::varImp(model)$importance
    } else {
      imp <- randomForest::importance(model)
    }
    
    data.frame(
      Feature = rownames(imp),
      Importance = imp[, 1], # 通常是第一列
      Model = "Random Forest"
    )
    
  } else if (model_name == "xgb") {
    # 确保加载xgboost包
    if (!require(xgboost, quietly = TRUE)) install.packages("xgboost")
    
    imp <- xgboost::xgb.importance(model = model$finalModel)
    data.frame(
      Feature = imp$Feature,
      Importance = imp$Gain,
      Model = "XGBoost"
    )
    
  } else if (model_name == "svm") {
    # 关键修改：从caret模型中提取最终模型
    if (inherits(model, "train")) {
      svm_final_model <- model$finalModel
    } else {
      svm_final_model <- model
    }
    
    # 确保加载DALEX包
    if (!require(DALEX, quietly = TRUE)) install.packages("DALEX")
    
    # 为ksvm对象创建自定义预测函数
    ksvm_predict <- function(object, newdata) {
      kernlab::predict(object, newdata = newdata, type = "probabilities")[, 2]
    }
    
    # 创建解释器
    explainer <- DALEX::explain(
      model = svm_final_model,
      data = as.data.frame(X_filtered),
      y = as.numeric(y) - 1,
      predict_function = ksvm_predict, # 使用自定义预测函数
      label = "SVM",
      type = "classification"
    )
    
    # 计算特征重要性
    imp <- DALEX::model_parts(explainer)
    
    imp %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(Importance = mean(abs(dropout_loss))) %>%
      dplyr::transmute(
        Feature = variable,
        Importance,
        Model = "SVM"
      )
    
  } else {
    stop("Unsupported model type")
  }
}
# 运行获取特征重要性的代码 (使用修正后的函数)
importance_list <- list(
  rf = get_feature_importance(rf_model, "rf"),
  xgb = get_feature_importance(xgb_model, "xgb"),
  svm = get_feature_importance(svm_model, "svm")
)

# 获取各模型特征重要性
importance_list <- list(
  rf = get_feature_importance(rf_model, "rf"),
  xgb = get_feature_importance(xgb_model, "xgb"),
  svm = get_feature_importance(svm_model, "svm")
)

# 合并重要性数据
importance_combined <- bind_rows(importance_list) %>%
  group_by(Feature) %>%
  mutate(Global_Importance = mean(Importance)) %>%
  ungroup()

### ---------------------------
### 2️⃣ 特征重要性可视化（优化版）
### ---------------------------

# 自定义科学配色
model_colors <- c(
  "Random Forest" = "#1F77B4",
  "XGBoost" = "#FF7F0E",
  "SVM" = "#2CA02C"
)

# 绘制特征重要性对比图
ggplot(importance_combined, aes(x = Importance, y = reorder(Feature, Global_Importance), fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Comparative Feature Importance Analysis",
    x = "Normalized Importance Score",
    y = "Immune Cell Type",
    caption = "Importance scores normalized per model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  ) +
  guides(fill = guide_legend(nrow = 1))

# 保存图形
ggsave("feature_importance_comparison.png", width = 10, height = 8, dpi = 300)










##############################
### 新增：WGCNA共表达网络 ###
##############################
library(WGCNA)


# 启用多线程计算
enableWGCNAThreads()

# 准备表达数据
datExpr <- t(exprSet_new)  # 转置：行是样本，列是基因
datExpr <- apply(datExpr, 2, as.numeric)  # 确保数值型
rownames(datExpr) <- colnames(exprSet_new)  # 样本名为行名

# 1. 确定软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "unsigned",
  corFnc = cor,  # 显式使用基础cor函数
  corOptions = list(use = "pairwise.complete.obs")
)

# 2. 绘制软阈值选择图
png(file.path(graphs_dir, "Soft_Threshold_Selection.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))

# 图1: Scale independence
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red")

# 图2: Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, cex = 0.9, col = "red")

dev.off()

# 3. 选择软阈值
if (is.na(sft$powerEstimate)) {
  power <- 6  # 默认值
  cat("No power reached R² > 0.9, using default power =", power, "\n")
} else {
  power <- sft$powerEstimate
  cat("Selected power =", power, "\n")
}

# 4. 构建共表达网络 (关键修复)
net <- blockwiseModules(
  datExpr,
  power = power,
  networkType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(working_dir, "EBV_TOM"),
  verbose = 3,
  
  # 关键修复参数
  corType = "pearson",
  corFnc = cor,  # 显式使用基础cor函数
  corOptions = list(use = "pairwise.complete.obs"),
  
  # 防止内部传递额外参数
  weights = NULL,
  TOMType = "unsigned"
)

# 5. 可视化模块
moduleColors <- labels2colors(net$colors)

png(file.path(graphs_dir, "Module_Dendrogram.png"), width = 1200, height = 800)
plotDendroAndColors(
  net$dendrograms[[1]], 
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# 保存模块分配结果
module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = net$colors,
  ModuleColor = moduleColors
)
write.csv(module_df, file.path(tables_dir, "WGCNA_Module_Assignments.csv"))



# 4. 模块-性状关联分析
# 准备性状数据（EBV拷贝数组）
traitData <- data.frame(
  EBV_CP = ifelse(sample_info$Group == "High_CP", 1, 0),
  row.names = rownames(sample_info)
)

# 确保性状数据与表达数据样本顺序一致
traitData <- traitData[rownames(datExpr), , drop = FALSE]

# 计算模块特征基因（module eigengenes）
MEs <- net$MEs

# 计算模块与性状的相关性
moduleTraitCor <- cor(MEs, traitData, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 可视化关联热图
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)

# 5. 关键模块分析（以与EBV_CP最相关的模块为例）
# 找出与EBV拷贝数最相关的模块
ebv_cor <- moduleTraitCor[, "EBV_CP"]
ebv_p <- moduleTraitPvalue[, "EBV_CP"]

# 选择最显著的模块
signif_modules <- which(ebv_p < 0.05)
if (length(signif_modules) > 0) {
  top_module <- names(which.max(abs(ebv_cor[signif_modules])))
  
  # 提取模块基因
  module_genes <- colnames(datExpr)[net$colors == as.integer(gsub("ME", "", top_module))]
  
  # 保存关键模块基因
  write.csv(module_genes, file.path(tables_dir, paste0("KeyModule_", top_module, "_Genes.csv")))
  
  # 关键模块的GO/KEGG分析
  entrez_ids <- mapIds(org.Hs.eg.db, keys = module_genes, 
                       column = "ENTREZID", keytype = "SYMBOL")
  entrez_ids <- na.omit(entrez_ids)
  
  # GO富集分析
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # 可视化GO结果
  dotplot(ego, showCategory=15, title=paste("GO Enrichment -", top_module))
  ggsave(file.path(graphs_dir, paste0("GO_Enrichment_", top_module, ".png")))
  
  # KEGG富集分析
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  # 可视化KEGG结果
  dotplot(ekegg, showCategory=15, title=paste("KEGG Enrichment -", top_module))
  ggsave(file.path(graphs_dir, paste0("KEGG_Enrichment_", top_module, ".png")))
}



##############################
### 增强版机器学习可解释性 ###
##############################
library(DALEX)
library(modelDown)

# 对最佳模型(RF)进行解释
explainer_rf <- DALEX::explain(
  model = rf_model$finalModel,
  data = X_filtered,
  y = y,
  label = "Random Forest"
)


#######################################
### ---------------------------
### 7. 免疫相关基因表达与临床表型关联分析
### ---------------------------

# 定义免疫激活基因（与前面一致）
immune_genes <- c("CD3D","CD3E","CD8A","CD4","HLA-DRA","HLA-DRB1","GZMB","PRF1","IFNG","IL2RA")
immune_genes <- immune_genes[immune_genes %in% rownames(exprSet_new)]

# 提取表达数据并标准化
immune_expression_data <- exprSet_new[immune_genes, ]
immune_expr_scaled <- t(scale(t(immune_expression_data)))

# 获取临床响应信息
response_groups <- pheno_data$characteristics_ch1.7

# 创建响应注释
response_annotation <- data.frame(
  geo_accession = rownames(pheno_data),
  response_group = ifelse(grepl("PR", response_groups), "PR",
                          ifelse(grepl("SD", response_groups), "SD",
                                 ifelse(grepl("PD", response_groups), "PD", NA)))
)

# 创建热图注释
annotation_col <- data.frame(Response = response_annotation$response_group)
rownames(annotation_col) <- response_annotation$geo_accession

# 定义颜色
response_colors <- list(Response = c(PR = "green", SD = "orange", PD = "red"))

# 绘制免疫激活基因热图
pheatmap(immune_expr_scaled,
         annotation_col = annotation_col,
         annotation_colors = response_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         main = "Immune Activation Genes by Response Group",
         fontsize_row = 10)

### ---------------------------
### 8. EBV感染通路分析
### ---------------------------

# 定义EBV感染通路基因
entrez_ids <- c("596", "919", "915", "916", "917", "3109", "3113", "3115")
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID")
print(gene_symbols)

# 提取表达数据
kegg_expression_data <- exprSet_new[rownames(exprSet_new) %in% gene_symbols, ]

# 计算通路活性并添加到CIBERSORT结果
pathway_activity <- colMeans(kegg_expression_data)
results$Pathway_Activity <- pathway_activity

# 计算相关性
correlation <- cor(results[, 1:22], results$Pathway_Activity)
correlation_df <- data.frame(
  Immune_Cell = rownames(correlation),
  Correlation = correlation[, 1]
)
correlation_df <- correlation_df[order(correlation_df$Correlation), ]

# 绘制相关性条形图
ggplot(correlation_df, aes(x = reorder(Immune_Cell, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Correlation Between Immune Cells and EBV Pathway Activity",
       x = "Immune Cell Type",
       y = "Pearson Correlation") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.y = element_text(size = 10))

# 标准化表达数据
kegg_expr_scaled <- t(scale(t(kegg_expression_data)))

# 创建EBV分组注释
ebv_annotation <- data.frame(
  EBV_Group = sample_info$Group,
  row.names = rownames(sample_info)
)

# 合并响应和EBV注释
combined_annotation <- merge(response_annotation, ebv_annotation, by.x = "geo_accession", by.y = "row.names")
rownames(combined_annotation) <- combined_annotation$geo_accession
combined_annotation <- combined_annotation[, -1]  # 移除geo_accession列

# 定义颜色方案
response_colors <- list(Response = c(PR = "green", SD = "orange", PD = "red"))
ebv_colors <- list(EBV_Group = c(High_CP = "blue", Low_CP = "yellow"))
annotation_colors <- c(response_colors, ebv_colors)

# 绘制EBV通路基因热图
pheatmap(kegg_expr_scaled,
         annotation_col = combined_annotation,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         main = "EBV Pathway Genes by Response and EBV Copy Number",
         fontsize_row = 10)

### ---------------------------
### 9. 关键免疫细胞识别 (随机森林)
### ---------------------------

# 准备数据
X <- results[, 1:22]  # 免疫细胞比例
y <- results$Pathway_Activity  # EBV通路活性

# 训练随机森林模型
set.seed(123)
rf_model <- randomForest(x = X, y = y, importance = TRUE)

# 提取特征重要性
importance_scores <- importance(rf_model)

# 获取重要性最高的5个免疫细胞
top_immune_cells <- rownames(importance_scores)[
  order(importance_scores[, "%IncMSE"], decreasing = TRUE)][1:5]

cat("Top 5 immune cells predicting EBV pathway activity:\n")
print(top_immune_cells)

# 保存结果
write.csv(data.frame(
  Immune_Cell = rownames(importance_scores),
  Importance = importance_scores[, "%IncMSE"]
), file.path(tables_dir, "RF_Importance_EBV_Pathway.csv"))

### ---------------------------
### 10. 中性粒细胞趋化基因分析
### ---------------------------

# 定义中性粒细胞趋化基因
neutrophil_genes <- c("S100A8", "S100A9", "CXCR2", "CSF3R", "CXCL1", "CXCL2")
neutrophil_genes <- neutrophil_genes[neutrophil_genes %in% rownames(exprSet_new)]

# 提取表达数据并标准化
neutrophil_expression_data <- exprSet_new[neutrophil_genes, ]
neutrophil_expr_scaled <- t(scale(t(neutrophil_expression_data)))

# 创建分组注释
annotation_col <- data.frame(Group = sample_info$Group)
rownames(annotation_col) <- colnames(neutrophil_expr_scaled)

# 定义颜色
ann_colors <- list(Group = c(High_CP = "#E64B35FF", Low_CP = "#4DBBD5FF"))

# 绘制热图
pheatmap(neutrophil_expr_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Neutrophil Chemotaxis Genes by EBV Copy Number Group",
         fontsize_row = 10)

# 计算相关性
correlation <- cor(t(neutrophil_expression_data), results[, 1:22])
correlation_df <- as.data.frame(correlation)
correlation_df$Gene <- rownames(correlation_df)

# 可视化相关性
correlation_melt <- reshape2::melt(correlation_df, id.vars = "Gene")
ggplot(correlation_melt, aes(x = variable, y = Gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation Between Neutrophil Genes and Immune Cells",
       x = "Immune Cell Type",
       y = "Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))



############################################################################################################################

# 1. 定义你的 Hub Genes 列表 (确保名字与 exprSet_new 的 rownames 完全一致)
hub_genes_list <- c("C5AR1", "CD247", "CD3E", "CSF3R", "CXCL8", "CXCR2", 
                    "HLA-DPB1", "IL1A", "IL1B", "IL2RG", "ITGAX", "LCK", 
                    "MS4A1", "NLRP3", "S100A8", "S100A9", "THBD", "TREM1")

# 2. 从 exprSet_new 中过滤出这些基因
# 检查哪些基因在矩阵中存在
existing_genes <- hub_genes_list[hub_genes_list %in% rownames(exprSet_new)]
hub_expr_matrix <- exprSet_new[existing_genes, ]

# 3. 数据转置：将行(基因)转为列(特征)，行变为样本名
hub_expr_df <- as.data.frame(t(hub_expr_matrix))

# 4. 加入 Sample_ID 列 (确保与你 CIBERSORT 结果的 ID 对齐)
hub_expr_df$Sample_ID <- rownames(hub_expr_df)

# 5. 调整列顺序，让 Sample_ID 在第一列
hub_expr_df <- hub_expr_df[, c("Sample_ID", existing_genes)]

# 6. 保存 CSV
write.csv(hub_expr_df, file.path(tables_dir, "Hub_Genes_Expression.csv"), row.names = FALSE)

cat("Hub_Genes_Expression.csv 已成功保存至:", tables_dir, "\n")


###############################################
### ---------------------------
### 11. Hub基因与免疫细胞关联分析 (新增)
### ---------------------------

# 定义hub基因列表 (根据您提供的热图)
hub_genes <- c("C5AR1", "CD247", "CD3E", "CSF3R", "CXCL8", "CXCR2", 
               "HLA-DPB1", "IL1A", "IL1B", "IL2RG", "ITGAX", "LCK", 
               "MS4A1", "NLRP3", "S100A8", "S100A9", "THBD", "TREM1")

# 确保基因在表达矩阵中
hub_genes <- hub_genes[hub_genes %in% rownames(exprSet_new)]

# 提取hub基因表达数据
hub_expr <- exprSet_new[hub_genes, ]

# 准备免疫细胞比例数据 (来自CIBERSORT)
immune_cells <- results[, 1:22]

# 创建存储重要性结果的列表
importance_list <- list()

# 对每个hub基因进行随机森林分析
for (gene in hub_genes) {
  # 提取当前基因的表达量
  gene_expr <- as.numeric(hub_expr[gene, ])
  
  # 训练随机森林模型
  rf_model <- randomForest(
    x = immune_cells,
    y = gene_expr,
    importance = TRUE,
    ntree = 1000
  )
  
  # 提取重要性分数
  imp_scores <- importance(rf_model)
  
  # 存储结果
  importance_list[[gene]] <- imp_scores[, "%IncMSE"]
}

# 转换结果为数据框
importance_df <- as.data.frame(importance_list)
rownames(importance_df) <- colnames(immune_cells)

# 保存结果
write.csv(importance_df, file.path(tables_dir, "Hub_Genes_Immune_Importance.csv"))

### ---------------------------
### 12. 生成热图可视化 (修复版)
### ---------------------------

# 数据准备
plot_data <- as.matrix(importance_df)

# 创建行注释数据框
row_annot <- data.frame(
  Immune_Cell_Type = rownames(plot_data)
  rownames(row_annot) <- rownames(plot_data)
  
  # 创建列注释数据框
  col_annot <- data.frame(
    Gene_Group = ifelse(
      grepl("CD|HLA|IL|ITG|LCK|MS4A|THBD", colnames(plot_data)), 
      "Immune Signaling",
      "Inflammatory Response")
  )
  rownames(col_annot) <- colnames(plot_data)
  
  # 定义颜色方案
  cell_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  group_colors <- list(
    Gene_Group = c("Immune Signaling" = "#1f77b4", "Inflammatory Response" = "#ff7f0e")
  )
  
  # 创建热图
  pheatmap(plot_data,
           color = cell_colors,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           annotation_row = row_annot,
           annotation_col = col_annot,
           annotation_colors = group_colors,
           show_colnames = TRUE,
           show_rownames = TRUE,
           fontsize_row = 10,
           fontsize_col = 10,
           main = "Feature Importance of Immune Cells for Hub Genes",
           filename = file.path(graphs_dir, "Hub_Genes_Immune_Importance_Heatmap.png"),
           width = 12,
           height = 8)
  
  ### ---------------------------
  ### 13. 交互式热图可视化 (新增)
  ### ---------------------------
  
  # 安装并加载plotly
  if (!require("plotly")) install.packages("plotly")
  library(plotly)
  
  # 创建交互式热图
  heatmap_plot <- plot_ly(
    x = colnames(plot_data),
    y = rownames(plot_data),
    z = plot_data,
    type = "heatmap",
    colors = cell_colors,
    hoverinfo = "x+y+z",
    showscale = TRUE
  ) %>%
    layout(
      title = "Feature Importance of Immune Cells for Hub Genes",
      xaxis = list(title = "Hub Genes"),
      yaxis = list(title = "Immune Cell Types"),
      margin = list(l = 200)  # 增加左边距以显示完整的免疫细胞名称
    )
  
  # 保存为HTML文件
  htmlwidgets::saveWidget(
    heatmap_plot,
    file.path(graphs_dir, "Hub_Genes_Immune_Importance_Interactive.html")
  )
  
  # 打印到R Viewer
  print(heatmap_plot)

#####################################################################################

  # ---------------------------
  # Survival Analysis with Cox and KM Curves (Publication Ready)
  # ---------------------------
  
  # Load required libraries
  library(GEOquery)
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  
  # Step 1: Load and clean data
  # ------------------------------------------
  # Load clinical phenotype data
  geo <- getGEO("GSE102349", GSEMatrix = TRUE)
  pheno_data <- pData(geo[[1]])
  
  # Load manually downloaded gene expression matrix
  expr_data <- read.delim("GSE102349_NPC_mRNA_processed.txt", sep = "\t", row.names = 1)
  colnames(expr_data) <- rownames(pheno_data)
  
  
  
  # Define the comprehensive list of genes from up_intersect and down_intersect
  up_intersect <- c("CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
  down_intersect <- c("CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX")
  
  # Combine all genes for survival analysis
  genes_for_survival_analysis <- unique(c(up_intersect, down_intersect))
  # Define target genes
  sig_genes <- c("CD247","CXCR2","CXCL1","CSF3R","S100A8","CD4", "CD8A", "CD3D", "CD3E", "IL2RA", "HLA-DRA","HLA-DPB1")
  sig_genes <- genes_for_survival_analysis
  
  matched_genes <- sig_genes[sig_genes %in% rownames(expr_data)]
  expr_mat <- t(expr_data[matched_genes, ])
  expr_df <- as.data.frame(expr_mat)
  expr_df$geo_accession <- rownames(expr_df)
  
  # Clean survival data
  pheno_data$`time to event:ch1`[pheno_data$`time to event:ch1` == "N/A"] <- NA
  pheno_data$`time to event:ch1` <- as.numeric(pheno_data$`time to event:ch1`)
  
  surv_df <- pheno_data %>%
    select(geo_accession, `time to event:ch1`, `event:ch1`) %>%
    rename(time = `time to event:ch1`, event = `event:ch1`) %>%
    mutate(event = ifelse(event == "Disease progression", 1, 0)) %>%
    filter(!is.na(time))
  
  # Merge expression and survival
  merged <- merge(surv_df, expr_df, by = "geo_accession")
  
  # Prepare results container
  surv_results <- data.frame(gene = character(), coef = numeric(), p_value = numeric())
  graphs_dir <- "KM_plots"
  tables_dir <- "results"
  if (!dir.exists(graphs_dir)) dir.create(graphs_dir)
  if (!dir.exists(tables_dir)) dir.create(tables_dir)
  
  # Loop through genes
  for (gene in matched_genes) {
    if (!gene %in% colnames(merged)) next
    
    merged$group <- ifelse(merged[[gene]] > median(merged[[gene]], na.rm = TRUE), "High", "Low")
    surv_obj <- Surv(merged$time, merged$event)
    fit <- survfit(surv_obj ~ group, data = merged)
    model <- coxph(surv_obj ~ merged[[gene]], data = merged)
    model_sum <- summary(model)
    
    # Save results
    coef <- model_sum$coefficients[1, "coef"]
    pval <- model_sum$coefficients[1, "Pr(>|z|)"]
    surv_results <- rbind(surv_results, data.frame(gene = gene, coef = coef, p_value = pval))
    
    # Plot KM
    p <- ggsurvplot(
      fit,
      data = merged,
      conf.int = TRUE,
      pval = TRUE,
      pval.method = TRUE,
      risk.table = TRUE,
      palette = c("#E41A1C", "#377EB8"),
      title = paste("Survival Curve -", gene),
      legend.title = "Expression",
      legend.labs = c("High", "Low"),
      xlab = "Time to Event (Days)",
      ylab = "Survival Probability",
      font.main = c(14, "bold"),
      font.x = c(12),
      font.y = c(12),
      font.legend = c(12),
      font.tickslab = 12,
      ggtheme = theme_minimal(base_size = 14)
    )
    
    ggsave(filename = file.path(graphs_dir, paste0("KM_", gene, ".png")), plot = print(p), width = 8, height = 6, dpi = 300)
  }
  
  # Save Cox result table
  write.csv(surv_results, file.path(tables_dir, "cox_results.csv"), row.names = FALSE)
  
  # Highlight significant genes
  sig_genes_df <- surv_results %>% filter(p_value < 0.05)
  print(sig_genes_df)
  
