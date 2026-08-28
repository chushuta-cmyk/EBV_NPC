# -------------------------------------------------------------------------
# Publication-Ready SCI Heatmap (Nature / Lancet style)
# Works with your data: rows = conditions, columns = cytokine replicates
# -------------------------------------------------------------------------

# 1. Load required libraries
library(pheatmap)
library(tidyverse)   # for data manipulation
library(RColorBrewer)

# 2. Use your existing data_raw (already loaded)
#    If not loaded, read CSV with:
data_raw <- read.csv("cytokine_replicate.csv", row.names = 1, check.names = FALSE)

# 3. Average replicates for each cytokine
#    Columns are like "IL5_R1", "IL5_R2", "IL5_R3", "IL17_R1", ...
#    We need to group by the cytokine name (before "_R")
#    Then average across replicate columns for each condition (row)

# First, get unique cytokine base names
col_names <- colnames(data_raw)
cytokine_bases <- unique(sub("_R[0-9]+$", "", col_names))

# Create an empty list to store averaged data
averaged_list <- list()

for (cyt in cytokine_bases) {
  # Select columns that start with this cytokine name
  rep_cols <- grep(paste0("^", cyt, "_R"), col_names, value = TRUE)
  if (length(rep_cols) > 0) {
    # Compute row-wise mean of replicates
    averaged_list[[cyt]] <- rowMeans(data_raw[, rep_cols, drop = FALSE], na.rm = TRUE)
  }
}

# Convert list to data frame (rows = conditions, cols = cytokines)
df_avg <- as.data.frame(averaged_list)
rownames(df_avg) <- rownames(data_raw)   # conditions: UT, DMSO, TPA, MOCK, Zta

# 4. Transpose so that cytokines are rows and conditions are columns (standard for heatmaps)
mat <- t(as.matrix(df_avg))

# 5. Row scaling (z-score) – shows relative up/down regulation per cytokine
mat_scaled <- t(scale(t(mat)))   # scale rows

# 6. Choose colour palette (Nature style: blue-white-red diverging)
palette_nature <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# 7. Generate heatmap in RStudio plot pane
p <- pheatmap(
  mat_scaled,
  color = palette_nature,
  cluster_rows = TRUE,           # cluster cytokines
  cluster_cols = TRUE,           # cluster conditions
  clustering_method = "complete",
  border_color = "white",        # thin white lines
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,                # rotate condition labels
  treeheight_row = 20,
  treeheight_col = 15,
  main = NA,
  silent = FALSE                 # show plot
)

# 8. Save as vector PDF (for Nature, Lancet, etc.)
pdf("heatmap_nature_style.pdf", width = 5.5, height = 6.5)
pheatmap(
  mat_scaled,
  color = palette_nature,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "white",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  treeheight_row = 20,
  treeheight_col = 15,
  main = NA
)
dev.off()

# 9. Also save as high-resolution TIFF (if raster required)
tiff("heatmap_nature_style.tiff", width = 5.5, height = 6.5, units = "in", res = 600, compression = "lzw")
pheatmap(
  mat_scaled,
  color = palette_nature,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "white",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  treeheight_row = 20,
  treeheight_col = 15,
  main = NA
)
dev.off()

cat("Heatmap displayed. PDF and TIFF saved to working directory.\n")