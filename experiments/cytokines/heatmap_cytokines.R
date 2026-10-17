# Create a publication-style cytokine heatmap.

required <- c("pheatmap", "RColorBrewer", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
input_file <- Sys.getenv(
  "CYTOKINE_INPUT",
  file.path(project_root, "EBV_recoding", "cytokine_replicate.csv")
)
output_dir <- Sys.getenv("CYTOKINE_OUTPUT_DIR", file.path(script_dir, "output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_raw <- suppressWarnings(utils::read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
))
data_raw[] <- lapply(data_raw, as.numeric)
column_names <- colnames(data_raw)
cytokines <- unique(sub("_R[0-9]+$", "", column_names))

averaged <- lapply(cytokines, function(cytokine) {
  columns <- column_names[sub("_R[0-9]+$", "", column_names) == cytokine]
  rowMeans(data_raw[, columns, drop = FALSE], na.rm = TRUE)
})
names(averaged) <- cytokines
matrix_raw <- t(as.matrix(as.data.frame(averaged, row.names = rownames(data_raw))))
matrix_scaled <- t(scale(t(matrix_raw)))
matrix_scaled[!is.finite(matrix_scaled)] <- 0
display_names <- rownames(matrix_scaled)
display_names <- gsub("α", "-alpha", display_names, fixed = TRUE)
display_names <- gsub("β", "-beta", display_names, fixed = TRUE)
display_names <- gsub("γ", "-gamma", display_names, fixed = TRUE)
rownames(matrix_scaled) <- display_names
palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

draw_heatmap <- function(filename = NA) {
  pheatmap::pheatmap(
    matrix_scaled,
    color = palette,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    border_color = "white",
    fontsize = 8,
    fontsize_row = 7,
    fontsize_col = 8,
    angle_col = 45,
    treeheight_row = 20,
    treeheight_col = 15,
    filename = filename,
    width = 5.5,
    height = 6.5,
    silent = !is.na(filename)
  )
}

draw_heatmap(file.path(output_dir, "heatmap_nature_style.pdf"))
ragg::agg_tiff(
  file.path(output_dir, "heatmap_nature_style.tiff"),
  width = 5.5,
  height = 6.5,
  units = "in",
  res = 600,
  compression = "lzw"
)
draw_heatmap()
grDevices::dev.off()
grDevices::png(
  file.path(output_dir, "heatmap_nature_style.png"),
  width = 1650,
  height = 1950,
  res = 300
)
draw_heatmap()
grDevices::dev.off()
utils::write.csv(matrix_raw, file.path(output_dir, "cytokine_mean_matrix.csv"))
cat("Output:", output_dir, "\n")
