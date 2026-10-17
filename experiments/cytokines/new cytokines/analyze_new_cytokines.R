# Analyze the co-culture cytokine workbook.

required <- c("readxl", "dplyr", "tidyr", "ggplot2", "pheatmap", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
arguments <- commandArgs(trailingOnly = TRUE)

argument_value <- function(name, default) {
  match <- grep(paste0("^--", name, "="), arguments, value = TRUE)
  if (length(match) == 0) return(default)
  sub(paste0("^--", name, "="), "", match[[1]])
}

input_file <- argument_value("input", file.path(script_dir, "共培养上清液.xlsx"))
output_dir <- argument_value("output-dir", file.path(script_dir, "output"))
sample_map_file <- argument_value("sample-map", file.path(script_dir, "sample_map.csv"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- readxl::read_excel(input_file, col_names = FALSE, .name_repair = "minimal")
if (nrow(raw) < 3 || ncol(raw) < 2) {
  stop("The workbook does not contain the expected cutoff, header, and sample rows.")
}

cytokines <- make.unique(as.character(unlist(raw[2, -1])), sep = "_")
cutoffs <- suppressWarnings(as.numeric(unlist(raw[1, -1])))
samples <- as.character(unlist(raw[-c(1, 2), 1]))
values <- as.data.frame(raw[-c(1, 2), -1], check.names = FALSE)
colnames(values) <- cytokines
rownames(values) <- samples
values[] <- lapply(values, function(column) suppressWarnings(as.numeric(column)))

if (any(!is.finite(cutoffs)) || any(cutoffs <= 0)) {
  stop("Each cytokine needs a positive numeric cutoff.")
}
if (anyNA(values)) {
  stop("The sample table contains missing or non-numeric measurements.")
}

sample_map <- data.frame(Sample = samples, Condition = samples)
if (file.exists(sample_map_file)) {
  supplied_map <- utils::read.csv(sample_map_file, stringsAsFactors = FALSE)
  if (!all(c("Sample", "Condition") %in% colnames(supplied_map))) {
    stop("The sample map needs Sample and Condition columns.")
  }
  if (!all(samples %in% supplied_map$Sample)) {
    stop("The sample map does not contain every workbook sample.")
  }
  sample_map <- supplied_map[match(samples, supplied_map$Sample), c("Sample", "Condition")]
}

long <- data.frame(Sample = samples, values, check.names = FALSE) |>
  tidyr::pivot_longer(-Sample, names_to = "Cytokine", values_to = "Concentration") |>
  dplyr::left_join(sample_map, by = "Sample")
cutoff_table <- data.frame(Cytokine = cytokines, Cutoff = cutoffs)
long <- long |>
  dplyr::left_join(cutoff_table, by = "Cytokine") |>
  dplyr::mutate(
    Fold_over_cutoff = Concentration / Cutoff,
    Log2_fold_over_cutoff = log2(Fold_over_cutoff),
    Above_cutoff = Concentration >= Cutoff
  )

summary_data <- long |>
  dplyr::group_by(Condition, Cytokine, Cutoff) |>
  dplyr::summarise(
    Mean = mean(Concentration),
    SD = if (dplyr::n() > 1) stats::sd(Concentration) else NA_real_,
    Mean_fold_over_cutoff = mean(Fold_over_cutoff),
    Samples = dplyr::n(),
    .groups = "drop"
  )

ranking <- long |>
  dplyr::group_by(Cytokine, Cutoff) |>
  dplyr::summarise(
    Mean_concentration = mean(Concentration),
    Mean_fold_over_cutoff = mean(Fold_over_cutoff),
    Above_cutoff_samples = sum(Above_cutoff),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(Mean_fold_over_cutoff))
top_cytokines <- head(ranking$Cytokine, 8)

plot_summary <- summary_data[summary_data$Cytokine %in% top_cytokines, ]
plot_points <- long[long$Cytokine %in% top_cytokines, ]
plot_cutoffs <- cutoff_table[cutoff_table$Cytokine %in% top_cytokines, ]
plot_summary$Cytokine <- factor(plot_summary$Cytokine, levels = top_cytokines)
plot_points$Cytokine <- factor(plot_points$Cytokine, levels = top_cytokines)
plot_cutoffs$Cytokine <- factor(plot_cutoffs$Cytokine, levels = top_cytokines)

condition_colors <- c(DMSO = "#B3B3B3", TPA = "#F28E00", UT = "#6AAED6")
unknown <- setdiff(unique(plot_summary$Condition), names(condition_colors))
if (length(unknown) > 0) {
  condition_colors <- c(
    condition_colors,
    stats::setNames(grDevices::hcl.colors(length(unknown), "Dark 3"), unknown)
  )
}

bar_plot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(Condition, Mean, fill = Condition)
) +
  ggplot2::geom_col(width = 0.66, color = "black", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
    width = 0.2,
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = plot_points,
    ggplot2::aes(Condition, Concentration),
    inherit.aes = FALSE,
    position = ggplot2::position_jitter(width = 0.08),
    size = 1.8,
    shape = 21,
    fill = "white"
  ) +
  ggplot2::geom_hline(
    data = plot_cutoffs,
    ggplot2::aes(yintercept = Cutoff),
    color = "#B2182B",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  ggplot2::facet_wrap(
    ~Cytokine,
    scales = "free_y",
    ncol = 4,
    labeller = ggplot2::labeller(Cytokine = function(x) {
      x <- gsub("α", "-alpha", x, fixed = TRUE)
      x <- gsub("β", "-beta", x, fixed = TRUE)
      gsub("γ", "-gamma", x, fixed = TRUE)
    })
  ) +
  ggplot2::scale_fill_manual(values = condition_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(
    title = "Co-culture cytokine concentrations",
    subtitle = "Dashed lines show assay cutoffs",
    x = NULL,
    y = "Concentration (pg/mL)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "#F7F7F7", color = NA),
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#E6E6E6", color = "black"),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = ggplot2::element_text(face = "italic", hjust = 0.5),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(output_dir, "Nature_style_new_cytokine_barplot.pdf"),
  bar_plot,
  width = 8.5,
  height = 5.8,
  useDingbats = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "Nature_style_new_cytokine_barplot.tiff"),
  bar_plot,
  width = 8.5,
  height = 5.8,
  dpi = 600,
  compression = "lzw",
  device = ragg::agg_tiff
)
ggplot2::ggsave(
  file.path(output_dir, "Nature_style_new_cytokine_barplot.png"),
  bar_plot,
  width = 8.5,
  height = 5.8,
  dpi = 300
)

heatmap_matrix <- as.matrix(values)
heatmap_matrix <- t(heatmap_matrix)
heatmap_matrix <- log2(sweep(heatmap_matrix, 1, cutoffs, "/"))
display_names <- rownames(heatmap_matrix)
display_names <- gsub("α", "-alpha", display_names, fixed = TRUE)
display_names <- gsub("β", "-beta", display_names, fixed = TRUE)
display_names <- gsub("γ", "-gamma", display_names, fixed = TRUE)
rownames(heatmap_matrix) <- display_names
palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
heatmap_args <- list(
  mat = heatmap_matrix,
  color = palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "white",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  main = "log2 concentration relative to cutoff"
)
do.call(
  pheatmap::pheatmap,
  c(heatmap_args, list(
    filename = file.path(output_dir, "heatmap_new_cytokines.pdf"),
    width = 6.5,
    height = 7.5
  ))
)
ragg::agg_tiff(
  file.path(output_dir, "heatmap_new_cytokines.tiff"),
  width = 6.5,
  height = 7.5,
  units = "in",
  res = 600,
  compression = "lzw"
)
do.call(pheatmap::pheatmap, heatmap_args)
grDevices::dev.off()
grDevices::png(
  file.path(output_dir, "heatmap_new_cytokines.png"),
  width = 1950,
  height = 2250,
  res = 300
)
do.call(pheatmap::pheatmap, heatmap_args)
grDevices::dev.off()

utils::write.csv(long, file.path(output_dir, "new_cytokines_processed.csv"), row.names = FALSE)
utils::write.csv(summary_data, file.path(output_dir, "new_cytokines_summary.csv"), row.names = FALSE)
utils::write.csv(ranking, file.path(output_dir, "new_cytokines_ranking.csv"), row.names = FALSE)

if (!file.exists(sample_map_file)) {
  message("No sample map was found. Samples were plotted as separate conditions.")
}
cat("Output:", output_dir, "\n")
