# Create a publication-style cytokine bar plot.

required <- c("ggplot2", "dplyr", "tidyr", "ragg")
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
rows <- lapply(cytokines, function(cytokine) {
  columns <- column_names[sub("_R[0-9]+$", "", column_names) == cytokine]
  values <- data_raw[, columns, drop = FALSE]
  data.frame(
    Condition = rep(rownames(values), each = ncol(values)),
    Cytokine = cytokine,
    Replicate = rep(seq_len(ncol(values)), times = nrow(values)),
    Value = as.numeric(t(values)),
    stringsAsFactors = FALSE
  )
})
replicates <- do.call(rbind, rows)

summary_data <- replicates |>
  dplyr::group_by(Condition, Cytokine) |>
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = stats::sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

fold_change <- summary_data |>
  dplyr::select(Condition, Cytokine, Mean) |>
  tidyr::pivot_wider(names_from = Condition, values_from = Mean)
if (!all(c("TPA", "UT") %in% names(fold_change))) {
  stop("The cytokine data need TPA and UT rows.")
}
fold_change$TPA_to_UT <- fold_change$TPA / fold_change$UT
fold_change <- fold_change[order(fold_change$TPA_to_UT, decreasing = TRUE), ]
top_cytokines <- head(fold_change$Cytokine, 8)

plot_data <- summary_data[summary_data$Cytokine %in% top_cytokines, ]
plot_data$Cytokine <- factor(plot_data$Cytokine, levels = top_cytokines)
point_data <- replicates[replicates$Cytokine %in% top_cytokines, ]
point_data$Cytokine <- factor(point_data$Cytokine, levels = top_cytokines)

condition_colors <- c(
  UT = "#6AAED6",
  DMSO = "#B3B3B3",
  TPA = "#F28E00",
  MOCK = "#2C7FB8",
  Zta = "#18A87D"
)
unknown <- setdiff(unique(plot_data$Condition), names(condition_colors))
if (length(unknown) > 0) {
  condition_colors <- c(
    condition_colors,
    stats::setNames(grDevices::hcl.colors(length(unknown), "Dark 3"), unknown)
  )
}

plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(Condition, Mean, fill = Condition)
) +
  ggplot2::geom_col(width = 0.66, color = "black", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
    width = 0.2,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    data = point_data,
    ggplot2::aes(Condition, Value),
    inherit.aes = FALSE,
    position = ggplot2::position_jitter(width = 0.08),
    size = 1.3,
    shape = 21,
    fill = "white"
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
    title = "Top cytokines responding to TPA treatment",
    subtitle = "Mean and SD with replicate values",
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
  file.path(output_dir, "Nature_style_cytokine_barplot.pdf"),
  plot,
  width = 8.5,
  height = 5.8,
  useDingbats = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "Nature_style_cytokine_barplot.tiff"),
  plot,
  width = 8.5,
  height = 5.8,
  dpi = 600,
  compression = "lzw",
  device = ragg::agg_tiff
)
ggplot2::ggsave(
  file.path(output_dir, "Nature_style_cytokine_barplot.png"),
  plot,
  width = 8.5,
  height = 5.8,
  dpi = 300
)
utils::write.csv(
  fold_change,
  file.path(output_dir, "cytokine_TPA_to_UT_fold_change.csv"),
  row.names = FALSE
)
cat("Output:", output_dir, "\n")
