# Load or calculate immune cell proportions.

run_immune_module <- function(state, config) {
  required_packages(c("ggplot2", "tidyr", "pheatmap", "RColorBrewer"))
  need_state(state, c("expression_matrix", "sample_info"))

  if (nzchar(config$cibersort_file)) {
    if (!file.exists(config$cibersort_file)) {
      stop("RNA_CIBERSORT_FILE does not exist: ", config$cibersort_file, call. = FALSE)
    }
    immune_results <- utils::read.csv(
      config$cibersort_file,
      row.names = 1,
      check.names = FALSE
    )
  } else {
    required_packages("CIBERSORT")
    signature <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
    if (!nzchar(signature)) {
      stop("The CIBERSORT LM22 signature file was not found.", call. = FALSE)
    }
    immune_results <- CIBERSORT::cibersort(
      signature,
      state$expression_matrix,
      perm = 0,
      QN = TRUE
    )
    immune_results <- as.data.frame(immune_results)
  }

  sample_ids <- rownames(state$sample_info)
  if (!all(sample_ids %in% rownames(immune_results))) {
    stop("The immune results do not contain all analysis samples.", call. = FALSE)
  }
  immune_results <- immune_results[sample_ids, , drop = FALSE]
  if (ncol(immune_results) < 22) {
    stop("The immune results need at least 22 cell columns.", call. = FALSE)
  }

  cell_data <- data.matrix(immune_results[, seq_len(22), drop = FALSE])
  state$immune_results <- immune_results
  state$immune_cells <- cell_data

  immune_long <- data.frame(
    sample = rownames(cell_data),
    group = state$sample_info[rownames(cell_data), "Group"],
    cell_data,
    check.names = FALSE
  )
  immune_long <- tidyr::pivot_longer(
    immune_long,
    cols = -c(sample, group),
    names_to = "cell_type",
    values_to = "proportion"
  )

  group_colors <- c(Low_CP = "#4E79A7", High_CP = "#E15759")
  boxplot <- ggplot2::ggplot(
    immune_long,
    ggplot2::aes(cell_type, proportion, fill = group)
  ) +
    ggplot2::geom_boxplot(width = 0.65, outlier.size = 1) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::labs(
      title = "Immune cell proportions",
      x = NULL,
      y = "Proportion",
      fill = "Group"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(
    file.path(config$figures_dir, "immune_boxplot.png"),
    boxplot,
    width = 12,
    height = 7,
    dpi = 300
  )

  stacked <- ggplot2::ggplot(
    immune_long,
    ggplot2::aes(sample, proportion, fill = cell_type)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::labs(
      title = "Immune cell composition",
      x = NULL,
      y = "Proportion",
      fill = "Cell type"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  ggplot2::ggsave(
    file.path(config$figures_dir, "immune_composition.png"),
    stacked,
    width = 13,
    height = 8,
    dpi = 300
  )

  annotation <- data.frame(Group = state$sample_info$Group)
  rownames(annotation) <- rownames(state$sample_info)
  pheatmap::pheatmap(
    t(as.matrix(cell_data)),
    annotation_col = annotation,
    annotation_colors = list(Group = group_colors),
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    show_colnames = FALSE,
    border_color = NA,
    filename = file.path(config$figures_dir, "immune_heatmap.png"),
    width = 10,
    height = 8
  )

  save_table(immune_results, file.path(config$tables_dir, "cibersort_results.csv"), row.names = TRUE)
  save_table(immune_long, file.path(config$tables_dir, "immune_cell_proportions_long.csv"))
  message("Immune module complete.")
  invisible(state)
}
