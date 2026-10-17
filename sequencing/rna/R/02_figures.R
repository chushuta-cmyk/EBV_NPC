# Create quality-control and expression figures.

run_figures_module <- function(state, config) {
  required_packages(c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))
  need_state(state, c("raw_counts", "expression_matrix", "sample_info", "de_results"))

  sample_colors <- c(Low_CP = "#4DBBD5", High_CP = "#E64B35")
  plot_colors <- unname(sample_colors[as.character(state$sample_info$Group)])

  grDevices::png(
    file.path(config$figures_dir, "normalization_boxplots.png"),
    width = 2400,
    height = 1200,
    res = 200
  )
  graphics::par(mfrow = c(1, 2), mar = c(8, 5, 3, 1))
  graphics::boxplot(
    state$raw_counts,
    col = plot_colors,
    las = 2,
    main = "Before normalization",
    ylab = "Count"
  )
  graphics::boxplot(
    state$expression_matrix,
    col = plot_colors,
    las = 2,
    main = "After normalization",
    ylab = "rlog expression"
  )
  grDevices::dev.off()

  variable_genes <- apply(state$expression_matrix, 1, stats::var) > 0
  pca <- stats::prcomp(
    t(state$expression_matrix[variable_genes, , drop = FALSE]),
    scale. = TRUE
  )
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_data <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = state$sample_info[rownames(pca$x), "Group"]
  )
  pca_plot <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(PC1, PC2, color = group)
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3, max.overlaps = Inf) +
    ggplot2::stat_ellipse(ggplot2::aes(fill = group), geom = "polygon", alpha = 0.12) +
    ggplot2::scale_color_manual(values = sample_colors) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::labs(
      title = "PCA of normalized expression",
      x = sprintf("PC1 (%.1f%%)", explained[[1]]),
      y = sprintf("PC2 (%.1f%%)", explained[[2]]),
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::theme_bw(base_size = 12)
  ggplot2::ggsave(
    file.path(config$figures_dir, "pca.png"),
    pca_plot,
    width = 9,
    height = 7,
    dpi = 300
  )

  volcano <- as.data.frame(state$de_results)
  volcano$gene <- rownames(volcano)
  volcano$status <- "Not significant"
  volcano$status[!is.na(volcano$padj) & volcano$padj < 0.01 & volcano$log2FoldChange > 2] <- "Up"
  volcano$status[!is.na(volcano$padj) & volcano$padj < 0.01 & volcano$log2FoldChange < -2] <- "Down"
  volcano$score <- -log10(pmax(volcano$padj, .Machine$double.xmin))

  labeled <- rbind(
    head(volcano[order(-volcano$log2FoldChange), ], 10),
    head(volcano[order(volcano$log2FoldChange), ], 10)
  )
  volcano_plot <- ggplot2::ggplot(
    volcano,
    ggplot2::aes(log2FoldChange, score, color = status)
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = labeled,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    ggplot2::scale_color_manual(
      values = c(Down = "#4DBBD5", `Not significant` = "grey70", Up = "#E64B35")
    ) +
    ggplot2::labs(
      title = "Differential expression",
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(config$figures_dir, "volcano.png"),
    volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  ranked <- volcano[order(volcano$padj, -abs(volcano$log2FoldChange), na.last = NA), ]
  heatmap_genes <- head(ranked$gene, 50)
  if (length(heatmap_genes) > 1) {
    annotation <- data.frame(Group = state$sample_info$Group)
    rownames(annotation) <- rownames(state$sample_info)
    pheatmap::pheatmap(
      safe_scale_rows(state$expression_matrix[heatmap_genes, , drop = FALSE]),
      annotation_col = annotation,
      annotation_colors = list(Group = sample_colors),
      color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      show_colnames = FALSE,
      fontsize_row = 7,
      border_color = NA,
      filename = file.path(config$figures_dir, "differential_expression_heatmap.png"),
      width = 9,
      height = 9
    )
  }

  save_table(volcano, file.path(config$tables_dir, "volcano_data.csv"))
  message("Figures module complete.")
  invisible(state)
}
