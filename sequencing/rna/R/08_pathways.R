# Compare pathway scores with immune cells.

run_pathways_module <- function(state, config) {
  required_packages(c("ggplot2", "pheatmap", "org.Hs.eg.db"))
  need_state(state, c("expression_matrix", "immune_cells", "sample_info"))

  expression <- state$expression_matrix
  immune <- state$immune_cells[colnames(expression), , drop = FALSE]

  score_pathway <- function(genes) {
    genes <- intersect(genes, rownames(expression))
    if (length(genes) == 0) {
      stop("No pathway genes are present in the expression matrix.", call. = FALSE)
    }
    colMeans(expression[genes, , drop = FALSE], na.rm = TRUE)
  }

  correlation_table <- function(score) {
    rows <- lapply(colnames(immune), function(cell) {
      test <- stats::cor.test(immune[, cell], score, use = "complete.obs")
      data.frame(
        immune_cell = cell,
        correlation = unname(test$estimate),
        p_value = test$p.value,
        row.names = NULL
      )
    })
    result <- do.call(rbind, rows)
    result$adjusted_p_value <- stats::p.adjust(result$p_value, method = "BH")
    result[order(result$correlation), ]
  }

  plot_correlation <- function(table, title, path) {
    plot <- ggplot2::ggplot(
      table,
      ggplot2::aes(reorder(immune_cell, correlation), correlation, fill = correlation)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
      ggplot2::labs(title = title, x = NULL, y = "Pearson correlation") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(path, plot, width = 9, height = 7, dpi = 300)
  }

  immune_genes <- c(
    "CD3D", "CD3E", "CD8A", "CD4", "HLA-DRA", "HLA-DRB1",
    "GZMB", "PRF1", "IFNG", "IL2RA"
  )
  immune_score <- score_pathway(immune_genes)
  immune_correlation <- correlation_table(immune_score)
  save_table(
    immune_correlation,
    file.path(config$tables_dir, "immune_activation_correlation.csv")
  )
  plot_correlation(
    immune_correlation,
    "Immune activation and immune cells",
    file.path(config$figures_dir, "immune_activation_correlation.png")
  )

  ebv_entrez <- c("596", "919", "915", "916", "917", "3109", "3113", "3115")
  ebv_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ebv_entrez,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  ebv_score <- score_pathway(stats::na.omit(unname(ebv_symbols)))
  ebv_correlation <- correlation_table(ebv_score)
  save_table(
    ebv_correlation,
    file.path(config$tables_dir, "ebv_pathway_correlation.csv")
  )
  plot_correlation(
    ebv_correlation,
    "EBV pathway and immune cells",
    file.path(config$figures_dir, "ebv_pathway_correlation.png")
  )

  neutrophil_genes <- intersect(
    c("S100A8", "S100A9", "CXCR2", "CSF3R", "CXCL1", "CXCL2"),
    rownames(expression)
  )
  if (length(neutrophil_genes) > 0) {
    neutrophil_correlation <- stats::cor(
      t(expression[neutrophil_genes, , drop = FALSE]),
      immune,
      use = "pairwise.complete.obs"
    )
    save_table(
      neutrophil_correlation,
      file.path(config$tables_dir, "neutrophil_gene_immune_correlation.csv"),
      row.names = TRUE
    )
    pheatmap::pheatmap(
      neutrophil_correlation,
      color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-1, 1, length.out = 101),
      border_color = NA,
      filename = file.path(config$figures_dir, "neutrophil_gene_immune_correlation.png"),
      width = 10,
      height = 6
    )
  }

  up_go <- c(
    "KLRK1", "SH2D1A", "BCL2", "CCL19", "CCL21", "VCAM1", "CD5",
    "CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1"
  )
  up_kegg <- c(
    "CD3D", "CD3G", "CD19", "CD34", "HLA-DQB1", "CD3E", "CD247",
    "LCK", "IL2RG", "HLA-DPB1", "MS4A1"
  )
  down_go <- c(
    "CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R",
    "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX"
  )
  down_kegg <- down_go
  overlap <- rbind(
    data.frame(direction = "up", gene = intersect(up_go, up_kegg)),
    data.frame(direction = "down", gene = intersect(down_go, down_kegg))
  )
  save_table(overlap, file.path(config$tables_dir, "pathway_gene_overlap.csv"))

  state$immune_activation_score <- immune_score
  state$ebv_pathway_score <- ebv_score
  state$pathway_overlap <- overlap
  message("Pathway module complete.")
  invisible(state)
}
