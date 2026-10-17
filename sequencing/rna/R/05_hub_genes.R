# Compare hub genes with immune cell proportions.

run_hub_module <- function(state, config) {
  required_packages(c("ggplot2", "pheatmap", "randomForest"))
  need_state(state, c("expression_matrix", "sample_info", "immune_cells"))

  up_hub_genes <- c("CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
  down_hub_genes <- c(
    "CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R",
    "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX"
  )
  hub_genes <- c(up_hub_genes, down_hub_genes)
  hub_genes <- hub_genes[hub_genes %in% rownames(state$expression_matrix)]
  if (length(hub_genes) == 0) {
    stop("None of the hub genes are present in the expression matrix.", call. = FALSE)
  }

  expression <- state$expression_matrix[hub_genes, , drop = FALSE]
  immune <- state$immune_cells[colnames(expression), , drop = FALSE]
  correlation <- stats::cor(t(expression), immune, use = "pairwise.complete.obs")

  importance <- lapply(hub_genes, function(gene) {
    model <- randomForest::randomForest(
      x = immune,
      y = as.numeric(expression[gene, ]),
      importance = TRUE,
      ntree = 1000
    )
    scores <- randomForest::importance(model)
    data.frame(
      gene = gene,
      immune_cell = rownames(scores),
      importance = scores[, "%IncMSE"],
      row.names = NULL
    )
  })
  importance <- do.call(rbind, importance)

  annotation <- data.frame(
    Group = ifelse(colnames(correlation) %in% colnames(immune), "Immune cell", "Other"),
    row.names = colnames(correlation)
  )
  pheatmap::pheatmap(
    correlation,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    border_color = NA,
    filename = file.path(config$figures_dir, "hub_gene_immune_correlation.png"),
    width = 11,
    height = 8
  )

  importance_plot <- ggplot2::ggplot(
    importance,
    ggplot2::aes(importance, reorder(immune_cell, importance), color = gene)
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::labs(
      title = "Immune cell importance for hub genes",
      x = "Increase in mean squared error",
      y = "Immune cell"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(
    file.path(config$figures_dir, "hub_gene_immune_importance.png"),
    importance_plot,
    width = 11,
    height = 8,
    dpi = 300
  )

  expression_table <- data.frame(
    sample = colnames(expression),
    t(expression),
    check.names = FALSE,
    row.names = NULL
  )
  save_table(as.data.frame(correlation), file.path(config$tables_dir, "hub_gene_immune_correlation.csv"), row.names = TRUE)
  save_table(importance, file.path(config$tables_dir, "hub_gene_immune_importance.csv"))
  save_table(expression_table, file.path(config$tables_dir, "hub_gene_expression.csv"))

  state$hub_genes <- hub_genes
  state$hub_expression <- expression
  state$hub_immune_correlation <- correlation
  message("Hub gene module complete.")
  invisible(state)
}
