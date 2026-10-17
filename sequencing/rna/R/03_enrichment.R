# Run GO and KEGG enrichment for changed genes.

run_enrichment_module <- function(state, config) {
  required_packages(c("clusterProfiler", "org.Hs.eg.db", "ggplot2"))
  need_state(state, "de_results")

  de_table <- as.data.frame(state$de_results)
  de_table$gene <- rownames(de_table)
  up_genes <- de_table$gene[
    !is.na(de_table$padj) & de_table$padj < 0.01 & de_table$log2FoldChange > 2
  ]
  down_genes <- de_table$gene[
    !is.na(de_table$padj) & de_table$padj < 0.01 & de_table$log2FoldChange < -2
  ]

  convert_ids <- function(genes) {
    if (length(genes) == 0) return(character())
    mapped <- clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
    unique(mapped$ENTREZID)
  }

  plot_enrichment <- function(table, title, path) {
    if (nrow(table) == 0) return(invisible(NULL))
    shown <- head(table[order(table$p.adjust, -table$Count), ], 15)
    shown$GeneRatioNumber <- ratio_to_number(shown$GeneRatio)
    shown$Description <- factor(shown$Description, levels = rev(shown$Description))
    plot <- ggplot2::ggplot(
      shown,
      ggplot2::aes(GeneRatioNumber, Description, size = Count, color = p.adjust)
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(title = title, x = "Gene ratio", y = NULL) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(path, plot, width = 10, height = 7, dpi = 300)
  }

  run_group <- function(genes, label) {
    ids <- convert_ids(genes)
    if (length(ids) == 0) {
      warning("No ENTREZ IDs were found for ", label, " genes.")
      return(list(go = data.frame(), kegg = data.frame()))
    }

    go <- clusterProfiler::enrichGO(
      gene = ids,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    go_table <- as.data.frame(go)

    kegg_table <- tryCatch({
      as.data.frame(clusterProfiler::enrichKEGG(
        gene = ids,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ))
    }, error = function(error) {
      warning("KEGG enrichment failed for ", label, ": ", conditionMessage(error))
      data.frame()
    })

    save_table(go_table, file.path(config$tables_dir, paste0("go_", label, ".csv")))
    save_table(kegg_table, file.path(config$tables_dir, paste0("kegg_", label, ".csv")))
    plot_enrichment(
      go_table,
      paste("GO enrichment:", label),
      file.path(config$figures_dir, paste0("go_", label, ".png"))
    )
    plot_enrichment(
      kegg_table,
      paste("KEGG enrichment:", label),
      file.path(config$figures_dir, paste0("kegg_", label, ".png"))
    )
    list(go = go_table, kegg = kegg_table)
  }

  state$upregulated_genes <- up_genes
  state$downregulated_genes <- down_genes
  state$enrichment <- list(
    up = run_group(up_genes, "upregulated"),
    down = run_group(down_genes, "downregulated")
  )

  message("Enrichment module complete.")
  invisible(state)
}
