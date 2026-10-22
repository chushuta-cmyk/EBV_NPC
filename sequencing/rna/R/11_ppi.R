# Build updated STRING protein interaction networks.

run_ppi_analysis <- function(state, config) {
  required_packages(c("STRINGdb", "igraph", "ggraph", "ggplot2", "patchwork"))
  need_state(state, "de_results")

  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  de <- as.data.frame(state$de_results)
  de$gene <- rownames(de)
  strict <- !is.na(de$padj) & de$padj < 0.01 & abs(de$log2FoldChange) > 2
  up <- de$gene[strict & de$log2FoldChange > 0]
  down <- de$gene[strict & de$log2FoldChange < 0]

  string_cache <- file.path(config$cache_dir, "stringdb")
  dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)
  string_db <- STRINGdb::STRINGdb$new(
    version = "12.0",
    species = 9606,
    score_threshold = 400,
    input_directory = string_cache
  )

  make_network <- function(genes, direction) {
    mapped <- string_db$map(
      data.frame(gene = genes, stringsAsFactors = FALSE),
      "gene",
      removeUnmappedRows = TRUE,
      takeFirst = TRUE
    )
    network <- string_db$get_subnetwork(mapped$STRING_id)
    if (igraph::vcount(network) == 0) {
      stop("STRING returned an empty network for ", direction, ".", call. = FALSE)
    }
    name_map <- setNames(mapped$gene, mapped$STRING_id)
    igraph::V(network)$label <- unname(name_map[igraph::V(network)$name])
    igraph::V(network)$degree <- igraph::degree(network)
    igraph::V(network)$closeness <- igraph::closeness(
      network,
      normalized = TRUE
    )
    nodes <- data.frame(
      direction = direction,
      string_id = igraph::V(network)$name,
      gene = igraph::V(network)$label,
      degree = igraph::V(network)$degree,
      closeness = igraph::V(network)$closeness,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    edges <- igraph::as_data_frame(network, what = "edges")
    edges$direction <- direction
    list(network = network, nodes = nodes, edges = edges)
  }

  up_network <- make_network(up, "Higher in High_CP")
  down_network <- make_network(down, "Lower in High_CP")

  plot_network <- function(result, title, color) {
    set.seed(42)
    ggraph::ggraph(result$network, layout = "fr") +
      ggraph::geom_edge_link(color = "grey75", alpha = 0.7) +
      ggraph::geom_node_point(
        ggplot2::aes(size = degree),
        color = color,
        alpha = 0.9
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = label),
        repel = TRUE,
        size = 3
      ) +
      ggplot2::scale_size_continuous(range = c(3, 10)) +
      ggplot2::labs(title = title, size = "Degree") +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(legend.position = "bottom")
  }

  figure <- plot_network(
    up_network,
    "Genes higher in High_CP",
    "#D55E00"
  ) + plot_network(
    down_network,
    "Genes lower in High_CP",
    "#0072B2"
  ) + patchwork::plot_layout(guides = "collect")
  ggplot2::ggsave(
    file.path(config$figures_dir, "ppi_updated.png"),
    figure,
    width = 14,
    height = 7,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(config$figures_dir, "ppi_updated.pdf"),
    figure,
    width = 14,
    height = 7
  )

  save_table(
    rbind(up_network$nodes, down_network$nodes),
    file.path(config$tables_dir, "ppi_node_centrality.csv")
  )
  save_table(
    rbind(up_network$edges, down_network$edges),
    file.path(config$tables_dir, "ppi_edges.csv")
  )
  invisible(list(up = up_network, down = down_network))
}
