# Build a gene co-expression network.

run_wgcna_module <- function(state, config) {
  required_packages("WGCNA")
  suppressPackageStartupMessages(library(WGCNA))
  need_state(state, c("expression_matrix", "sample_info"))

  expression <- t(state$expression_matrix)
  storage.mode(expression) <- "numeric"
  good <- WGCNA::goodSamplesGenes(expression, verbose = 0)
  if (!good$allOK) {
    expression <- expression[good$goodSamples, good$goodGenes, drop = FALSE]
  }
  if (ncol(expression) < 30) {
    stop("WGCNA needs at least 30 usable genes.", call. = FALSE)
  }

  powers <- 1:20
  soft_fit <- WGCNA::pickSoftThreshold(
    expression,
    powerVector = powers,
    networkType = "unsigned",
    verbose = 0
  )
  power <- soft_fit$powerEstimate
  if (!is.finite(power)) power <- 6

  grDevices::png(
    file.path(config$figures_dir, "wgcna_soft_threshold.png"),
    width = 1800,
    height = 900,
    res = 150
  )
  graphics::par(mfrow = c(1, 2))
  graphics::plot(
    soft_fit$fitIndices[, 1],
    -sign(soft_fit$fitIndices[, 3]) * soft_fit$fitIndices[, 2],
    xlab = "Power",
    ylab = "Scale-free fit",
    type = "n"
  )
  graphics::text(
    soft_fit$fitIndices[, 1],
    -sign(soft_fit$fitIndices[, 3]) * soft_fit$fitIndices[, 2],
    labels = powers,
    col = "red"
  )
  graphics::abline(h = 0.9, col = "grey40", lty = 2)
  graphics::plot(
    soft_fit$fitIndices[, 1],
    soft_fit$fitIndices[, 5],
    xlab = "Power",
    ylab = "Mean connectivity",
    type = "b"
  )
  grDevices::dev.off()

  network <- WGCNA::blockwiseModules(
    expression,
    power = power,
    networkType = "unsigned",
    TOMType = "unsigned",
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
  module_colors <- WGCNA::labels2colors(network$colors)

  grDevices::png(
    file.path(config$figures_dir, "wgcna_modules.png"),
    width = 1800,
    height = 1200,
    res = 150
  )
  WGCNA::plotDendroAndColors(
    network$dendrograms[[1]],
    module_colors[network$blockGenes[[1]]],
    "Module",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  grDevices::dev.off()

  trait <- data.frame(
    EBV_CP = as.numeric(state$sample_info[rownames(expression), "Group"] == "High_CP"),
    row.names = rownames(expression)
  )
  eigengenes <- WGCNA::orderMEs(network$MEs)
  module_correlation <- stats::cor(eigengenes, trait, use = "pairwise.complete.obs")
  module_p <- WGCNA::corPvalueStudent(module_correlation, nrow(expression))

  grDevices::png(
    file.path(config$figures_dir, "wgcna_module_trait.png"),
    width = 1200,
    height = 1200,
    res = 150
  )
  labels <- paste0(signif(module_correlation, 2), "\n(", signif(module_p, 2), ")")
  dim(labels) <- dim(module_correlation)
  WGCNA::labeledHeatmap(
    Matrix = module_correlation,
    xLabels = colnames(trait),
    yLabels = rownames(module_correlation),
    ySymbols = rownames(module_correlation),
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = labels,
    zlim = c(-1, 1),
    main = "Module-trait relationships"
  )
  grDevices::dev.off()

  module_table <- data.frame(
    gene = colnames(expression),
    module_number = network$colors,
    module_color = module_colors
  )
  trait_table <- data.frame(
    module = rownames(module_correlation),
    correlation = module_correlation[, 1],
    p_value = module_p[, 1],
    row.names = NULL
  )
  save_table(module_table, file.path(config$tables_dir, "wgcna_modules.csv"))
  save_table(trait_table, file.path(config$tables_dir, "wgcna_module_trait.csv"))

  state$wgcna_network <- network
  state$wgcna_trait <- trait_table
  message("WGCNA module complete. Soft power: ", power, ".")
  invisible(state)
}
