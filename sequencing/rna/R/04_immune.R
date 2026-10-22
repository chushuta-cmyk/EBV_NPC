# Run several immune deconvolution methods.

run_immune_module <- function(state, config) {
  required_packages(c("ggplot2", "tidyr", "pheatmap", "RColorBrewer"))
  need_state(state, c(
    "normalized_expression_matrix", "expression_matrix", "sample_info"
  ))

  sample_ids <- rownames(state$sample_info)
  expression <- state$normalized_expression_matrix[, sample_ids, drop = FALSE]
  group_colors <- c(Low_CP = "#4E79A7", High_CP = "#E15759")
  method_scores <- list()
  audit_rows <- list()

  # Capture warnings without stopping other methods.
  run_safely <- function(code) {
    warning_messages <- character()
    value <- tryCatch(
      withCallingHandlers(
        force(code),
        warning = function(warning) {
          warning_messages <<- c(warning_messages, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) error
    )
    list(value = value, warnings = unique(warning_messages))
  }

  # Keep audit notes short.
  short_note <- function(values, limit = 1500) {
    value <- paste(unique(values[nzchar(values)]), collapse = " | ")
    if (nchar(value) > limit) {
      value <- paste0(substr(value, 1, limit), "...")
    }
    value
  }

  add_audit <- function(
    method, status, cells = NA_integer_, matched = NA_integer_,
    signature_total = NA_integer_, note = ""
  ) {
    coverage <- if (
      is.finite(matched) && is.finite(signature_total) && signature_total > 0
    ) {
      100 * matched / signature_total
    } else {
      NA_real_
    }
    audit_rows[[length(audit_rows) + 1]] <<- data.frame(
      method = method,
      status = status,
      input_genes = nrow(expression),
      output_cell_types = cells,
      matched_signature_genes = matched,
      signature_genes = signature_total,
      signature_coverage_percent = coverage,
      note = note,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  # Convert immunedeconv output to a score matrix.
  as_score_matrix <- function(result) {
    result <- as.data.frame(result, check.names = FALSE)
    if (ncol(result) < 2) {
      stop("Immune result has too few columns.", call. = FALSE)
    }
    cell_column <- if ("cell_type" %in% names(result)) "cell_type" else names(result)[1]
    missing_samples <- setdiff(sample_ids, names(result))
    if (length(missing_samples) > 0) {
      stop(
        "Immune result is missing samples: ",
        paste(missing_samples, collapse = ", "),
        call. = FALSE
      )
    }
    scores <- data.matrix(result[, sample_ids, drop = FALSE])
    rownames(scores) <- make.unique(as.character(result[[cell_column]]))
    scores
  }

  # Read the CIBERSORT signature coverage.
  cibersort_signature <- ""
  cibersort_matched <- NA_integer_
  cibersort_total <- NA_integer_
  if (requireNamespace("CIBERSORT", quietly = TRUE)) {
    cibersort_signature <- system.file(
      "extdata", "LM22.txt", package = "CIBERSORT"
    )
    if (nzchar(cibersort_signature)) {
      signature_table <- utils::read.delim(
        cibersort_signature,
        row.names = 1,
        check.names = FALSE
      )
      cibersort_total <- nrow(signature_table)
      cibersort_matched <- length(intersect(
        rownames(expression), rownames(signature_table)
      ))
    }
  }

  # Load or run CIBERSORT as the primary immune method.
  if (nzchar(config$cibersort_file)) {
    if (!file.exists(config$cibersort_file)) {
      stop("RNA_CIBERSORT_FILE does not exist: ", config$cibersort_file, call. = FALSE)
    }
    immune_results <- utils::read.csv(
      config$cibersort_file,
      row.names = 1,
      check.names = FALSE
    )
    if (!all(sample_ids %in% rownames(immune_results))) {
      stop("The CIBERSORT table does not contain all samples.", call. = FALSE)
    }
    immune_results <- immune_results[sample_ids, , drop = FALSE]
    add_audit(
      "CIBERSORT", "loaded", min(22L, ncol(immune_results)),
      cibersort_matched, cibersort_total,
      paste("Precomputed file:", config$cibersort_file)
    )
  } else {
    required_packages("CIBERSORT")
    if (!nzchar(cibersort_signature)) {
      stop("The CIBERSORT LM22 file was not found.", call. = FALSE)
    }

    old_parallel_options <- options(
      parallelly.availableCores.methods = "mc.cores",
      mc.cores = 3L
    )
    cibersort_run <- run_safely(CIBERSORT::cibersort(
      cibersort_signature,
      expression,
      perm = config$cibersort_permutations,
      QN = FALSE
    ))
    options(old_parallel_options)
    if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::sequential)
    }
    if (inherits(cibersort_run$value, "error")) {
      stop(
        "CIBERSORT failed: ", conditionMessage(cibersort_run$value),
        call. = FALSE
      )
    }
    immune_results <- as.data.frame(cibersort_run$value, check.names = FALSE)
    add_audit(
      "CIBERSORT", "completed", min(22L, ncol(immune_results)),
      cibersort_matched, cibersort_total,
      short_note(c(
        paste("Permutations:", config$cibersort_permutations),
        "Quantile normalization: FALSE",
        cibersort_run$warnings
      ))
    )
  }

  if (ncol(immune_results) < 22) {
    stop("CIBERSORT needs at least 22 cell columns.", call. = FALSE)
  }
  cell_data <- data.matrix(immune_results[, seq_len(22), drop = FALSE])
  cell_data <- cell_data[sample_ids, , drop = FALSE]
  method_scores[["CIBERSORT"]] <- t(cell_data)

  # Prepare signature coverage for other methods.
  coverage_info <- list(
    MCPcounter = c(matched = NA_integer_, total = NA_integer_),
    EPIC = c(matched = NA_integer_, total = NA_integer_),
    quanTIseq = c(matched = NA_integer_, total = NA_integer_),
    xCell = c(matched = NA_integer_, total = NA_integer_),
    ESTIMATE = c(matched = NA_integer_, total = NA_integer_)
  )

  if (requireNamespace("EPIC", quietly = TRUE)) {
    epic_data <- new.env(parent = emptyenv())
    utils::data("BRef", package = "EPIC", envir = epic_data)
    epic_genes <- epic_data$BRef$sigGenes
    coverage_info$EPIC <- c(
      matched = length(intersect(rownames(expression), epic_genes)),
      total = length(epic_genes)
    )
  }
  if (requireNamespace("quantiseqr", quietly = TRUE)) {
    quantiseq_file <- system.file(
      "extdata", "TIL10_signature.txt", package = "quantiseqr"
    )
    remove_files <- c(
      system.file("extdata", "TIL10_rmgenes.txt", package = "quantiseqr"),
      system.file(
        "extdata", "TIL10_TCGA_aberrant_immune_genes.txt",
        package = "quantiseqr"
      )
    )
    quantiseq_table <- utils::read.delim(
      quantiseq_file,
      row.names = 1,
      check.names = FALSE
    )
    remove_genes <- unique(unlist(lapply(remove_files, readLines, warn = FALSE)))
    quantiseq_genes <- setdiff(rownames(quantiseq_table), remove_genes)
    coverage_info$quanTIseq <- c(
      matched = length(intersect(rownames(expression), quantiseq_genes)),
      total = length(quantiseq_genes)
    )
  }
  if (requireNamespace("xCell", quietly = TRUE)) {
    xcell_data <- new.env(parent = emptyenv())
    utils::data("xCell.data", package = "xCell", envir = xcell_data)
    xcell_genes <- xcell_data$xCell.data$genes
    coverage_info$xCell <- c(
      matched = length(intersect(rownames(expression), xcell_genes)),
      total = length(xcell_genes)
    )
  }

  # Run methods provided by immunedeconv.
  additional_methods <- c(
    MCPcounter = "mcp_counter",
    EPIC = "epic",
    quanTIseq = "quantiseq",
    xCell = "xcell",
    ESTIMATE = "estimate"
  )
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    for (method in names(additional_methods)) {
      info <- coverage_info[[method]]
      add_audit(
        method, "not available", NA_integer_, info[["matched"]],
        info[["total"]], "immunedeconv is not installed"
      )
    }
  } else {
    for (method in names(additional_methods)) {
      package_name <- switch(
        method,
        MCPcounter = "MCPcounter",
        EPIC = "EPIC",
        quanTIseq = "quantiseqr",
        xCell = "xCell",
        ESTIMATE = "immunedeconv"
      )
      info <- coverage_info[[method]]
      if (!requireNamespace(package_name, quietly = TRUE)) {
        add_audit(
          method, "not available", NA_integer_, info[["matched"]],
          info[["total"]], paste(package_name, "is not installed")
        )
        next
      }

      suppressPackageStartupMessages(
        library(package_name, character.only = TRUE)
      )
      if (method == "xCell") {
        utils::data("xCell.data", package = "xCell", envir = .GlobalEnv)
      }
      if (method == "MCPcounter") {
        method_run <- run_safely(MCPcounter::MCPcounter.estimate(
          expression,
          featuresType = "HUGO_symbols"
        ))
      } else {
        method_run <- run_safely(immunedeconv::deconvolute(
          expression,
          method = additional_methods[[method]],
          tumor = TRUE,
          arrays = FALSE
        ))
      }

      if (inherits(method_run$value, "error")) {
        add_audit(
          method, "failed", NA_integer_, info[["matched"]], info[["total"]],
          short_note(c(conditionMessage(method_run$value), method_run$warnings))
        )
        next
      }

      score_matrix <- if (method == "MCPcounter") {
        as.matrix(method_run$value)
      } else {
        as_score_matrix(method_run$value)
      }
      method_scores[[method]] <- score_matrix
      add_audit(
        method, "completed", nrow(score_matrix), info[["matched"]],
        info[["total"]], short_note(method_run$warnings)
      )
      file_name <- paste0(
        "immune_", tolower(gsub("[^A-Za-z0-9]+", "_", method)), "_scores.csv"
      )
      save_table(
        data.frame(
          cell_type = rownames(score_matrix), score_matrix,
          check.names = FALSE
        ),
        file.path(config$tables_dir, file_name)
      )
    }
  }

  # Build one long table for all successful methods.
  method_long <- do.call(rbind, lapply(names(method_scores), function(method) {
    scores <- method_scores[[method]]
    wide <- data.frame(
      method = method,
      cell_type = rownames(scores),
      scores,
      check.names = FALSE,
      row.names = NULL
    )
    tidyr::pivot_longer(
      wide,
      cols = tidyselect::all_of(sample_ids),
      names_to = "sample",
      values_to = "score"
    )
  }))
  method_long$group <- state$sample_info[method_long$sample, "Group"]

  # Compare groups within each method and cell type.
  split_scores <- split(
    method_long,
    interaction(method_long$method, method_long$cell_type, drop = TRUE)
  )
  comparison_rows <- lapply(split_scores, function(table) {
    low <- table$score[table$group == "Low_CP"]
    high <- table$score[table$group == "High_CP"]
    p_value <- tryCatch(
      stats::wilcox.test(high, low, exact = FALSE)$p.value,
      error = function(error) NA_real_
    )
    data.frame(
      method = table$method[[1]],
      cell_type = table$cell_type[[1]],
      low_mean = mean(low, na.rm = TRUE),
      high_mean = mean(high, na.rm = TRUE),
      high_minus_low = mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE),
      p_value = p_value,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  group_comparison <- do.call(rbind, comparison_rows)
  group_comparison$padj <- ave(
    group_comparison$p_value,
    group_comparison$method,
    FUN = function(value) stats::p.adjust(value, method = "BH")
  )
  group_comparison <- group_comparison[
    order(group_comparison$method, group_comparison$padj),
    ,
    drop = FALSE
  ]

  # Save one heatmap for each successful method.
  annotation <- data.frame(Group = state$sample_info$Group)
  rownames(annotation) <- rownames(state$sample_info)
  for (method in names(method_scores)) {
    scores <- method_scores[[method]][, sample_ids, drop = FALSE]
    safe_name <- tolower(gsub("[^A-Za-z0-9]+", "_", method))
    pheatmap::pheatmap(
      safe_scale_rows(scores),
      annotation_col = annotation,
      annotation_colors = list(Group = group_colors),
      color = grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(11, "RdBu"))
      )(100),
      show_colnames = FALSE,
      border_color = NA,
      filename = file.path(
        config$figures_dir,
        paste0("immune_", safe_name, "_heatmap.png")
      ),
      width = 10,
      height = max(5, min(14, 3 + 0.18 * nrow(scores)))
    )
  }

  # Keep the original CIBERSORT figures.
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

  boxplot <- ggplot2::ggplot(
    immune_long,
    ggplot2::aes(cell_type, proportion, fill = group)
  ) +
    ggplot2::geom_boxplot(width = 0.65, outlier.size = 1) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::labs(
      title = "CIBERSORT immune cell proportions",
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
      title = "CIBERSORT immune cell composition",
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

  pheatmap::pheatmap(
    t(as.matrix(cell_data)),
    annotation_col = annotation,
    annotation_colors = list(Group = group_colors),
    color = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    )(100),
    show_colnames = FALSE,
    border_color = NA,
    filename = file.path(config$figures_dir, "immune_heatmap.png"),
    width = 10,
    height = 8
  )

  immune_audit <- do.call(rbind, audit_rows)
  save_table(
    immune_results,
    file.path(config$tables_dir, "cibersort_results.csv"),
    row.names = TRUE
  )
  save_table(
    immune_long,
    file.path(config$tables_dir, "immune_cell_proportions_long.csv")
  )
  save_table(
    method_long,
    file.path(config$tables_dir, "immune_multimethod_scores_long.csv")
  )
  save_table(
    group_comparison,
    file.path(config$tables_dir, "immune_multimethod_group_comparison.csv")
  )
  save_table(
    immune_audit,
    file.path(config$tables_dir, "immune_method_audit.csv")
  )

  state$immune_results <- immune_results
  state$immune_cells <- cell_data
  state$immune_method_scores <- method_scores
  state$immune_method_audit <- immune_audit
  message(
    "Immune module complete: ", length(method_scores),
    " methods completed."
  )
  invisible(state)
}
