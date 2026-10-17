# Load counts and run differential expression.

run_data_module <- function(state, config) {
  required_packages("DESeq2")

  low_ids <- c(
    "GSM7024418", "GSM7024402", "GSM7024424", "GSM7024411", "GSM7024416"
  )
  high_ids <- c(
    "GSM7024387", "GSM7024389", "GSM7024391", "GSM7024392", "GSM7024394",
    "GSM7024396", "GSM7024398", "GSM7024404", "GSM7024407", "GSM7024413",
    "GSM7024414", "GSM7024422", "GSM7024385", "GSM7024400", "GSM7024420",
    "GSM7024409", "GSM7024405"
  )

  if (nzchar(config$counts_file)) {
    if (!file.exists(config$counts_file)) {
      stop("RNA_COUNTS_FILE does not exist: ", config$counts_file, call. = FALSE)
    }
    counts <- utils::read.csv(
      config$counts_file,
      row.names = 1,
      check.names = FALSE
    )
    counts <- as.matrix(counts)
    storage.mode(counts) <- "numeric"
    pheno_data <- NULL
  } else {
    required_packages("GEOquery")
    geo <- GEOquery::getGEO(
      "GSE224450",
      GSEMatrix = TRUE,
      destdir = config$cache_dir
    )
    counts <- Biobase::exprs(geo[[1]])
    pheno_data <- Biobase::pData(geo[[1]])
  }

  sample_ids <- c(low_ids, high_ids)
  missing_samples <- setdiff(sample_ids, colnames(counts))
  if (length(missing_samples) > 0) {
    stop(
      "The count matrix is missing samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  counts <- counts[, sample_ids, drop = FALSE]
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("The count matrix must contain finite non-negative values.", call. = FALSE)
  }
  counts <- round(counts)

  sample_info <- data.frame(
    Group = factor(
      c(rep("Low_CP", length(low_ids)), rep("High_CP", length(high_ids))),
      levels = c("Low_CP", "High_CP")
    ),
    row.names = sample_ids
  )

  keep <- rowSums(counts >= 10) >= 2
  filtered_counts <- counts[keep, , drop = FALSE]
  if (nrow(filtered_counts) == 0) {
    stop("No genes passed the count filter.", call. = FALSE)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = sample_info,
    design = ~Group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  de_results <- DESeq2::results(dds, contrast = c("Group", "High_CP", "Low_CP"))
  rld <- DESeq2::rlog(dds, blind = FALSE)
  expression_matrix <- SummarizedExperiment::assay(rld)

  state$raw_counts <- counts
  state$filtered_counts <- filtered_counts
  state$sample_info <- sample_info
  state$pheno_data <- pheno_data
  state$dds <- dds
  state$de_results <- de_results
  state$expression_matrix <- expression_matrix

  de_table <- as.data.frame(de_results)
  de_table$gene <- rownames(de_table)
  de_table <- de_table[, c("gene", setdiff(names(de_table), "gene"))]
  save_table(de_table, file.path(config$tables_dir, "differential_expression.csv"))
  save_table(
    data.frame(sample = rownames(sample_info), group = sample_info$Group),
    file.path(config$tables_dir, "sample_groups.csv")
  )
  saveRDS(state, file.path(config$cache_dir, "data_state.rds"))

  message("Data module complete: ", nrow(filtered_counts), " genes and ", ncol(filtered_counts), " samples.")
  invisible(state)
}
