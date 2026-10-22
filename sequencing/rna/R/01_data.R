# Read raw RCC counts and run DESeq2.

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
  sample_ids <- c(low_ids, high_ids)

  if (!file.exists(config$rcc_tar_file)) {
    stop("Raw RCC archive does not exist: ", config$rcc_tar_file, call. = FALSE)
  }

  # Read one RCC file.
  read_rcc <- function(path) {
    lines <- readLines(gzfile(path), warn = FALSE)
    start <- which(trimws(lines) == "<Code_Summary>")
    end <- which(trimws(lines) == "</Code_Summary>")
    if (length(start) != 1 || length(end) != 1 || start >= end) {
      stop("Code_Summary block not found in ", path, call. = FALSE)
    }
    block <- lines[(start + 1):(end - 1)]
    block <- block[nzchar(trimws(block))]
    table <- utils::read.csv(
      text = paste(block, collapse = "\n"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    required <- c("Name", "CodeClass", "Count")
    if (!all(required %in% names(table))) {
      stop("Required RCC columns are missing in ", path, call. = FALSE)
    }
    table$Count <- as.numeric(table$Count)
    table
  }

  # Extract the selected RCC files.
  archive_members <- utils::untar(config$rcc_tar_file, list = TRUE)
  selected_members <- vapply(sample_ids, function(sample_id) {
    hits <- archive_members[
      startsWith(basename(archive_members), paste0(sample_id, "_")) &
        grepl("\\.RCC\\.gz$", archive_members, ignore.case = TRUE)
    ]
    if (length(hits) != 1) {
      stop(
        "Expected one RCC file for ", sample_id, ", found ", length(hits),
        call. = FALSE
      )
    }
    hits
  }, character(1))

  extract_dir <- tempfile("nanostring_rcc_")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(
    config$rcc_tar_file,
    files = unname(selected_members),
    exdir = extract_dir
  )

  rcc_tables <- lapply(selected_members, function(member) {
    read_rcc(file.path(extract_dir, member))
  })
  names(rcc_tables) <- names(selected_members)

  # Check probe names and classes.
  reference_names <- rcc_tables[[1]]$Name
  reference_classes <- rcc_tables[[1]]$CodeClass
  if (anyDuplicated(reference_names)) {
    stop("Raw RCC probe names are not unique.", call. = FALSE)
  }
  for (sample_id in sample_ids) {
    if (!identical(rcc_tables[[sample_id]]$Name, reference_names)) {
      stop("Probe order differs in ", sample_id, call. = FALSE)
    }
    if (!identical(rcc_tables[[sample_id]]$CodeClass, reference_classes)) {
      stop("Probe classes differ in ", sample_id, call. = FALSE)
    }
  }

  probe_class <- setNames(reference_classes, reference_names)
  raw_counts <- do.call(cbind, lapply(sample_ids, function(sample_id) {
    rcc_tables[[sample_id]]$Count
  }))
  rownames(raw_counts) <- reference_names
  colnames(raw_counts) <- sample_ids

  # Check raw counts without rounding.
  if (any(!is.finite(raw_counts)) || any(raw_counts < 0)) {
    stop("Raw RCC counts must be finite and non-negative.", call. = FALSE)
  }
  if (any(raw_counts != round(raw_counts))) {
    stop("Raw RCC counts must be integers.", call. = FALSE)
  }
  storage.mode(raw_counts) <- "integer"

  sample_info <- data.frame(
    Group = factor(
      c(rep("Low_CP", length(low_ids)), rep("High_CP", length(high_ids))),
      levels = c("Low_CP", "High_CP")
    ),
    row.names = sample_ids
  )

  # Keep endogenous probes that pass the original filter.
  endogenous <- names(probe_class)[probe_class == "Endogenous"]
  housekeeping <- names(probe_class)[probe_class == "Housekeeping"]
  keep <- rowSums(raw_counts[endogenous, , drop = FALSE] >= 10) >= 2
  filtered_genes <- endogenous[keep]
  if (length(filtered_genes) == 0) {
    stop("No endogenous probes passed the count filter.", call. = FALSE)
  }
  if (length(housekeeping) == 0) {
    stop("No housekeeping probes were found.", call. = FALSE)
  }

  # Use housekeeping probes for DESeq2 size factors.
  model_genes <- unique(c(filtered_genes, housekeeping))
  model_counts <- raw_counts[model_genes, sample_ids, drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = model_counts,
    colData = sample_info,
    design = ~ Group
  )
  housekeeping_index <- which(rownames(dds) %in% housekeeping)
  dds <- DESeq2::estimateSizeFactors(dds, controlGenes = housekeeping_index)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # Keep only endogenous probes in the reported results.
  all_results <- DESeq2::results(
    dds,
    contrast = c("Group", "High_CP", "Low_CP")
  )
  de_results <- all_results[filtered_genes, ]

  # Make normalized counts for immune analysis.
  normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  normalized_matrix <- normalized_matrix[filtered_genes, sample_ids, drop = FALSE]

  # Make rlog values only after differential analysis.
  rld <- DESeq2::rlog(dds, blind = FALSE)
  rlog_matrix <- SummarizedExperiment::assay(rld)
  expression_matrix <- rlog_matrix[filtered_genes, sample_ids, drop = FALSE]

  state$raw_counts <- raw_counts[filtered_genes, sample_ids, drop = FALSE]
  state$filtered_counts <- state$raw_counts
  state$sample_info <- sample_info
  state$pheno_data <- NULL
  state$probe_class <- probe_class
  state$housekeeping_genes <- housekeeping
  state$dds <- dds
  state$de_results <- de_results
  state$expression_matrix <- expression_matrix
  state$normalized_expression_matrix <- normalized_matrix

  de_table <- as.data.frame(de_results)
  de_table$gene <- rownames(de_table)
  de_table <- de_table[, c("gene", setdiff(names(de_table), "gene"))]
  audit <- data.frame(
    sample = sample_ids,
    group = as.character(sample_info$Group),
    raw_library_sum = colSums(raw_counts[, sample_ids, drop = FALSE]),
    raw_housekeeping_sum = colSums(
      raw_counts[housekeeping, sample_ids, drop = FALSE]
    ),
    deseq2_housekeeping_size_factor = DESeq2::sizeFactors(dds)[sample_ids],
    row.names = NULL,
    check.names = FALSE
  )
  probe_audit <- data.frame(
    metric = c(
      "Raw RCC probes", "Endogenous probes", "Housekeeping probes",
      "Filtered endogenous probes", "Low_CP samples", "High_CP samples",
      "Raw RCC values are integers", "Rounding applied"
    ),
    value = c(
      nrow(raw_counts), length(endogenous), length(housekeeping),
      length(filtered_genes), length(low_ids), length(high_ids), "TRUE", "FALSE"
    ),
    check.names = FALSE
  )

  save_table(de_table, file.path(config$tables_dir, "differential_expression.csv"))
  save_table(
    data.frame(gene = rownames(state$raw_counts), state$raw_counts, check.names = FALSE),
    file.path(config$tables_dir, "filtered_raw_rcc_counts.csv")
  )
  save_table(
    data.frame(gene = rownames(expression_matrix), expression_matrix, check.names = FALSE),
    file.path(config$tables_dir, "rlog_expression.csv")
  )
  save_table(
    data.frame(gene = rownames(normalized_matrix), normalized_matrix, check.names = FALSE),
    file.path(config$tables_dir, "deseq2_normalized_counts.csv")
  )
  save_table(audit, file.path(config$tables_dir, "raw_rcc_normalization_audit.csv"))
  save_table(probe_audit, file.path(config$tables_dir, "raw_rcc_probe_audit.csv"))
  save_table(
    data.frame(sample = rownames(sample_info), group = sample_info$Group),
    file.path(config$tables_dir, "sample_groups.csv")
  )
  saveRDS(state, file.path(config$cache_dir, "data_state.rds"))

  message(
    "Data module complete: ", length(filtered_genes),
    " endogenous probes and ", length(sample_ids), " samples."
  )
  invisible(state)
}
