# Run metadata and differential-expression sensitivity checks.

run_sensitivity_checks <- function(state, config, series_matrix_file) {
  required_packages(c("DESeq2", "GEOquery"))
  need_state(state, c("dds", "housekeeping_genes", "sample_info", "de_results"))

  if (!file.exists(series_matrix_file)) {
    stop("The GEO series matrix file was not found.", call. = FALSE)
  }

  series <- GEOquery::getGEO(filename = series_matrix_file)
  metadata <- Biobase::pData(series)
  sample_ids <- rownames(state$sample_info)
  metadata <- metadata[sample_ids, , drop = FALSE]

  clean_value <- function(value) {
    sub("^[^:]+:[[:space:]]*", "", as.character(value))
  }
  clinical <- data.frame(
    sample = sample_ids,
    group = as.character(state$sample_info$Group),
    age = as.numeric(metadata[sample_ids, "age:ch1"]),
    sex = metadata[sample_ids, "Sex:ch1"],
    stage = metadata[sample_ids, "staging (ajcc):ch1"],
    metastasis = metadata[sample_ids, "metastasis?:ch1"],
    ebv_dna_record = metadata[sample_ids, "baseline ebv dna (copies/ml):ch1"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  clinical$sex <- clean_value(clinical$sex)
  clinical$stage <- clean_value(clinical$stage)
  clinical$metastasis <- clean_value(clinical$metastasis)
  clinical$ebv_dna_record <- clean_value(clinical$ebv_dna_record)
  clinical$stage4c <- ifelse(
    clinical$stage == "4C", "Stage4C",
    ifelse(clinical$stage == "NA", NA_character_, "OtherStage")
  )
  save_table(clinical, file.path(config$tables_dir, "sample_clinical_metadata.csv"))

  model_counts <- DESeq2::counts(state$dds, normalized = FALSE)
  housekeeping_index <- which(rownames(model_counts) %in% state$housekeeping_genes)
  filtered_genes <- rownames(state$de_results)

  run_model <- function(samples, design, label) {
    model_metadata <- clinical[match(samples, clinical$sample), , drop = FALSE]
    rownames(model_metadata) <- model_metadata$sample
    model_metadata$Group <- factor(
      model_metadata$group,
      levels = c("Low_CP", "High_CP")
    )
    model_metadata$Stage4C <- factor(
      model_metadata$stage4c,
      levels = c("OtherStage", "Stage4C")
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = model_counts[, samples, drop = FALSE],
      colData = model_metadata,
      design = design
    )
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = housekeeping_index)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    result <- DESeq2::results(
      dds,
      contrast = c("Group", "High_CP", "Low_CP")
    )
    result <- as.data.frame(result[filtered_genes, ])
    result$gene <- rownames(result)
    result$analysis <- label
    result <- result[, c("analysis", "gene", setdiff(names(result), c("analysis", "gene")))]
    result
  }

  known_samples <- clinical$sample[clinical$ebv_dna_record != "NA"]
  stage_samples <- clinical$sample[
    clinical$ebv_dna_record != "NA" & !is.na(clinical$stage4c)
  ]
  known_result <- run_model(known_samples, ~ Group, "Exclude missing EBV DNA")
  stage_result <- run_model(
    stage_samples,
    ~ Stage4C + Group,
    "Exclude missing EBV DNA and adjust for stage 4C"
  )
  sensitivity <- rbind(known_result, stage_result)
  save_table(
    sensitivity,
    file.path(config$tables_dir, "differential_expression_sensitivity.csv")
  )

  primary <- as.data.frame(state$de_results)
  primary$gene <- rownames(primary)
  strict_count <- function(table) {
    keep <- !is.na(table$padj) & table$padj < 0.01 & abs(table$log2FoldChange) > 2
    c(
      total = sum(keep),
      up = sum(keep & table$log2FoldChange > 0),
      down = sum(keep & table$log2FoldChange < 0)
    )
  }
  summary <- rbind(
    data.frame(
      analysis = "Primary original-script groups",
      samples = ncol(state$dds),
      as.list(strict_count(primary)),
      check.names = FALSE
    ),
    data.frame(
      analysis = "Exclude missing EBV DNA",
      samples = length(known_samples),
      as.list(strict_count(known_result)),
      check.names = FALSE
    ),
    data.frame(
      analysis = "Exclude missing EBV DNA and adjust for stage 4C",
      samples = length(stage_samples),
      as.list(strict_count(stage_result)),
      check.names = FALSE
    )
  )
  save_table(
    summary,
    file.path(config$tables_dir, "differential_expression_sensitivity_summary.csv")
  )

  # Summarize CIBERSORT fit quality by group.
  if (!is.null(state$immune_results)) {
    immune <- as.data.frame(state$immune_results, check.names = FALSE)
    immune$sample <- rownames(immune)
    immune$group <- clinical$group[match(immune$sample, clinical$sample)]
    qc_rows <- do.call(rbind, lapply(c("Low_CP", "High_CP"), function(group) {
      table <- immune[immune$group == group, , drop = FALSE]
      data.frame(
        group = group,
        samples = nrow(table),
        permutation_p_below_0_05 = sum(table$`P-value` < 0.05),
        correlation_at_least_0_30 = sum(table$Correlation >= 0.30),
        rmse_at_most_1 = sum(table$RMSE <= 1),
        combined_correlation_and_rmse = sum(
          table$Correlation >= 0.30 & table$RMSE <= 1
        ),
        check.names = FALSE
      )
    }))
    save_table(
      qc_rows,
      file.path(config$tables_dir, "cibersort_qc_by_group.csv")
    )

    # Compare non-stage-4C samples as a small sensitivity check.
    non_stage4c <- clinical$sample[
      clinical$ebv_dna_record != "NA" & clinical$stage4c == "OtherStage"
    ]
    immune_subset <- immune[match(non_stage4c, immune$sample), , drop = FALSE]
    cell_columns <- names(immune)[seq_len(22)]
    immune_rows <- lapply(cell_columns, function(cell) {
      low <- immune_subset[immune_subset$group == "Low_CP", cell]
      high <- immune_subset[immune_subset$group == "High_CP", cell]
      data.frame(
        cell_type = cell,
        low_samples = length(low),
        high_samples = length(high),
        low_mean = mean(low),
        high_mean = mean(high),
        p_value = stats::wilcox.test(high, low, exact = FALSE)$p.value,
        check.names = FALSE
      )
    })
    immune_sensitivity <- do.call(rbind, immune_rows)
    immune_sensitivity$adjusted_p_value <- stats::p.adjust(
      immune_sensitivity$p_value,
      method = "BH"
    )
    save_table(
      immune_sensitivity,
      file.path(config$tables_dir, "cibersort_non_stage4c_sensitivity.csv")
    )
  }

  invisible(list(clinical = clinical, results = sensitivity, summary = summary))
}
