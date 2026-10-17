# Run Cox models and Kaplan-Meier plots.

run_survival_module <- function(state, config) {
  required_packages("survival")
  if (!nzchar(config$survival_expression_file)) {
    stop(
      "Set RNA_SURVIVAL_EXPRESSION_FILE to the GSE102349 expression matrix.",
      call. = FALSE
    )
  }
  if (!file.exists(config$survival_expression_file)) {
    stop("The survival expression file does not exist.", call. = FALSE)
  }

  if (nzchar(config$survival_phenotype_file)) {
    required_packages("jsonlite")
    phenotype_json <- jsonlite::fromJSON(config$survival_phenotype_file)
    phenotype <- data.frame(
      geo_accession = phenotype_json$geo_accession,
      expression_id = phenotype_json$description,
      time = phenotype_json[["time to event"]],
      event_text = phenotype_json$event,
      stringsAsFactors = FALSE
    )
  } else {
    required_packages("GEOquery")
    geo <- GEOquery::getGEO(
      "GSE102349",
      GSEMatrix = TRUE,
      destdir = config$cache_dir
    )
    geo_phenotype <- Biobase::pData(geo[[1]])
    description_column <- grep("^description", colnames(geo_phenotype), value = TRUE)[1]
    phenotype <- data.frame(
      geo_accession = geo_phenotype$geo_accession,
      expression_id = geo_phenotype[[description_column]],
      time = geo_phenotype[["time to event:ch1"]],
      event_text = geo_phenotype[["event:ch1"]],
      stringsAsFactors = FALSE
    )
  }
  expression <- utils::read.delim(
    config$survival_expression_file,
    row.names = 1,
    check.names = FALSE
  )

  hub_genes <- c(
    "CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1", "CXCL8", "IL1B",
    "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3",
    "IL1A", "THBD", "ITGAX"
  )
  genes <- intersect(hub_genes, rownames(expression))
  if (length(genes) == 0) {
    stop("No survival genes were found in the expression matrix.", call. = FALSE)
  }

  survival_data <- data.frame(
    sample = phenotype$expression_id,
    geo_accession = phenotype$geo_accession,
    time = suppressWarnings(as.numeric(phenotype$time)),
    event = as.numeric(phenotype$event_text == "Disease progression"),
    stringsAsFactors = FALSE
  )
  expression_data <- data.frame(
    sample = colnames(expression),
    t(expression[genes, , drop = FALSE]),
    check.names = FALSE,
    row.names = NULL
  )
  merged <- merge(survival_data, expression_data, by = "sample")
  merged <- merged[is.finite(merged$time), , drop = FALSE]

  result_rows <- lapply(genes, function(gene) {
    model <- survival::coxph(
      survival::Surv(time, event) ~ merged[[gene]],
      data = merged
    )
    model_summary <- summary(model)
    group <- ifelse(merged[[gene]] > stats::median(merged[[gene]], na.rm = TRUE), "High", "Low")
    fit <- survival::survfit(survival::Surv(time, event) ~ group, data = merged)

    grDevices::png(
      file.path(config$figures_dir, paste0("survival_", gene, ".png")),
      width = 1400,
      height = 1000,
      res = 160
    )
    graphics::plot(
      fit,
      col = c("#377EB8", "#E41A1C"),
      lwd = 2,
      xlab = "Time to event",
      ylab = "Survival probability",
      main = paste("Survival by", gene, "expression")
    )
    graphics::legend(
      "bottomleft",
      legend = names(fit$strata),
      col = c("#377EB8", "#E41A1C"),
      lwd = 2,
      bty = "n"
    )
    grDevices::dev.off()

    data.frame(
      gene = gene,
      coefficient = model_summary$coefficients[1, "coef"],
      hazard_ratio = model_summary$coefficients[1, "exp(coef)"],
      p_value = model_summary$coefficients[1, "Pr(>|z|)"],
      row.names = NULL
    )
  })

  survival_results <- do.call(rbind, result_rows)
  save_table(survival_results, file.path(config$tables_dir, "survival_cox_results.csv"))
  state$survival_results <- survival_results
  message("Survival module complete.")
  invisible(state)
}
