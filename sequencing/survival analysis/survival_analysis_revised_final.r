# Use bind_rows to avoid column mismatch and keep one row per gene

surv_results <- tibble::tibble(
  gene = character(),
  n = integer(),
  events = integer(),
  cutoff_median = numeric(),
  coef = numeric(),
  HR = numeric(),
  lower95 = numeric(),
  upper95 = numeric(),
  p_value = numeric(),
  ph_p_value = numeric()
)

for (gene in matched_genes) {
  if (!gene %in% colnames(merged)) next

  gene_df <- merged %>%
    dplyr::select(geo_accession, time, event, all_of(gene)) %>%
    filter(!is.na(.data[[gene]]))

  if (nrow(gene_df) < 5) next
  if (sd(gene_df[[gene]], na.rm = TRUE) == 0) next
  if (dplyr::n_distinct(gene_df[[gene]]) < 2) next

  cutoff <- median(gene_df[[gene]], na.rm = TRUE)
  gene_df$group <- factor(ifelse(gene_df[[gene]] > cutoff, "High", "Low"), levels = c("Low", "High"))
  surv_obj <- Surv(gene_df$time, gene_df$event)

  gene_q <- paste0("`", gene, "`")
  cox_model <- coxph(as.formula(paste0("surv_obj ~ ", gene_q)), data = gene_df)
  cox_sum <- summary(cox_model)
  ph_p <- as.numeric(cox.zph(cox_model)$table[1, "p"])

  row <- tibble::tibble(
    gene = gene,
    n = nrow(gene_df),
    events = sum(gene_df$event, na.rm = TRUE),
    cutoff_median = cutoff,
    coef = as.numeric(cox_sum$coefficients[1, "coef"]),
    HR = as.numeric(cox_sum$coefficients[1, "exp(coef)"]),
    lower95 = as.numeric(cox_sum$conf.int[1, "lower .95"]),
    upper95 = as.numeric(cox_sum$conf.int[1, "upper .95"]),
    p_value = as.numeric(cox_sum$coefficients[1, "Pr(>|z|)"]),
    ph_p_value = ph_p
  )
  surv_results <- dplyr::bind_rows(surv_results, row)

  fit <- survfit(surv_obj ~ group, data = gene_df)
  p <- ggsurvplot(
    fit,
    data = gene_df,
    conf.int = TRUE,
    conf.int.alpha = 0.18,
    pval = TRUE,
    pval.method = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.y.text = FALSE,
    risk.table.y.text.col = TRUE,
    break.time.by = 10,
    censor.shape = 124,
    censor.size = 3.2,
    linetype = 1,
    size = 1.0,
    palette = c("#4C78A8", "#D64B4B"),
    title = gene,
    legend.title = NULL,
    legend.labs = c("Low", "High"),
    legend = c(0.80, 0.92),
    xlab = "Time to Event (Months)",
    ylab = "Progression-free survival probability",
    font.main = c(13, "bold"),
    font.x = c(12),
    font.y = c(12),
    font.legend = c(10),
    font.tickslab = c(10),
    risk.table.height = 0.20,
    risk.table.fontsize = 3.8,
    ggtheme = theme_classic(base_size = 12),
    tables.theme = theme_cleantable()
  )

  p$plot <- p$plot + theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.80, 0.92),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

  p$table <- p$table + theme(
    axis.title.x = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

  ggsave(file.path(graphs_dir, paste0("KM_", gene, ".pdf")), plot = arrange_ggsurvplots(list(p), print = FALSE), width = 7.2, height = 6.0, device = cairo_pdf)
  ggsave(file.path(graphs_dir, paste0("KM_", gene, ".svg")), plot = arrange_ggsurvplots(list(p), print = FALSE), width = 7.2, height = 6.0)
  ggsave(file.path(graphs_dir, paste0("KM_", gene, ".png")), plot = arrange_ggsurvplots(list(p), print = FALSE), width = 7.2, height = 6.0, dpi = 600)
}

surv_results$padj_BH <- p.adjust(surv_results$p_value, method = "BH")
write.csv(surv_results, file.path(tables_dir, "cox_results_revised.csv"), row.names = FALSE)
print(surv_results)
