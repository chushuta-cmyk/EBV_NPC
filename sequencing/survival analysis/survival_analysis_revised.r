################################################################################
# Survival Analysis with Cox and KM Curves (Publication-Ready Revised)
################################################################################

library(GEOquery)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# ---------------------------
# 1. Load data
# ---------------------------
geo <- getGEO("GSE102349", GSEMatrix = TRUE)
pheno_data <- pData(geo[[1]])
expr_data <- read.delim("GSE102349_NPC_mRNA_processed.txt", sep = "\t", row.names = 1, check.names = FALSE)

if (ncol(expr_data) != nrow(pheno_data)) {
  stop("Column count of expression matrix does not match phenotype rows.")
}
colnames(expr_data) <- rownames(pheno_data)

# ---------------------------
# 2. Define candidate genes
# ---------------------------
up_intersect <- c("CD3E", "CD247", "LCK", "IL2RG", "HLA-DPB1", "MS4A1")
down_intersect <- c("CXCL8", "IL1B", "C5AR1", "TREM1", "CXCR2", "CSF3R", "S100A9", "S100A8", "NLRP3", "IL1A", "THBD", "ITGAX")
genes_for_survival_analysis <- unique(c(up_intersect, down_intersect))

matched_genes <- genes_for_survival_analysis[genes_for_survival_analysis %in% rownames(expr_data)]
if (length(matched_genes) == 0) {
  stop("None of the candidate genes were found in the expression matrix.")
}

expr_df <- as.data.frame(t(expr_data[matched_genes, , drop = FALSE]))
expr_df$geo_accession <- rownames(expr_df)

# ---------------------------
# 3. Clean survival data
# ---------------------------
pheno_data$`time to event:ch1`[pheno_data$`time to event:ch1` %in% c("N/A", "NA", "")] <- NA
pheno_data$`time to event:ch1` <- as.numeric(pheno_data$`time to event:ch1`)

surv_df <- pheno_data %>%
  select(geo_accession, `time to event:ch1`, `event:ch1`) %>%
  rename(time = `time to event:ch1`, event = `event:ch1`) %>%
  mutate(event = ifelse(event == "Disease progression", 1, 0)) %>%
  filter(!is.na(time), !is.na(event))

merged <- merge(surv_df, expr_df, by = "geo_accession")

if (nrow(merged) == 0) {
  stop("No samples remained after merging expression and survival data.")
}

# ---------------------------
# 4. Output directories
# ---------------------------
graphs_dir <- "KM_plots"
tables_dir <- "results"
dir.create(graphs_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

# ---------------------------
# 5. Univariate Cox + KM
# ---------------------------
surv_results <- data.frame(
  gene = character(),
  n = integer(),
  events = integer(),
  cutoff_median = numeric(),
  coef = numeric(),
  HR = numeric(),
  lower95 = numeric(),
  upper95 = numeric(),
  p_value = numeric(),
  ph_p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in matched_genes) {
  if (!gene %in% colnames(merged)) next

  gene_df <- merged %>%
    select(geo_accession, time, event, all_of(gene)) %>%
    filter(!is.na(.data[[gene]]))

  if (nrow(gene_df) < 5) next
  if (sd(gene_df[[gene]], na.rm = TRUE) == 0) next

  cutoff <- median(gene_df[[gene]], na.rm = TRUE)
  gene_df$group <- ifelse(gene_df[[gene]] > cutoff, "High", "Low")
  gene_df$group <- factor(gene_df$group, levels = c("Low", "High"))

  surv_obj <- Surv(gene_df$time, gene_df$event)

  gene_q <- paste0("`", gene, "`") 
  cox_formula <- as.formula(paste0("surv_obj ~ ", gene_q))
  cox_model <- coxph(cox_formula, data = gene_df)
  cox_sum <- summary(cox_model)

  zph <- cox.zph(cox_model)
  ph_p <- as.numeric(zph$table[1, "p"])

  coef_val <- as.numeric(cox_sum$coefficients[1, "coef"])
  hr_val <- as.numeric(cox_sum$coefficients[1, "exp(coef)"])
  pval <- as.numeric(cox_sum$coefficients[1, "Pr(>|z|)"])
  lower95 <- as.numeric(cox_sum$conf.int[1, "lower .95"])
  upper95 <- as.numeric(cox_sum$conf.int[1, "upper .95"])

  surv_results <- rbind(
    surv_results,
    data.frame(
      gene = gene,
      n = nrow(gene_df),
      events = sum(gene_df$event, na.rm = TRUE),
      cutoff_median = cutoff,
      coef = coef_val,
      HR = hr_val,
      lower95 = lower95,
      upper95 = upper95,
      p_value = pval,
      ph_p_value = ph_p
    )
  )

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
  
  p$plot <- p$plot +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.80, 0.92),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5)
    )
  
  p$table <- p$table +
    theme(
      axis.title.x = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      plot.title = element_text(size = 11, face = "bold"),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = file.path(graphs_dir, paste0("KM_", gene, ".pdf")),
    plot = arrange_ggsurvplots(list(p), print = FALSE),
    width = 7.2,
    height = 6.0,
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(graphs_dir, paste0("KM_", gene, ".svg")),
    plot = arrange_ggsurvplots(list(p), print = FALSE),
    width = 7.2,
    height = 6.0
  )
  
  ggsave(
    filename = file.path(graphs_dir, paste0("KM_", gene, ".png")),
    plot = arrange_ggsurvplots(list(p), print = FALSE),
    width = 7.2,
    height = 6.0,
    dpi = 600
  )
}

# ---------------------------
# 6. Multiple-testing correction and export
# ---------------------------
surv_results$padj_BH <- p.adjust(surv_results$p_value, method = "BH")
surv_results <- surv_results %>% arrange(p_value)

write.csv(surv_results, file.path(tables_dir, "cox_results_revised.csv"), row.names = FALSE)
write.csv(subset(surv_results, p_value < 0.05), file.path(tables_dir, "cox_results_nominally_significant.csv"), row.names = FALSE)

print(surv_results)
