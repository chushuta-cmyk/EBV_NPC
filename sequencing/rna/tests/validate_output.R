# Validate the complete pipeline output.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- sub("^--file=", "", file_argument[[1]])
project_dir <- dirname(dirname(normalizePath(test_path, mustWork = TRUE)))
output_dir <- Sys.getenv("RNA_OUTPUT_DIR", file.path(project_dir, "output"))
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

read_table <- function(name) {
  utils::read.csv(file.path(tables_dir, name), check.names = FALSE)
}

de <- read_table("differential_expression.csv")
status <- read_table("module_run_status.csv")
go_up <- read_table("go_upregulated.csv")
go_down <- read_table("go_downregulated.csv")
kegg_up <- read_table("kegg_upregulated.csv")
kegg_down <- read_table("kegg_downregulated.csv")
ml <- read_table("model_performance.csv")
wgcna <- read_table("wgcna_module_trait.csv")
survival <- read_table("survival_cox_results.csv")
probe_audit <- read_table("raw_rcc_probe_audit.csv")
immune_audit <- read_table("immune_method_audit.csv")

strict <- !is.na(de$padj) & de$padj < 0.01 & abs(de$log2FoldChange) > 2
png_files <- list.files(figures_dir, pattern = "\\.png$", full.names = TRUE)

checks <- data.frame(
  check = c(
    "Nine modules completed",
    "DE table has 730 probes",
    "Strict DE count is 78",
    "Strict DE directions are 11 up and 67 down",
    "Raw RCC values are integers",
    "No rounding was applied",
    "At least four immune methods completed",
    "GO tables are not empty",
    "KEGG uses the targeted panel background",
    "ML ROC values are finite",
    "WGCNA table is not empty",
    "Survival table is not empty",
    "All PNG files are not empty"
  ),
  passed = c(
    nrow(status) == 9 && all(status$status == "complete"),
    nrow(de) == 730,
    sum(strict) == 78,
    sum(strict & de$log2FoldChange > 0) == 11 &&
      sum(strict & de$log2FoldChange < 0) == 67,
    probe_audit$value[probe_audit$metric == "Raw RCC values are integers"] == "TRUE",
    probe_audit$value[probe_audit$metric == "Rounding applied"] == "FALSE",
    sum(immune_audit$status == "completed") >= 4,
    nrow(go_up) > 0 && nrow(go_down) > 0,
    nrow(kegg_down) > 0 &&
      all(as.numeric(sub(".*/", "", kegg_down$BgRatio)) <= 730),
    all(is.finite(ml$ROC)),
    nrow(wgcna) > 0,
    nrow(survival) > 0,
    length(png_files) > 0 && all(file.info(png_files)$size > 0)
  ),
  check.names = FALSE
)

summary <- data.frame(
  metric = c(
    "Low_CP samples", "High_CP samples", "Filtered endogenous probes",
    "Strict significant probes", "Strict upregulated probes",
    "Strict downregulated probes", "GO up terms", "GO down terms",
    "KEGG up pathways", "KEGG down pathways", "WGCNA modules",
    "Survival genes", "Machine-learning models completed",
    "Immune methods completed"
  ),
  value = c(
    5, 17, nrow(de), sum(strict),
    sum(strict & de$log2FoldChange > 0),
    sum(strict & de$log2FoldChange < 0),
    nrow(go_up), nrow(go_down), nrow(kegg_up), nrow(kegg_down),
    nrow(wgcna), nrow(survival), nrow(ml),
    sum(immune_audit$status == "completed")
  ),
  check.names = FALSE
)

utils::write.csv(
  checks,
  file.path(tables_dir, "pipeline_validation_checks.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary,
  file.path(tables_dir, "pipeline_summary.csv"),
  row.names = FALSE
)

print(checks, row.names = FALSE)
if (!all(checks$passed)) {
  stop("One or more pipeline checks failed.", call. = FALSE)
}
cat("All pipeline output checks passed.\n")
