# Run selected analysis modules in the correct order.

find_project_dir <- function() {
  configured <- Sys.getenv("RNA_PROJECT_DIR", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))

  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) > 0) {
    script_path <- sub("^--file=", "", file_argument[[1]])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

project_dir <- find_project_dir()
source(file.path(project_dir, "R", "00_setup.R"))

module_files <- c(
  data = "01_data.R",
  figures = "02_figures.R",
  enrichment = "03_enrichment.R",
  immune = "04_immune.R",
  hub = "05_hub_genes.R",
  ml = "06_machine_learning.R",
  wgcna = "07_wgcna.R",
  pathways = "08_pathways.R",
  survival = "09_survival.R"
)
module_functions <- c(
  data = "run_data_module",
  figures = "run_figures_module",
  enrichment = "run_enrichment_module",
  immune = "run_immune_module",
  hub = "run_hub_module",
  ml = "run_ml_module",
  wgcna = "run_wgcna_module",
  pathways = "run_pathways_module",
  survival = "run_survival_module"
)
module_dependencies <- list(
  data = character(),
  figures = "data",
  enrichment = "data",
  immune = "data",
  hub = c("data", "immune"),
  ml = c("data", "immune"),
  wgcna = "data",
  pathways = c("data", "immune"),
  survival = character()
)

resolve_modules <- function(requested) {
  resolved <- character()
  add_module <- function(name) {
    if (!name %in% names(module_files)) {
      stop("Unknown module: ", name, call. = FALSE)
    }
    for (dependency in module_dependencies[[name]]) add_module(dependency)
    resolved <<- unique(c(resolved, name))
  }
  for (name in requested) add_module(name)
  resolved
}

arguments <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% arguments) {
  cat(
    "Usage: Rscript run_analysis.R [--modules=core|all|name,name] [--dry-run]\n",
    "\n",
    "Environment variables:\n",
    "  RNA_COUNTS_FILE                 Local GSE224450 count matrix\n",
    "  RNA_CIBERSORT_FILE              Local CIBERSORT result table\n",
    "  RNA_SURVIVAL_PHENOTYPE_FILE     Local GSE102349 metadata JSON\n",
    "  RNA_SURVIVAL_EXPRESSION_FILE    Local GSE102349 expression matrix\n",
    "  RNA_OUTPUT_DIR                  Output directory\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}
if ("--list-modules" %in% arguments) {
  cat(paste(names(module_files), collapse = "\n"), "\n")
  quit(save = "no", status = 0)
}

module_argument <- grep("^--modules=", arguments, value = TRUE)
selection <- if (length(module_argument) == 0) {
  "core"
} else {
  sub("^--modules=", "", module_argument[[1]])
}
requested <- switch(
  selection,
  core = c("data", "figures"),
  all = names(module_files),
  strsplit(selection, ",", fixed = TRUE)[[1]]
)
selected_modules <- resolve_modules(trimws(requested))

cat("Modules:", paste(selected_modules, collapse = " -> "), "\n")
if ("--dry-run" %in% arguments) {
  quit(save = "no", status = 0)
}

package_requirements <- list(
  data = "DESeq2",
  figures = c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer"),
  enrichment = c("clusterProfiler", "org.Hs.eg.db", "ggplot2"),
  immune = c("ggplot2", "tidyr", "pheatmap", "RColorBrewer"),
  hub = c("ggplot2", "pheatmap", "randomForest"),
  ml = c("caret", "pROC", "randomForest", "xgboost", "glmnet", "kernlab", "ggplot2"),
  wgcna = "WGCNA",
  pathways = c("ggplot2", "pheatmap", "org.Hs.eg.db", "AnnotationDbi"),
  survival = "survival"
)
packages <- unique(unlist(package_requirements[selected_modules], use.names = FALSE))
if ("data" %in% selected_modules && !nzchar(Sys.getenv("RNA_COUNTS_FILE", unset = ""))) {
  packages <- c(packages, "GEOquery")
}
if ("immune" %in% selected_modules && !nzchar(Sys.getenv("RNA_CIBERSORT_FILE", unset = ""))) {
  packages <- c(packages, "CIBERSORT")
}
if ("survival" %in% selected_modules) {
  if (nzchar(Sys.getenv("RNA_SURVIVAL_PHENOTYPE_FILE", unset = ""))) {
    packages <- c(packages, "jsonlite")
  } else {
    packages <- c(packages, "GEOquery")
  }
}
required_packages(unique(packages))

config <- make_config(project_dir)
state <- new.env(parent = emptyenv())
set.seed(config$random_seed)

for (module in selected_modules) {
  source(file.path(project_dir, "R", module_files[[module]]))
  cat("Running", module, "module...\n")
  runner <- get(module_functions[[module]], mode = "function")
  tryCatch(
    runner(state, config),
    error = function(error) {
      stop("Module '", module, "' failed: ", conditionMessage(error), call. = FALSE)
    }
  )
}

saveRDS(state, file.path(config$cache_dir, "pipeline_state.rds"))
cat("Pipeline complete. Output:", config$output_dir, "\n")
