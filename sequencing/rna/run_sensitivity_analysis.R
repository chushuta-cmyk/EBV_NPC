# Run sensitivity checks with GEO clinical metadata.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_argument[[1]])
project_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

source(file.path(project_dir, "R", "00_setup.R"))
source(file.path(project_dir, "R", "10_sensitivity.R"))

external_library <- Sys.getenv("RNA_EXTERNAL_LIBRARY", unset = "")
if (nzchar(external_library) && dir.exists(external_library)) {
  .libPaths(unique(c(external_library, .libPaths())))
}

config <- make_config(project_dir)
state <- readRDS(file.path(config$cache_dir, "pipeline_state.rds"))
series_matrix_file <- Sys.getenv(
  "RNA_SERIES_MATRIX",
  file.path(config$input_dir, "GSE224450_series_matrix.txt.gz")
)
run_sensitivity_checks(state, config, series_matrix_file)
cat("Sensitivity analysis complete.\n")
