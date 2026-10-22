# Run the updated STRING network analysis.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_argument[[1]])
project_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

source(file.path(project_dir, "R", "00_setup.R"))
source(file.path(project_dir, "R", "11_ppi.R"))

external_library <- Sys.getenv("RNA_EXTERNAL_LIBRARY", unset = "")
if (nzchar(external_library) && dir.exists(external_library)) {
  .libPaths(unique(c(external_library, .libPaths())))
}

config <- make_config(project_dir)
state <- readRDS(file.path(config$cache_dir, "pipeline_state.rds"))
run_ppi_analysis(state, config)
cat("PPI analysis complete.\n")
