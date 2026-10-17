# Keep the old file name as a pipeline entry point.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_argument[[1]])
project_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
Sys.setenv(RNA_PROJECT_DIR = project_dir)
source(file.path(project_dir, "run_analysis.R"))
