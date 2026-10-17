# Test file structure and small helper functions.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- sub("^--file=", "", file_argument[[1]])
project_dir <- dirname(dirname(normalizePath(test_path, mustWork = TRUE)))

r_files <- c(
  file.path(project_dir, "run_analysis.R"),
  file.path(project_dir, "EBV_resequencing.R"),
  list.files(file.path(project_dir, "R"), pattern = "\\.R$", full.names = TRUE)
)
for (path in r_files) {
  parse(file = path)
  comments <- grep("^[[:space:]]*#", readLines(path, warn = FALSE), value = TRUE)
  stopifnot(!any(grepl("[^ -~]", comments)))
}

source(file.path(project_dir, "R", "00_setup.R"))
stopifnot(
  identical(ratio_to_number(c("1/2", "3/4")), c(0.5, 0.75)),
  is.na(ratio_to_number("bad")),
  all(is.finite(safe_scale_rows(matrix(c(1, 1, 1, 1, 2, 3), nrow = 2, byrow = TRUE))))
)

test_output <- tempfile("rna-test-")
old_output <- Sys.getenv("RNA_OUTPUT_DIR", unset = NA_character_)
Sys.setenv(RNA_OUTPUT_DIR = test_output)
config <- make_config(project_dir)
stopifnot(
  dir.exists(config$output_dir),
  dir.exists(config$tables_dir),
  dir.exists(config$figures_dir),
  dir.exists(config$cache_dir)
)
unlink(test_output, recursive = TRUE)
if (is.na(old_output)) {
  Sys.unsetenv("RNA_OUTPUT_DIR")
} else {
  Sys.setenv(RNA_OUTPUT_DIR = old_output)
}

cat("All structure and helper tests passed.\n")
