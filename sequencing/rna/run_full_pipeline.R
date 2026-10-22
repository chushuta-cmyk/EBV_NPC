# Run the complete raw RCC analysis.

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", file_argument[[1]])
project_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

input_dir <- Sys.getenv(
  "RNA_INPUT_DIR",
  file.path(project_dir, "data", "inputs")
)
if (!dir.exists(input_dir)) {
  stop(
    "RNA_INPUT_DIR is not set and the default '", input_dir,
    "' does not exist. Set RNA_INPUT_DIR to the external input-data ",
    "directory (see README).",
    call. = FALSE
  )
}
Sys.setenv(
  RNA_PROJECT_DIR = project_dir,
  RNA_MODULES = "all",
  RNA_INPUT_DIR = input_dir,
  RNA_CIBERSORT_PERM = "1000",
  RNA_SURVIVAL_PHENOTYPE_FILE = file.path(
    input_dir,
    "GSE102349_samples.json"
  ),
  RNA_SURVIVAL_EXPRESSION_FILE = file.path(
    input_dir,
    "GSE102349_NPC_mRNA_processed.txt.gz"
  )
)

source(file.path(project_dir, "run_analysis.R"))
