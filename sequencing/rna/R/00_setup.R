# Shared settings and small helper functions.

required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install them before running this module.",
      call. = FALSE
    )
  }
}

need_state <- function(state, names) {
  missing <- names[!vapply(names, exists, logical(1), envir = state, inherits = FALSE)]
  if (length(missing) > 0) {
    stop("Missing pipeline data: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

ratio_to_number <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(value) {
    if (length(value) != 2) return(NA_real_)
    as.numeric(value[[1]]) / as.numeric(value[[2]])
  }, numeric(1))
}

safe_scale_rows <- function(x) {
  scaled <- t(scale(t(as.matrix(x))))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

save_table <- function(x, path, row.names = FALSE) {
  utils::write.csv(as.data.frame(x), path, row.names = row.names)
  invisible(path)
}

get_env_path <- function(name) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return("")
  normalizePath(value, mustWork = FALSE)
}

make_config <- function(project_dir) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  output_dir <- Sys.getenv("RNA_OUTPUT_DIR", file.path(project_dir, "output"))

  list(
    project_dir = project_dir,
    output_dir = make_dir(output_dir),
    tables_dir = make_dir(file.path(output_dir, "tables")),
    figures_dir = make_dir(file.path(output_dir, "figures")),
    cache_dir = make_dir(file.path(output_dir, "cache")),
    counts_file = get_env_path("RNA_COUNTS_FILE"),
    cibersort_file = get_env_path("RNA_CIBERSORT_FILE"),
    survival_phenotype_file = get_env_path("RNA_SURVIVAL_PHENOTYPE_FILE"),
    survival_expression_file = get_env_path("RNA_SURVIVAL_EXPRESSION_FILE"),
    random_seed = 123L
  )
}
