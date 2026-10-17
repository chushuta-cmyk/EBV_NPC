# Compare models that predict the sample group.

run_ml_module <- function(state, config) {
  required_packages(c(
    "caret", "pROC", "randomForest", "xgboost", "glmnet", "kernlab", "ggplot2"
  ))
  need_state(state, c("immune_cells", "sample_info"))
  set.seed(config$random_seed)

  predictors <- state$immune_cells[rownames(state$sample_info), , drop = FALSE]
  outcome <- factor(state$sample_info$Group, levels = c("Low_CP", "High_CP"))
  correlation <- stats::cor(predictors, use = "pairwise.complete.obs")
  remove <- caret::findCorrelation(correlation, cutoff = 0.8, names = TRUE)
  predictors <- predictors[, !colnames(predictors) %in% remove, drop = FALSE]
  if (ncol(predictors) < 2) {
    stop("Too few immune features remain after correlation filtering.", call. = FALSE)
  }

  control <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )

  train_model <- function(method, grid = NULL, tune_length = NULL, ...) {
    arguments <- list(
      x = predictors,
      y = outcome,
      method = method,
      metric = "ROC",
      trControl = control,
      preProcess = c("center", "scale")
    )
    if (!is.null(grid)) arguments$tuneGrid <- grid
    if (!is.null(tune_length)) arguments$tuneLength <- tune_length
    arguments <- c(arguments, list(...))
    do.call(caret::train, arguments)
  }

  models <- list(
    `Random Forest` = train_model(
      "rf",
      grid = expand.grid(mtry = unique(pmin(c(2, 3, 5), ncol(predictors))))
    ),
    XGBoost = train_model("xgbTree", tune_length = 3),
    `Elastic Net` = train_model(
      "glmnet",
      grid = expand.grid(
        alpha = seq(0, 1, 0.25),
        lambda = 10^seq(-3, 0, length.out = 12)
      )
    ),
    Lasso = train_model(
      "glmnet",
      grid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length.out = 12))
    ),
    SVM = train_model("svmRadial", tune_length = 5),
    KNN = train_model("knn", tune_length = 8),
    Logistic = train_model("glm", family = stats::binomial())
  )

  performance <- data.frame(
    model = names(models),
    ROC = vapply(models, function(model) max(model$results$ROC, na.rm = TRUE), numeric(1)),
    row.names = NULL
  )
  performance <- performance[order(performance$ROC, decreasing = TRUE), ]

  performance_plot <- ggplot2::ggplot(
    performance,
    ggplot2::aes(reorder(model, ROC), ROC, fill = model)
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", ROC)), vjust = -0.3) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Cross-validated model ROC", x = NULL, y = "ROC") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  ggplot2::ggsave(
    file.path(config$figures_dir, "model_roc_comparison.png"),
    performance_plot,
    width = 10,
    height = 7,
    dpi = 300
  )

  importance <- lapply(names(models), function(name) {
    value <- caret::varImp(models[[name]], scale = TRUE)$importance
    data.frame(
      feature = rownames(value),
      importance = apply(value, 1, max, na.rm = TRUE),
      model = name,
      row.names = NULL
    )
  })
  importance <- do.call(rbind, importance)

  save_table(performance, file.path(config$tables_dir, "model_performance.csv"))
  save_table(importance, file.path(config$tables_dir, "model_feature_importance.csv"))
  saveRDS(models, file.path(config$cache_dir, "trained_models.rds"))

  state$models <- models
  state$model_performance <- performance
  state$removed_immune_features <- remove
  message("Machine-learning module complete.")
  invisible(state)
}
