source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

xgboost_has_config_key <- function(config, name) {
  !is.null(config) && !is.null(names(config)) && name %in% names(config)
}

is_missing_xgboost_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

xgboost_scalar <- function(value, default = NULL) {
  if (is_missing_xgboost_value(value)) {
    return(default)
  }
  as.character(value[[1]])
}

xgboost_numeric <- function(value, default, field, lower = NULL, upper = NULL) {
  if (is_missing_xgboost_value(value)) {
    return(default)
  }
  numeric_value <- suppressWarnings(as.numeric(value[[1]]))
  if (is.na(numeric_value) || !is.finite(numeric_value)) {
    stop(field, " must be numeric.", call. = FALSE)
  }
  if (!is.null(lower) && numeric_value <= lower) {
    stop(field, " must be greater than ", lower, ".", call. = FALSE)
  }
  if (!is.null(upper) && numeric_value >= upper) {
    stop(field, " must be less than ", upper, ".", call. = FALSE)
  }
  numeric_value
}

xgboost_integer <- function(value, default, field, lower = NULL) {
  numeric_value <- xgboost_numeric(value, default, field, lower = lower)
  integer_value <- as.integer(numeric_value)
  if (!identical(as.numeric(integer_value), as.numeric(numeric_value))) {
    stop(field, " must be an integer.", call. = FALSE)
  }
  integer_value
}

xgboost_character_vector <- function(value, field) {
  values <- as.character(unlist(value, use.names = FALSE))
  values <- values[!is.na(values) & nzchar(trimws(values))]
  values <- unique(trimws(values))
  if (length(values) == 0) {
    stop(field, " must contain at least one value.", call. = FALSE)
  }
  values
}

sanitize_xgboost_path_component <- function(x, fallback = "model") {
  value <- xgboost_scalar(x, fallback)
  value <- trimws(value)
  value <- gsub("[[:space:]/\\\\:;,*?\"<>|]+", "_", value)
  value <- gsub("[^A-Za-z0-9._+=@-]", "_", value)
  value <- gsub("_+", "_", value)
  value <- sub("^_+", "", value)
  value <- sub("_+$", "", value)
  if (!nzchar(value)) fallback else value
}

normalize_xgboost_params <- function(global_params = list(), model_params = list()) {
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 4,
    eta = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    gamma = 0
  )
  params <- utils::modifyList(params, global_params %||% list(), keep.null = FALSE)
  params <- utils::modifyList(params, model_params %||% list(), keep.null = FALSE)
  params$objective <- "binary:logistic"
  params$eval_metric <- "auc"
  params
}

normalize_xgboost_target <- function(target, model_name) {
  if (is.null(target) || length(target) == 0) {
    stop("XGBoost model '", model_name, "' must define target.", call. = FALSE)
  }
  column <- xgboost_scalar(target$column)
  if (is.null(column)) {
    stop("XGBoost model '", model_name, "' must define target.column.", call. = FALSE)
  }
  negative <- xgboost_character_vector(
    target$negative,
    paste0("XGBoost model '", model_name, "' target.negative")
  )
  positive <- xgboost_character_vector(
    target$positive,
    paste0("XGBoost model '", model_name, "' target.positive")
  )
  overlap <- intersect(negative, positive)
  if (length(overlap) > 0) {
    stop(
      "XGBoost model '",
      model_name,
      "' target.negative and target.positive overlap: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }
  list(column = column, negative = negative, positive = positive)
}

normalize_xgboost_model <- function(spec, global_config, project_config) {
  spec <- spec %||% list()
  model_name <- xgboost_scalar(spec$name)
  if (is.null(model_name)) {
    stop("Each xgboost_classification.models entry must define a non-empty name.", call. = FALSE)
  }

  split_group_col <- if (xgboost_has_config_key(spec, "split_group_col")) {
    xgboost_scalar(spec$split_group_col, default = NULL)
  } else {
    global_config$split_group_col
  }

  tax_agg_level <- xgboost_scalar(
    spec$tax_agg_level,
    default = global_config$tax_agg_level %||% project_config$tax_agg_level %||% "Genus"
  )
  norm_method <- xgboost_scalar(
    spec$norm_method,
    default = global_config$norm_method %||% project_config$norm_method %||% "log2HostMapped"
  )
  pseudocount <- xgboost_numeric(
    spec$pseudocount,
    default = global_config$pseudocount %||% project_config$pseudocount %||% 1,
    field = paste0("XGBoost model '", model_name, "' pseudocount")
  )

  list(
    name = model_name,
    safe_name = sanitize_xgboost_path_component(model_name),
    plot_title = xgboost_scalar(spec$plot_title, default = model_name),
    sample_filter = spec$sample_filter,
    target = normalize_xgboost_target(spec$target, model_name),
    tax_agg_level = tax_agg_level,
    norm_method = norm_method,
    pseudocount = pseudocount,
    batch_adj_covar = if (xgboost_has_config_key(spec, "batch_adj_covar")) spec$batch_adj_covar else global_config$batch_adj_covar,
    batch_adj_formula = if (xgboost_has_config_key(spec, "batch_adj_formula")) spec$batch_adj_formula else global_config$batch_adj_formula,
    batch_adj_method = xgboost_scalar(spec$batch_adj_method, default = global_config$batch_adj_method),
    split_group_col = split_group_col,
    test_fraction = xgboost_numeric(
      spec$test_fraction,
      default = global_config$test_fraction,
      field = paste0("XGBoost model '", model_name, "' test_fraction"),
      lower = 0,
      upper = 1
    ),
    cv_folds = xgboost_integer(
      spec$cv_folds,
      default = global_config$cv_folds,
      field = paste0("XGBoost model '", model_name, "' cv_folds"),
      lower = 1
    ),
    seed = xgboost_integer(
      spec$seed,
      default = global_config$seed,
      field = paste0("XGBoost model '", model_name, "' seed")
    ),
    nrounds = xgboost_integer(
      spec$nrounds,
      default = global_config$nrounds,
      field = paste0("XGBoost model '", model_name, "' nrounds"),
      lower = 0
    ),
    early_stopping_rounds = xgboost_integer(
      spec$early_stopping_rounds,
      default = global_config$early_stopping_rounds,
      field = paste0("XGBoost model '", model_name, "' early_stopping_rounds"),
      lower = -1
    ),
    prediction_threshold = xgboost_numeric(
      spec$prediction_threshold,
      default = global_config$prediction_threshold,
      field = paste0("XGBoost model '", model_name, "' prediction_threshold"),
      lower = 0,
      upper = 1
    ),
    split_max_attempts = xgboost_integer(
      spec$split_max_attempts,
      default = global_config$split_max_attempts,
      field = paste0("XGBoost model '", model_name, "' split_max_attempts"),
      lower = 0
    ),
    xgb_params = normalize_xgboost_params(global_config$xgb_params, spec$xgb_params)
  )
}

normalize_xgboost_config <- function(config = list(), project_config = list()) {
  config <- config %||% list()

  global_config <- list(
    tax_agg_level = xgboost_scalar(config$tax_agg_level, default = project_config$tax_agg_level %||% "Genus"),
    norm_method = xgboost_scalar(config$norm_method, default = project_config$norm_method %||% "log2HostMapped"),
    pseudocount = xgboost_numeric(config$pseudocount, default = project_config$pseudocount %||% 1, field = "xgboost_classification.pseudocount"),
    batch_adj_covar = config$batch_adj_covar,
    batch_adj_formula = config$batch_adj_formula,
    batch_adj_method = xgboost_scalar(config$batch_adj_method, default = "LimmaRemoveBatchEffect"),
    split_group_col = if (xgboost_has_config_key(config, "split_group_col")) {
      xgboost_scalar(config$split_group_col, default = NULL)
    } else {
      "PatientID"
    },
    test_fraction = xgboost_numeric(config$test_fraction, default = 0.2, field = "xgboost_classification.test_fraction", lower = 0, upper = 1),
    cv_folds = xgboost_integer(config$cv_folds, default = 5L, field = "xgboost_classification.cv_folds", lower = 1),
    seed = xgboost_integer(config$seed, default = 42L, field = "xgboost_classification.seed"),
    nrounds = xgboost_integer(config$nrounds, default = 500L, field = "xgboost_classification.nrounds", lower = 0),
    early_stopping_rounds = xgboost_integer(
      config$early_stopping_rounds,
      default = 20L,
      field = "xgboost_classification.early_stopping_rounds",
      lower = -1
    ),
    prediction_threshold = xgboost_numeric(
      config$prediction_threshold,
      default = 0.5,
      field = "xgboost_classification.prediction_threshold",
      lower = 0,
      upper = 1
    ),
    split_max_attempts = xgboost_integer(config$split_max_attempts, default = 100L, field = "xgboost_classification.split_max_attempts", lower = 0),
    xgb_params = config$xgb_params %||% list()
  )

  if (is.null(config$models) || length(config$models) == 0) {
    stop("xgboost_classification.models must define at least one model.", call. = FALSE)
  }

  models <- lapply(config$models, normalize_xgboost_model, global_config = global_config, project_config = project_config)
  model_names <- vapply(models, function(model) model$name, character(1))
  duplicate_names <- unique(model_names[duplicated(model_names)])
  if (length(duplicate_names) > 0) {
    stop(
      "xgboost_classification.models contains duplicate model name(s): ",
      paste(duplicate_names, collapse = ", "),
      call. = FALSE
    )
  }

  global_config$models <- models
  global_config
}

warn_xgboost_target_leakage <- function(spec) {
  if (is_missing_xgboost_value(spec$batch_adj_formula)) {
    return(invisible(FALSE))
  }

  formula_text <- as.character(spec$batch_adj_formula[[1]])
  formula_vars <- tryCatch(
    all.vars(stats::as.formula(formula_text)),
    error = function(e) {
      warning(
        "Could not parse xgboost_classification model '",
        spec$name,
        "' batch_adj_formula for target-leakage screening: ",
        conditionMessage(e),
        call. = FALSE
      )
      character()
    }
  )
  risky_vars <- intersect(formula_vars, spec$target$column)
  if (length(risky_vars) > 0) {
    warning(
      "xgboost_classification model '",
      spec$name,
      "' batch_adj_formula references target variable(s): ",
      paste(risky_vars, collapse = ", "),
      ". This can inflate classifier performance because labels influence preprocessing.",
      call. = FALSE
    )
  }
  invisible(length(risky_vars) > 0)
}

prepare_xgboost_model_metadata <- function(physeq, spec) {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  metadata_df$.sample_name <- rownames(metadata_df)
  target_col <- spec$target$column
  fail_missing_columns(names(metadata_df), target_col, "phyloseq sample_data")

  target_values <- standardize_metadata_missing(metadata_df[[target_col]])
  target_values <- as.character(target_values)
  negative <- as.character(spec$target$negative)
  positive <- as.character(spec$target$positive)
  negative_idx <- !is.na(target_values) & target_values %in% negative
  positive_idx <- !is.na(target_values) & target_values %in% positive
  keep_idx <- negative_idx | positive_idx

  dropped <- sum(!keep_idx)
  if (dropped > 0) {
    message(
      "XGBoost model '",
      spec$name,
      "' dropped ",
      dropped,
      " sample(s) outside configured target values for ",
      target_col,
      "."
    )
  }

  metadata_df <- metadata_df[keep_idx, , drop = FALSE]
  metadata_df$.xgb_class <- ifelse(positive_idx[keep_idx], 1L, 0L)

  if (!is.null(spec$split_group_col)) {
    fail_missing_columns(names(metadata_df), spec$split_group_col, "phyloseq sample_data")
    group_values <- standardize_metadata_missing(metadata_df[[spec$split_group_col]])
    group_keep <- !is_metadata_missing_like(group_values)
    if (any(!group_keep)) {
      message(
        "XGBoost model '",
        spec$name,
        "' dropped ",
        sum(!group_keep),
        " sample(s) with missing split group column ",
        spec$split_group_col,
        "."
      )
    }
    metadata_df <- metadata_df[group_keep, , drop = FALSE]
    metadata_df[[spec$split_group_col]] <- as.character(group_values[group_keep])
  }

  if (nrow(metadata_df) == 0 || length(unique(metadata_df$.xgb_class)) < 2) {
    stop("XGBoost model '", spec$name, "' does not leave at least two classes.", call. = FALSE)
  }

  metadata_df
}

make_xgboost_split <- function(meta, spec) {
  if (nrow(meta) < 2) {
    stop("XGBoost model '", spec$name, "' needs at least two samples for train/test splitting.", call. = FALSE)
  }

  if (is.null(spec$split_group_col)) {
    return(make_xgboost_sample_split(meta, spec))
  }
  make_xgboost_grouped_split(meta, spec)
}

make_xgboost_grouped_split <- function(meta, spec) {
  group_values <- as.character(meta[[spec$split_group_col]])
  groups <- unique(group_values)
  if (length(groups) < 2) {
    stop(
      "XGBoost model '",
      spec$name,
      "' needs at least two groups in ",
      spec$split_group_col,
      " for grouped train/test splitting.",
      call. = FALSE
    )
  }

  test_size <- min(max(1L, ceiling(spec$test_fraction * length(groups))), length(groups) - 1L)
  for (attempt in seq_len(spec$split_max_attempts)) {
    test_groups <- sample(groups, size = test_size)
    train_idx <- !(group_values %in% test_groups)
    test_idx <- group_values %in% test_groups
    if (xgboost_split_has_two_classes(meta$.xgb_class, train_idx, test_idx)) {
      return(list(
        train_idx = train_idx,
        test_idx = test_idx,
        split_group_col = spec$split_group_col,
        test_groups = test_groups
      ))
    }
  }

  stop(
    "Could not create grouped XGBoost train/test split for model '",
    spec$name,
    "' with both classes in train and test after ",
    spec$split_max_attempts,
    " attempts.",
    call. = FALSE
  )
}

make_xgboost_sample_split <- function(meta, spec) {
  sample_indices <- seq_len(nrow(meta))
  test_size <- min(max(1L, ceiling(spec$test_fraction * length(sample_indices))), length(sample_indices) - 1L)
  for (attempt in seq_len(spec$split_max_attempts)) {
    test_samples <- sample(sample_indices, size = test_size)
    test_idx <- sample_indices %in% test_samples
    train_idx <- !test_idx
    if (xgboost_split_has_two_classes(meta$.xgb_class, train_idx, test_idx)) {
      return(list(
        train_idx = train_idx,
        test_idx = test_idx,
        split_group_col = NULL,
        test_samples = meta$.sample_name[test_idx]
      ))
    }
  }

  stop(
    "Could not create sample-level XGBoost train/test split for model '",
    spec$name,
    "' with both classes in train and test after ",
    spec$split_max_attempts,
    " attempts.",
    call. = FALSE
  )
}

xgboost_split_has_two_classes <- function(labels, train_idx, test_idx) {
  length(unique(labels[train_idx])) == 2 && length(unique(labels[test_idx])) == 2
}

make_xgboost_cv_folds <- function(meta, train_idx, spec) {
  train_labels <- meta$.xgb_class[train_idx]
  train_rows <- which(train_idx)
  if (length(unique(train_labels)) < 2) {
    stop("XGBoost model '", spec$name, "' needs both classes in training data.", call. = FALSE)
  }

  if (is.null(spec$split_group_col)) {
    return(make_xgboost_sample_cv_folds(train_labels, spec))
  }
  train_groups <- as.character(meta[[spec$split_group_col]][train_rows])
  make_xgboost_grouped_cv_folds(train_labels, train_groups, spec)
}

make_xgboost_grouped_cv_folds <- function(train_labels, group_values, spec) {
  groups <- unique(group_values)
  k <- min(spec$cv_folds, length(groups))
  if (k < 2) {
    stop("XGBoost model '", spec$name, "' needs at least two training groups for cross-validation.", call. = FALSE)
  }

  train_positions <- seq_along(train_labels)
  for (attempt in seq_len(spec$split_max_attempts)) {
    group_fold <- sample(rep(seq_len(k), length.out = length(groups)))
    names(group_fold) <- groups
    folds <- lapply(seq_len(k), function(fold_id) {
      train_positions[group_values %in% names(group_fold)[group_fold == fold_id]]
    })
    if (xgboost_cv_folds_are_usable(train_labels, train_positions, folds)) {
      return(folds)
    }
  }

  stop(
    "Could not create grouped XGBoost cross-validation folds for model '",
    spec$name,
    "' with both classes in each validation fold after ",
    spec$split_max_attempts,
    " attempts.",
    call. = FALSE
  )
}

make_xgboost_sample_cv_folds <- function(train_labels, spec) {
  train_positions <- seq_along(train_labels)
  k <- min(spec$cv_folds, length(train_positions))
  if (k < 2) {
    stop("XGBoost model '", spec$name, "' needs at least two training samples for cross-validation.", call. = FALSE)
  }

  for (attempt in seq_len(spec$split_max_attempts)) {
    fold_ids <- sample(rep(seq_len(k), length.out = length(train_positions)))
    folds <- lapply(seq_len(k), function(fold_id) train_positions[fold_ids == fold_id])
    if (xgboost_cv_folds_are_usable(train_labels, train_positions, folds)) {
      return(folds)
    }
  }

  stop(
    "Could not create sample-level XGBoost cross-validation folds for model '",
    spec$name,
    "' with both classes in each validation fold after ",
    spec$split_max_attempts,
    " attempts.",
    call. = FALSE
  )
}

xgboost_cv_folds_are_usable <- function(labels, train_rows, folds) {
  all(vapply(folds, function(fold_rows) {
    holdout_labels <- labels[fold_rows]
    fit_labels <- labels[setdiff(train_rows, fold_rows)]
    length(unique(holdout_labels)) == 2 && length(unique(fit_labels)) == 2
  }, logical(1)))
}
