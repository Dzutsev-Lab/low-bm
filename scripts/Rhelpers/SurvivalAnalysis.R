source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

sanitize_survival_path_component <- function(x, fallback = "feature") {
  vapply(as.character(x), function(value) {
    if (length(value) == 0 || is.na(value) || !nzchar(trimws(value))) {
      value <- fallback
    }
    value <- trimws(value)
    value <- gsub("[[:space:]/\\\\:;,*?\"<>|]+", "_", value)
    value <- gsub("[^A-Za-z0-9._+=@-]", "_", value)
    value <- gsub("_+", "_", value)
    value <- sub("^_+", "", value)
    value <- sub("_+$", "", value)
    if (!nzchar(value)) fallback else value
  }, character(1), USE.NAMES = FALSE)
}

normalize_survival_config <- function(config = list(), project_config = list()) {
  config <- config %||% list()
  if (is.null(config$time_col)) config$time_col <- "SurvivalDays"
  if (is.null(config$status_col)) config$status_col <- "SurvivalStatus"
  if (is.null(config$patient_id_col)) config$patient_id_col <- "PatientID"
  if (is.null(config$event_code)) config$event_code <- 1
  if (is.null(config$censor_code)) config$censor_code <- 0
  if (is.null(config$min_n)) config$min_n <- 30
  if (is.null(config$min_events)) config$min_events <- 10
  if (is.null(config$taxa_pseudocount)) config$taxa_pseudocount <- 1e-6

  if (is.null(config$analyses) || length(config$analyses) == 0) {
    stop("survival_analysis must define at least one analysis.", call. = FALSE)
  }

  config$analyses <- lapply(config$analyses, function(spec) {
    spec <- spec %||% list()
    if (is.null(spec$name) || !nzchar(trimws(as.character(spec$name)))) {
      stop("Each survival analysis spec must define a non-empty name.", call. = FALSE)
    }
    if (is.null(spec$covariates)) spec$covariates <- character(0)
    if (is.null(spec$tax_agg_level)) spec$tax_agg_level <- project_config$tax_agg_level %||% "Genus"
    if (is.null(spec$abundance_norm_method)) spec$abundance_norm_method <- "RelAbund"
    if (is.null(spec$cox_feature_transform)) spec$cox_feature_transform <- "log2_pseudocount"
    if (is.null(spec$pcoa_distance)) spec$pcoa_distance <- "bray"
    if (is.null(spec$pcoa_axes)) spec$pcoa_axes <- 3
    if (is.null(spec$alpha_measures)) spec$alpha_measures <- c("Shannon", "Simpson")
    if (is.null(spec$taxa_min_prevalence)) spec$taxa_min_prevalence <- 0.10
    if (is.null(spec$taxa_min_mean_relative_abundance)) spec$taxa_min_mean_relative_abundance <- 1e-5
    if (is.null(spec$km_cutpoint)) spec$km_cutpoint <- "median"
    if (is.null(spec$km_top_n)) spec$km_top_n <- 20
    if (is.null(spec$min_n)) spec$min_n <- config$min_n
    if (is.null(spec$min_events)) spec$min_events <- config$min_events
    if (is.null(spec$taxa_pseudocount)) spec$taxa_pseudocount <- config$taxa_pseudocount
    spec
  })

  config
}

coerce_survival_numeric <- function(x, column, allow_missing = TRUE) {
  x <- standardize_metadata_missing(x)
  missing <- is_metadata_missing_like(x)
  numeric_x <- suppressWarnings(as.numeric(as.character(x)))
  invalid <- !missing & (is.na(numeric_x) | !is.finite(numeric_x))
  if (any(invalid)) {
    stop(
      column,
      " contains non-numeric value(s): ",
      paste(unique(as.character(x[invalid])), collapse = ", "),
      call. = FALSE
    )
  }
  if (!allow_missing && any(missing)) {
    stop(column, " contains missing value(s).", call. = FALSE)
  }
  numeric_x[missing] <- NA_real_
  numeric_x
}

prepare_survival_metadata <- function(metadata_df,
                                      time_col = "SurvivalDays",
                                      status_col = "SurvivalStatus",
                                      patient_id_col = "PatientID",
                                      event_code = 1,
                                      censor_code = 0,
                                      covariates = character(0)) {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  required <- unique(c(patient_id_col, time_col, status_col, covariates))
  fail_missing_columns(names(metadata_df), required, "phyloseq sample_data")

  metadata_df <- standardize_metadata_missing_df(metadata_df, required)
  metadata_df[[patient_id_col]] <- as.character(metadata_df[[patient_id_col]])
  missing_patient <- is_metadata_missing_like(metadata_df[[patient_id_col]])
  if (any(missing_patient)) {
    stop("Patient ID column contains missing value(s): ", patient_id_col, call. = FALSE)
  }

  metadata_df$.survival_time <- coerce_survival_numeric(metadata_df[[time_col]], time_col)
  metadata_df$.survival_status <- coerce_survival_numeric(metadata_df[[status_col]], status_col)

  bad_time <- !is.na(metadata_df$.survival_time) & metadata_df$.survival_time <= 0
  if (any(bad_time)) {
    stop(time_col, " must contain positive survival times.", call. = FALSE)
  }

  allowed_status <- c(as.numeric(event_code), as.numeric(censor_code))
  bad_status <- !is.na(metadata_df$.survival_status) & !metadata_df$.survival_status %in% allowed_status
  if (any(bad_status)) {
    stop(
      status_col,
      " must contain only event/censor code(s): ",
      paste(allowed_status, collapse = ", "),
      call. = FALSE
    )
  }
  metadata_df$.survival_status <- ifelse(
    is.na(metadata_df$.survival_status),
    NA_real_,
    ifelse(metadata_df$.survival_status == as.numeric(event_code), 1, 0)
  )

  metadata_df
}

collapse_metadata_value <- function(x) {
  x <- standardize_metadata_missing(x)
  nonmissing <- x[!is_metadata_missing_like(x)]
  if (length(nonmissing) == 0) {
    return(NA)
  }
  unique_values <- unique(as.character(nonmissing))
  unique_values[[1]]
}

collapse_physeq_by_patient <- function(physeq, patient_id_col = "PatientID") {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), patient_id_col, "phyloseq sample_data")
  metadata_df[[patient_id_col]] <- as.character(standardize_metadata_missing(metadata_df[[patient_id_col]]))
  if (any(is_metadata_missing_like(metadata_df[[patient_id_col]]))) {
    stop("Cannot collapse phyloseq object: patient ID contains missing value(s).", call. = FALSE)
  }

  otu_mat <- otu_samples_by_taxa(physeq)
  group_factor <- factor(metadata_df[[patient_id_col]], levels = unique(metadata_df[[patient_id_col]]))
  group_counts <- table(group_factor)
  patient_otu <- rowsum(otu_mat, group = group_factor, reorder = FALSE)
  patient_otu <- sweep(patient_otu, 1, as.numeric(group_counts[rownames(patient_otu)]), FUN = "/")

  patient_metadata <- do.call(rbind, lapply(levels(group_factor), function(patient) {
    rows <- metadata_df[group_factor == patient, , drop = FALSE]
    values <- lapply(rows, collapse_metadata_value)
    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  rownames(patient_metadata) <- levels(group_factor)
  patient_metadata[[patient_id_col]] <- rownames(patient_metadata)

  phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(patient_otu), taxa_are_rows = FALSE),
    phyloseq::sample_data(patient_metadata),
    phyloseq::tax_table(physeq)
  )
}

make_survival_base_df <- function(patient_physeq,
                                  time_col = "SurvivalDays",
                                  status_col = "SurvivalStatus",
                                  patient_id_col = "PatientID",
                                  event_code = 1,
                                  censor_code = 0,
                                  covariates = character(0)) {
  metadata_df <- as.data.frame(phyloseq::sample_data(patient_physeq), stringsAsFactors = FALSE)
  metadata_df <- prepare_survival_metadata(
    metadata_df,
    time_col = time_col,
    status_col = status_col,
    patient_id_col = patient_id_col,
    event_code = event_code,
    censor_code = censor_code,
    covariates = covariates
  )
  base_cols <- unique(c(patient_id_col, ".survival_time", ".survival_status", covariates))
  base_df <- metadata_df[, base_cols, drop = FALSE]
  colnames(base_df)[colnames(base_df) == patient_id_col] <- "PatientID"
  base_df
}

build_alpha_features <- function(patient_physeq, measures = c("Shannon", "Simpson")) {
  if (is.null(measures) || length(measures) == 0) {
    return(list(matrix = data.frame(row.names = phyloseq::sample_names(patient_physeq)), map = data.frame()))
  }
  alpha_df <- phyloseq::estimate_richness(patient_physeq, measures = measures)
  alpha_df <- as.data.frame(alpha_df, stringsAsFactors = FALSE)
  keep <- intersect(measures, names(alpha_df))
  alpha_df <- alpha_df[, keep, drop = FALSE]
  colnames(alpha_df) <- paste0("alpha_", colnames(alpha_df))
  feature_map <- data.frame(
    feature = colnames(alpha_df),
    feature_family = "alpha",
    label = sub("^alpha_", "", colnames(alpha_df)),
    stringsAsFactors = FALSE
  )
  list(matrix = alpha_df, map = feature_map)
}

build_pcoa_features <- function(patient_physeq,
                                tax_agg_level = "Genus",
                                abundance_norm_method = "RelAbund",
                                pcoa_distance = "bray",
                                pcoa_axes = 3,
                                pseudocount = 1) {
  if (is.null(pcoa_axes) || as.integer(pcoa_axes) <= 0 || phyloseq::nsamples(patient_physeq) < 3) {
    return(list(matrix = data.frame(row.names = phyloseq::sample_names(patient_physeq)), map = data.frame()))
  }
  norm_physeq <- patient_physeq |>
    tax_glom_rename(tax_agg_level) |>
    counts_normalization(norm_method = abundance_norm_method, pseudocount = pseudocount) |>
    prune_empty_physeq()

  dist_obj <- phyloseq::distance(norm_physeq, method = as.character(pcoa_distance))
  axes <- min(as.integer(pcoa_axes), phyloseq::nsamples(norm_physeq) - 1)
  pcoa <- stats::cmdscale(as.dist(dist_obj), eig = TRUE, k = axes)
  pcoa_df <- as.data.frame(pcoa$points, stringsAsFactors = FALSE)
  colnames(pcoa_df) <- paste0("pcoa_", sanitize_survival_path_component(pcoa_distance), "_axis", seq_len(ncol(pcoa_df)))
  rownames(pcoa_df) <- rownames(as.matrix(dist_obj))
  feature_map <- data.frame(
    feature = colnames(pcoa_df),
    feature_family = "pcoa",
    label = paste0("PCoA ", seq_len(ncol(pcoa_df)), " (", pcoa_distance, ")"),
    stringsAsFactors = FALSE
  )
  list(matrix = pcoa_df, map = feature_map)
}

relative_abundance_physeq <- function(patient_physeq, tax_agg_level = "Genus") {
  patient_physeq |>
    tax_glom_rename(tax_agg_level) |>
    counts_normalization(norm_method = "RelAbund") |>
    prune_empty_physeq()
}

filter_taxa_by_relative_abundance <- function(rel_physeq,
                                              min_prevalence = 0.10,
                                              min_mean_relative_abundance = 1e-5) {
  rel_mat <- otu_samples_by_taxa(rel_physeq)
  prevalence <- colMeans(rel_mat > 0)
  mean_rel_abund <- colMeans(rel_mat)
  keep <- prevalence >= as.numeric(min_prevalence) &
    mean_rel_abund >= as.numeric(min_mean_relative_abundance)
  list(
    taxa = colnames(rel_mat)[keep],
    stats = data.frame(
      taxon = colnames(rel_mat),
      prevalence = prevalence,
      mean_relative_abundance = mean_rel_abund,
      retained = keep,
      stringsAsFactors = FALSE
    )
  )
}

transform_taxa_features <- function(rel_mat,
                                    transform = "log2_pseudocount",
                                    pseudocount = 1e-6) {
  transform <- as.character(transform %||% "log2_pseudocount")
  if (identical(transform, "none")) {
    return(rel_mat)
  }
  if (identical(transform, "log2_pseudocount")) {
    return(log2(rel_mat + as.numeric(pseudocount)))
  }
  stop("Unknown Cox feature transform: ", transform, call. = FALSE)
}

build_taxa_features <- function(patient_physeq,
                                tax_agg_level = "Genus",
                                min_prevalence = 0.10,
                                min_mean_relative_abundance = 1e-5,
                                cox_feature_transform = "log2_pseudocount",
                                taxa_pseudocount = 1e-6) {
  rel_physeq <- relative_abundance_physeq(patient_physeq, tax_agg_level)
  filter_info <- filter_taxa_by_relative_abundance(
    rel_physeq,
    min_prevalence = min_prevalence,
    min_mean_relative_abundance = min_mean_relative_abundance
  )
  if (length(filter_info$taxa) == 0) {
    empty <- data.frame(row.names = phyloseq::sample_names(patient_physeq))
    return(list(matrix = empty, map = data.frame(), filter_stats = filter_info$stats))
  }

  rel_mat <- otu_samples_by_taxa(rel_physeq)[, filter_info$taxa, drop = FALSE]
  taxa_mat <- transform_taxa_features(
    rel_mat,
    transform = cox_feature_transform,
    pseudocount = taxa_pseudocount
  )
  taxa_df <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  feature_names <- paste0("taxa_", make.unique(sanitize_survival_path_component(colnames(taxa_df))))
  feature_map <- data.frame(
    feature = feature_names,
    feature_family = "taxa",
    label = colnames(taxa_df),
    stringsAsFactors = FALSE
  )
  colnames(taxa_df) <- feature_names
  list(matrix = taxa_df, map = feature_map, filter_stats = filter_info$stats)
}

build_patient_feature_matrix <- function(patient_physeq,
                                         spec,
                                         survival_config,
                                         project_config = list()) {
  covariates <- as.character(unlist(spec$covariates %||% character(0), use.names = FALSE))
  base_df <- make_survival_base_df(
    patient_physeq,
    time_col = survival_config$time_col,
    status_col = survival_config$status_col,
    patient_id_col = survival_config$patient_id_col,
    event_code = survival_config$event_code,
    censor_code = survival_config$censor_code,
    covariates = covariates
  )
  rownames(base_df) <- base_df$PatientID

  alpha <- build_alpha_features(patient_physeq, measures = spec$alpha_measures)
  pcoa <- build_pcoa_features(
    patient_physeq,
    tax_agg_level = spec$tax_agg_level,
    abundance_norm_method = spec$abundance_norm_method,
    pcoa_distance = spec$pcoa_distance,
    pcoa_axes = spec$pcoa_axes,
    pseudocount = project_config$pseudocount %||% 1
  )
  taxa <- build_taxa_features(
    patient_physeq,
    tax_agg_level = spec$tax_agg_level,
    min_prevalence = spec$taxa_min_prevalence,
    min_mean_relative_abundance = spec$taxa_min_mean_relative_abundance,
    cox_feature_transform = spec$cox_feature_transform,
    taxa_pseudocount = spec$taxa_pseudocount
  )

  feature_parts <- list(alpha$matrix, pcoa$matrix, taxa$matrix)
  aligned_parts <- lapply(feature_parts, function(part) {
    part[rownames(base_df), , drop = FALSE]
  })
  feature_df <- cbind(base_df, do.call(cbind, aligned_parts))
  feature_map <- rbind(alpha$map, pcoa$map, taxa$map)
  rownames(feature_df) <- NULL

  list(
    patient_features = feature_df,
    feature_map = feature_map,
    taxa_filter_stats = taxa$filter_stats
  )
}

coerce_survival_model_columns <- function(model_df, covariates = character(0)) {
  model_df$.feature_z <- as.numeric(scale(as.numeric(model_df$.feature)))
  if (any(is.na(model_df$.feature_z)) || stats::sd(model_df$.feature, na.rm = TRUE) == 0) {
    stop("feature has zero variance after complete-case filtering", call. = FALSE)
  }

  alias_map <- data.frame(
    original = covariates,
    alias = if (length(covariates) == 0) character(0) else paste0(".cov", seq_along(covariates)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(covariates)) {
    covariate <- covariates[[i]]
    alias <- alias_map$alias[[i]]
    value <- standardize_metadata_missing(model_df[[covariate]])
    missing <- is_metadata_missing_like(value)
    if (any(missing)) {
      stop("internal error: missing covariate after complete-case filtering", call. = FALSE)
    }
    numeric_value <- suppressWarnings(as.numeric(as.character(value)))
    numeric_invalid <- is.na(numeric_value) | !is.finite(numeric_value)
    if (!any(numeric_invalid)) {
      if (stats::sd(numeric_value) == 0) {
        stop("covariate has fewer than two usable values: ", covariate, call. = FALSE)
      }
      model_df[[alias]] <- numeric_value
    } else {
      factor_value <- factor(as.character(value))
      if (length(levels(factor_value)) < 2) {
        stop("covariate has fewer than two usable levels: ", covariate, call. = FALSE)
      }
      model_df[[alias]] <- factor_value
    }
  }

  list(data = model_df, covariate_aliases = alias_map$alias)
}

fit_cox_feature <- function(patient_features,
                            feature,
                            feature_family,
                            feature_label = feature,
                            covariates = character(0),
                            min_n = 30,
                            min_events = 10) {
  required <- unique(c(".survival_time", ".survival_status", feature, covariates))
  fail_missing_columns(names(patient_features), required, "patient feature matrix")
  model_df <- patient_features[, required, drop = FALSE]
  names(model_df)[names(model_df) == feature] <- ".feature"
  model_df <- standardize_metadata_missing_df(model_df)

  complete <- stats::complete.cases(model_df)
  n_total <- nrow(model_df)
  n_used <- sum(complete)
  dropped_missing <- n_total - n_used
  if (n_used < as.integer(min_n)) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, "too_few_complete_cases", n_total, n_used, NA, dropped_missing)
    ))
  }

  model_df <- model_df[complete, , drop = FALSE]
  model_df$.survival_time <- as.numeric(model_df$.survival_time)
  model_df$.survival_status <- as.numeric(model_df$.survival_status)
  events_used <- sum(model_df$.survival_status == 1)
  if (events_used < as.integer(min_events)) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, "too_few_events", n_total, n_used, events_used, dropped_missing)
    ))
  }

  coerced <- tryCatch(
    coerce_survival_model_columns(model_df, covariates = covariates),
    error = function(e) e
  )
  if (inherits(coerced, "error")) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, conditionMessage(coerced), n_total, n_used, events_used, dropped_missing)
    ))
  }

  model_df <- coerced$data
  terms <- c(".feature_z", coerced$covariate_aliases)
  formula <- stats::as.formula(paste("survival::Surv(.survival_time, .survival_status) ~", paste(terms, collapse = " + ")))
  fit <- tryCatch(
    survival::coxph(formula, data = model_df),
    error = function(e) e,
    warning = function(w) w
  )
  if (inherits(fit, "error") || inherits(fit, "warning")) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, paste0("coxph failed: ", conditionMessage(fit)), n_total, n_used, events_used, dropped_missing)
    ))
  }

  fit_summary <- summary(fit)
  if (!".feature_z" %in% rownames(fit_summary$coefficients)) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, "coxph omitted feature coefficient", n_total, n_used, events_used, dropped_missing)
    ))
  }

  coef_row <- fit_summary$coefficients[".feature_z", , drop = FALSE]
  conf_row <- fit_summary$conf.int[".feature_z", , drop = FALSE]
  ph_p <- NA_real_
  global_ph_p <- NA_real_
  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (!is.null(zph) && ".feature_z" %in% rownames(zph$table)) {
    ph_p <- zph$table[".feature_z", "p"]
    global_ph_p <- zph$table["GLOBAL", "p"]
  }

  result <- data.frame(
    feature = feature,
    feature_family = feature_family,
    label = feature_label,
    hazard_ratio = unname(conf_row[, "exp(coef)"]),
    conf_low = unname(conf_row[, "lower .95"]),
    conf_high = unname(conf_row[, "upper .95"]),
    coef = unname(coef_row[, "coef"]),
    se = unname(coef_row[, "se(coef)"]),
    p = unname(coef_row[, "Pr(>|z|)"]),
    n_total = n_total,
    n_used = n_used,
    events_used = events_used,
    dropped_missing = dropped_missing,
    ph_p = ph_p,
    global_ph_p = global_ph_p,
    stringsAsFactors = FALSE
  )
  list(result = result, skip = NULL)
}

survival_skip_row <- function(feature,
                              feature_family,
                              feature_label,
                              reason,
                              n_total,
                              n_used,
                              events_used,
                              dropped_missing) {
  data.frame(
    feature = feature,
    feature_family = feature_family,
    label = feature_label,
    reason = reason,
    n_total = n_total,
    n_used = n_used,
    events_used = events_used,
    dropped_missing = dropped_missing,
    stringsAsFactors = FALSE
  )
}

run_survival_cox_models <- function(patient_features,
                                    feature_map,
                                    covariates = character(0),
                                    min_n = 30,
                                    min_events = 10) {
  if (nrow(feature_map) == 0) {
    return(list(results = data.frame(), skipped = data.frame()))
  }

  fits <- lapply(seq_len(nrow(feature_map)), function(i) {
    fit_cox_feature(
      patient_features = patient_features,
      feature = feature_map$feature[[i]],
      feature_family = feature_map$feature_family[[i]],
      feature_label = feature_map$label[[i]],
      covariates = covariates,
      min_n = min_n,
      min_events = min_events
    )
  })

  results <- do.call(rbind, lapply(fits, `[[`, "result"))
  skipped <- do.call(rbind, lapply(fits, `[[`, "skip"))
  if (is.null(results)) results <- data.frame()
  if (is.null(skipped)) skipped <- data.frame()

  if (nrow(results) > 0) {
    results$fdr <- NA_real_
    for (family in unique(results$feature_family)) {
      idx <- results$feature_family == family
      results$fdr[idx] <- stats::p.adjust(results$p[idx], method = "BH")
    }
    results <- results[order(results$fdr, results$p), , drop = FALSE]
  }

  list(results = results, skipped = skipped)
}

model_missingness_table <- function(results, skipped) {
  keep_cols <- c("feature", "feature_family", "label", "n_total", "n_used", "events_used", "dropped_missing")
  parts <- list()
  if (nrow(results) > 0) {
    parts[[length(parts) + 1]] <- results[, keep_cols, drop = FALSE]
  }
  if (nrow(skipped) > 0) {
    parts[[length(parts) + 1]] <- skipped[, keep_cols, drop = FALSE]
  }
  if (length(parts) == 0) {
    return(data.frame())
  }
  out <- do.call(rbind, parts)
  out[order(out$feature_family, out$feature), , drop = FALSE]
}

km_group_data <- function(patient_features,
                          feature,
                          covariates = character(0),
                          cutpoint = "median") {
  required <- unique(c(".survival_time", ".survival_status", feature, covariates))
  model_df <- patient_features[, required, drop = FALSE]
  names(model_df)[names(model_df) == feature] <- ".feature"
  model_df <- standardize_metadata_missing_df(model_df)
  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  if (nrow(model_df) < 2 || stats::sd(as.numeric(model_df$.feature), na.rm = TRUE) == 0) {
    stop("feature cannot be split into Kaplan-Meier groups", call. = FALSE)
  }

  if (!identical(as.character(cutpoint), "median")) {
    stop("Only median Kaplan-Meier cutpoints are supported in v1.", call. = FALSE)
  }
  median_value <- stats::median(as.numeric(model_df$.feature), na.rm = TRUE)
  model_df$.km_group <- ifelse(as.numeric(model_df$.feature) >= median_value, "High", "Low")
  model_df$.km_group <- factor(model_df$.km_group, levels = c("Low", "High"))
  if (length(unique(model_df$.km_group)) < 2) {
    stop("median split produced fewer than two Kaplan-Meier groups", call. = FALSE)
  }
  model_df$.survival_time <- as.numeric(model_df$.survival_time)
  model_df$.survival_status <- as.numeric(model_df$.survival_status)
  attr(model_df, "cutpoint") <- median_value
  model_df
}

survfit_to_plot_df <- function(fit) {
  fit_summary <- summary(fit)
  plot_df <- data.frame(
    time = fit_summary$time,
    survival = fit_summary$surv,
    lower = fit_summary$lower,
    upper = fit_summary$upper,
    strata = fit_summary$strata,
    stringsAsFactors = FALSE
  )
  if (nrow(plot_df) == 0) {
    return(plot_df)
  }
  starts <- do.call(rbind, lapply(unique(plot_df$strata), function(stratum) {
    data.frame(time = 0, survival = 1, lower = 1, upper = 1, strata = stratum, stringsAsFactors = FALSE)
  }))
  rbind(starts, plot_df)
}

km_logrank_p <- function(km_df) {
  sd <- survival::survdiff(survival::Surv(.survival_time, .survival_status) ~ .km_group, data = km_df)
  stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
}

plot_km_feature <- function(patient_features,
                            feature,
                            feature_label,
                            covariates = character(0),
                            cutpoint = "median",
                            cox_result = NULL) {
  km_df <- km_group_data(patient_features, feature, covariates = covariates, cutpoint = cutpoint)
  fit <- survival::survfit(survival::Surv(.survival_time, .survival_status) ~ .km_group, data = km_df)
  plot_df <- survfit_to_plot_df(fit)
  logrank_p <- km_logrank_p(km_df)
  group_counts <- aggregate(
    km_df$.survival_status,
    by = list(group = km_df$.km_group),
    FUN = function(x) paste0(length(x), " n / ", sum(x == 1), " events")
  )
  group_label <- paste(paste(group_counts$group, group_counts$x, sep = ": "), collapse = "\n")
  hr_label <- ""
  if (!is.null(cox_result) && nrow(cox_result) > 0) {
    hr_label <- sprintf(
      "\nCox HR per 1 SD: %.2f (95%% CI %.2f-%.2f)",
      cox_result$hazard_ratio[[1]],
      cox_result$conf_low[[1]],
      cox_result$conf_high[[1]]
    )
  }
  subtitle <- paste0(
    "Median split at ",
    signif(attr(km_df, "cutpoint"), 4),
    " | log-rank p = ",
    signif(logrank_p, 3),
    hr_label,
    "\n",
    group_label
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = survival, color = strata)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0("Kaplan-Meier: ", feature_label),
      subtitle = subtitle,
      x = "Survival days",
      y = "Survival probability",
      color = "Group"
    ) +
    ggplot2::theme_bw()
}
