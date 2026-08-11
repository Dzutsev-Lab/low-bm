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

unique_survival_column_name <- function(existing, proposed) {
  make.unique(c(existing, proposed))[[length(existing) + 1]]
}

survival_strata_enabled <- function(sample_strata_col) {
  !is.null(sample_strata_col) &&
    length(sample_strata_col) > 0 &&
    !is.na(sample_strata_col[[1]]) &&
    nzchar(trimws(as.character(sample_strata_col[[1]])))
}

normalize_survival_methods <- function(methods = NULL) {
  if (is.null(methods)) {
    return("cox")
  }
  methods <- as.character(unlist(methods, use.names = FALSE))
  methods <- methods[!is.na(methods) & nzchar(trimws(methods))]
  methods <- tolower(trimws(methods))
  methods[methods %in% c("coda", "coda4microbiome")] <- "coda4microbiome"
  methods[methods %in% c("standard_cox", "standard-cox", "per_taxon_cox", "per-taxon-cox")] <- "cox"
  methods <- unique(methods)
  valid <- c("cox", "coda4microbiome")
  invalid <- setdiff(methods, valid)
  if (length(invalid) > 0) {
    stop("Unsupported survival analysis method(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (length(methods) == 0) {
    stop("survival_analysis.methods must include at least one method.", call. = FALSE)
  }
  methods
}

survival_method_enabled <- function(spec, method) {
  method <- normalize_survival_methods(method)
  any(normalize_survival_methods(spec$methods) %in% method)
}

normalize_coda_options <- function(options = list()) {
  defaults <- list(
    lambda = "lambda.1se",
    alpha = 0.9,
    nfolds = 10,
    nvar = NULL,
    coef_threshold = 0,
    seed = 12345,
    show_plots = FALSE
  )
  options <- utils::modifyList(defaults, options %||% list(), keep.null = FALSE)
  options$alpha <- as.numeric(options$alpha)
  options$nfolds <- as.integer(options$nfolds)
  options$coef_threshold <- as.numeric(options$coef_threshold)
  options$seed <- as.integer(options$seed)
  if (!is.null(options$nvar)) options$nvar <- as.integer(options$nvar)
  options$show_plots <- if (is.logical(options$show_plots)) {
    isTRUE(options$show_plots)
  } else {
    !tolower(trimws(as.character(options$show_plots))) %in% c("0", "false", "f", "no", "n", "")
  }
  if (!is.finite(options$alpha) || options$alpha < 0 || options$alpha > 1) {
    stop("survival_analysis.coda4microbiome.alpha must be between 0 and 1.", call. = FALSE)
  }
  if (length(options$nfolds) == 0 || is.na(options$nfolds) || options$nfolds < 1) {
    stop("survival_analysis.coda4microbiome.nfolds must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(options$coef_threshold) || options$coef_threshold < 0) {
    stop("survival_analysis.coda4microbiome.coef_threshold must be non-negative.", call. = FALSE)
  }
  if (length(options$seed) == 0 || is.na(options$seed)) {
    stop("survival_analysis.coda4microbiome.seed must be an integer.", call. = FALSE)
  }
  options
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
  if (is.null(config$norm_method)) config$norm_method <- project_config$norm_method %||% "noNorm"
  if (is.null(config$pseudocount)) config$pseudocount <- project_config$pseudocount %||% 1
  if (is.null(config$tax_agg_level)) config$tax_agg_level <- project_config$tax_agg_level %||% "Genus"
  if (!"sample_strata_col" %in% names(config)) config$sample_strata_col <- "SampleType"
  config$methods <- normalize_survival_methods(config$methods)
  config$coda4microbiome <- normalize_coda_options(config$coda4microbiome)

  if (is.null(config$analyses) || length(config$analyses) == 0) {
    stop("survival_analysis must define at least one analysis.", call. = FALSE)
  }

  config$analyses <- lapply(config$analyses, function(spec) {
    spec <- spec %||% list()
    if (is.null(spec$name) || !nzchar(trimws(as.character(spec$name)))) {
      stop("Each survival analysis spec must define a non-empty name.", call. = FALSE)
    }
    if (is.null(spec$covariates)) spec$covariates <- character(0)
    if (is.null(spec$tax_agg_level)) spec$tax_agg_level <- config$tax_agg_level
    if (is.null(spec$norm_method)) spec$norm_method <- config$norm_method
    if (is.null(spec$pseudocount)) spec$pseudocount <- config$pseudocount
    if (!"sample_strata_col" %in% names(spec)) spec$sample_strata_col <- config$sample_strata_col
    if (is.null(spec$methods)) spec$methods <- config$methods
    spec$methods <- normalize_survival_methods(spec$methods)
    spec$coda4microbiome <- normalize_coda_options(utils::modifyList(
      config$coda4microbiome,
      spec$coda4microbiome %||% list(),
      keep.null = FALSE
    ))
    if (is.null(spec$taxa_min_prevalence)) spec$taxa_min_prevalence <- 0.10
    if (is.null(spec$taxa_min_mean_relative_abundance)) spec$taxa_min_mean_relative_abundance <- 1e-5
    if (is.null(spec$min_n)) spec$min_n <- config$min_n
    if (is.null(spec$min_events)) spec$min_events <- config$min_events
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

make_survival_base_df_from_metadata <- function(metadata_df,
                                                time_col = "SurvivalDays",
                                                status_col = "SurvivalStatus",
                                                patient_id_col = "PatientID",
                                                event_code = 1,
                                                censor_code = 0,
                                                covariates = character(0)) {
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

make_survival_base_df <- function(patient_physeq,
                                  time_col = "SurvivalDays",
                                  status_col = "SurvivalStatus",
                                  patient_id_col = "PatientID",
                                  event_code = 1,
                                  censor_code = 0,
                                  covariates = character(0)) {
  metadata_df <- as.data.frame(phyloseq::sample_data(patient_physeq), stringsAsFactors = FALSE)
  make_survival_base_df_from_metadata(
    metadata_df,
    time_col = time_col,
    status_col = status_col,
    patient_id_col = patient_id_col,
    event_code = event_code,
    censor_code = censor_code,
    covariates = covariates
  )
}

build_patient_sample_summary <- function(metadata_df,
                                         patient_id_col = "PatientID",
                                         sample_strata_col = "SampleType") {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), patient_id_col, "phyloseq sample_data")
  metadata_df[[patient_id_col]] <- as.character(standardize_metadata_missing(metadata_df[[patient_id_col]]))
  if (any(is_metadata_missing_like(metadata_df[[patient_id_col]]))) {
    stop("Cannot summarize survival samples: patient ID contains missing value(s).", call. = FALSE)
  }

  patient_factor <- factor(metadata_df[[patient_id_col]], levels = unique(metadata_df[[patient_id_col]]))
  patient_ids <- levels(patient_factor)
  patient_counts <- table(patient_factor)
  patient_metadata <- do.call(rbind, lapply(patient_ids, function(patient) {
    rows <- metadata_df[patient_factor == patient, , drop = FALSE]
    values <- lapply(rows, collapse_metadata_value)
    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  rownames(patient_metadata) <- patient_ids
  patient_metadata[[patient_id_col]] <- rownames(patient_metadata)

  audit_df <- data.frame(
    n_samples_collapsed = as.integer(patient_counts[patient_ids]),
    stringsAsFactors = FALSE,
    row.names = patient_ids
  )
  observed_strata <- character(0)
  strata_col <- NULL

  if (survival_strata_enabled(sample_strata_col)) {
    strata_col <- as.character(sample_strata_col[[1]])
    fail_missing_columns(names(metadata_df), strata_col, "phyloseq sample_data")
    strata_values <- as.character(standardize_metadata_missing(metadata_df[[strata_col]]))
    missing_strata <- is_metadata_missing_like(strata_values)
    strata_values[missing_strata] <- NA_character_
    observed_strata <- sort(unique(strata_values[!is.na(strata_values)]))
    strata_prefix <- sanitize_survival_path_component(strata_col, fallback = "strata")
    collapsed_col <- paste0("collapsed_", strata_prefix)
    audit_df[[collapsed_col]] <- vapply(patient_ids, function(patient) {
      values <- sort(unique(strata_values[patient_factor == patient & !is.na(strata_values)]))
      paste(values, collapse = ";")
    }, character(1), USE.NAMES = FALSE)

    for (stratum in observed_strata) {
      proposed_col <- paste0("n_samples_", strata_prefix, "_", sanitize_survival_path_component(stratum, fallback = "stratum"))
      count_col <- unique_survival_column_name(names(audit_df), proposed_col)
      audit_df[[count_col]] <- vapply(patient_ids, function(patient) {
        sum(patient_factor == patient & !is.na(strata_values) & strata_values == stratum)
      }, integer(1), USE.NAMES = FALSE)
    }
  }

  list(
    metadata = patient_metadata,
    audit = audit_df,
    patient_ids = patient_ids,
    sample_strata_col = strata_col,
    observed_strata = observed_strata
  )
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

empty_survival_feature_map <- function() {
  data.frame(
    feature = character(0),
    feature_family = character(0),
    label = character(0),
    sample_strata_col = character(0),
    sample_stratum = character(0),
    stringsAsFactors = FALSE
  )
}

empty_taxa_filter_stats <- function() {
  data.frame(
    sample_strata_col = character(0),
    sample_stratum = character(0),
    taxon = character(0),
    prevalence = numeric(0),
    mean_relative_abundance = numeric(0),
    retained = logical(0),
    stringsAsFactors = FALSE
  )
}

align_taxa_matrix_to_patients <- function(taxa_df, patient_ids) {
  if (ncol(taxa_df) == 0) {
    return(data.frame(row.names = patient_ids))
  }
  aligned <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = ncol(taxa_df),
    dimnames = list(patient_ids, colnames(taxa_df))
  )
  shared <- intersect(patient_ids, rownames(taxa_df))
  if (length(shared) > 0) {
    aligned[shared, ] <- as.matrix(taxa_df[shared, , drop = FALSE])
  }
  as.data.frame(aligned, stringsAsFactors = FALSE, check.names = FALSE)
}

combine_taxa_feature_parts <- function(parts, patient_ids) {
  if (length(parts) == 0) {
    return(list(
      matrix = data.frame(row.names = patient_ids),
      map = empty_survival_feature_map(),
      filter_stats = empty_taxa_filter_stats()
    ))
  }

  matrices <- lapply(parts, `[[`, "matrix")
  feature_maps <- lapply(parts, `[[`, "map")
  filter_stats <- lapply(parts, `[[`, "filter_stats")
  feature_df <- do.call(cbind, matrices)
  feature_map <- do.call(rbind, feature_maps)
  taxa_stats <- do.call(rbind, filter_stats)

  if (ncol(feature_df) > 0) {
    feature_names <- make.unique(colnames(feature_df))
    colnames(feature_df) <- feature_names
    feature_map$feature <- feature_names
  }

  list(matrix = feature_df, map = feature_map, filter_stats = taxa_stats)
}

build_taxa_features <- function(patient_physeq,
                                tax_agg_level = "Genus",
                                norm_method = "noNorm",
                                pseudocount = 1,
                                min_prevalence = 0.10,
                                min_mean_relative_abundance = 1e-5,
                                sample_strata_col = NA_character_,
                                sample_stratum = NA_character_,
                                feature_prefix = "taxa",
                                label_prefix = NULL) {
  rel_physeq <- relative_abundance_physeq(patient_physeq, tax_agg_level)
  filter_info <- filter_taxa_by_relative_abundance(
    rel_physeq,
    min_prevalence = min_prevalence,
    min_mean_relative_abundance = min_mean_relative_abundance
  )
  filter_stats <- filter_info$stats
  filter_stats <- cbind(
    data.frame(
      sample_strata_col = rep(as.character(sample_strata_col)[[1]], nrow(filter_stats)),
      sample_stratum = rep(as.character(sample_stratum)[[1]], nrow(filter_stats)),
      stringsAsFactors = FALSE
    ),
    filter_stats
  )
  if (length(filter_info$taxa) == 0) {
    empty <- data.frame(row.names = phyloseq::sample_names(patient_physeq))
    return(list(matrix = empty, map = empty_survival_feature_map(), filter_stats = filter_stats))
  }

  model_physeq <- patient_physeq |>
    tax_glom_rename(tax_agg_level) |>
    prune_empty_physeq() |>
    counts_normalization(norm_method = norm_method, pseudocount = pseudocount)
  model_mat <- otu_samples_by_taxa(model_physeq)
  taxa_mat <- model_mat[, filter_info$taxa, drop = FALSE]
  taxa_df <- as.data.frame(taxa_mat, stringsAsFactors = FALSE)
  feature_names <- paste0(feature_prefix, "_", make.unique(sanitize_survival_path_component(colnames(taxa_df))))
  labels <- colnames(taxa_df)
  if (!is.null(label_prefix) && nzchar(label_prefix)) {
    labels <- paste0(label_prefix, " | ", labels)
  }
  feature_map <- data.frame(
    feature = feature_names,
    feature_family = "taxa",
    label = labels,
    sample_strata_col = rep(as.character(sample_strata_col)[[1]], length(feature_names)),
    sample_stratum = rep(as.character(sample_stratum)[[1]], length(feature_names)),
    stringsAsFactors = FALSE
  )
  colnames(taxa_df) <- feature_names
  list(matrix = taxa_df, map = feature_map, filter_stats = filter_stats)
}

build_stratified_taxa_features <- function(sample_physeq,
                                           patient_ids,
                                           patient_id_col = "PatientID",
                                           sample_strata_col = "SampleType",
                                           observed_strata = NULL,
                                           tax_agg_level = "Genus",
                                           norm_method = "noNorm",
                                           pseudocount = 1,
                                           min_prevalence = 0.10,
                                           min_mean_relative_abundance = 1e-5) {
  if (!survival_strata_enabled(sample_strata_col)) {
    patient_physeq <- collapse_physeq_by_patient(sample_physeq, patient_id_col = patient_id_col)
    taxa <- build_taxa_features(
      patient_physeq,
      tax_agg_level = tax_agg_level,
      norm_method = norm_method,
      pseudocount = pseudocount,
      min_prevalence = min_prevalence,
      min_mean_relative_abundance = min_mean_relative_abundance
    )
    taxa$matrix <- align_taxa_matrix_to_patients(taxa$matrix, patient_ids)
    return(taxa)
  }

  strata_col <- as.character(sample_strata_col[[1]])
  metadata_df <- as.data.frame(phyloseq::sample_data(sample_physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), c(patient_id_col, strata_col), "phyloseq sample_data")
  strata_values <- as.character(standardize_metadata_missing(metadata_df[[strata_col]]))
  strata_values[is_metadata_missing_like(strata_values)] <- NA_character_
  if (is.null(observed_strata)) {
    observed_strata <- sort(unique(strata_values[!is.na(strata_values)]))
  }
  if (length(observed_strata) == 0) {
    return(list(
      matrix = data.frame(row.names = patient_ids),
      map = empty_survival_feature_map(),
      filter_stats = empty_taxa_filter_stats()
    ))
  }

  strata_prefix <- sanitize_survival_path_component(strata_col, fallback = "strata")
  parts <- lapply(observed_strata, function(stratum) {
    keep <- !is.na(strata_values) & strata_values == stratum
    stratum_physeq <- phyloseq::prune_samples(keep, sample_physeq)
    patient_stratum_physeq <- collapse_physeq_by_patient(stratum_physeq, patient_id_col = patient_id_col)
    stratum_prefix <- paste0("taxa_", strata_prefix, "_", sanitize_survival_path_component(stratum, fallback = "stratum"))
    taxa <- build_taxa_features(
      patient_stratum_physeq,
      tax_agg_level = tax_agg_level,
      norm_method = norm_method,
      pseudocount = pseudocount,
      min_prevalence = min_prevalence,
      min_mean_relative_abundance = min_mean_relative_abundance,
      sample_strata_col = strata_col,
      sample_stratum = stratum,
      feature_prefix = stratum_prefix,
      label_prefix = paste0(strata_col, "=", stratum)
    )
    taxa$matrix <- align_taxa_matrix_to_patients(taxa$matrix, patient_ids)
    taxa
  })

  combine_taxa_feature_parts(parts, patient_ids)
}

build_patient_feature_matrix <- function(sample_physeq,
                                         spec,
                                         survival_config,
                                         project_config = list()) {
  covariates <- as.character(unlist(spec$covariates %||% character(0), use.names = FALSE))
  sample_summary <- build_patient_sample_summary(
    as.data.frame(phyloseq::sample_data(sample_physeq), stringsAsFactors = FALSE),
    patient_id_col = survival_config$patient_id_col,
    sample_strata_col = spec$sample_strata_col
  )
  base_df <- make_survival_base_df_from_metadata(
    sample_summary$metadata,
    time_col = survival_config$time_col,
    status_col = survival_config$status_col,
    patient_id_col = survival_config$patient_id_col,
    event_code = survival_config$event_code,
    censor_code = survival_config$censor_code,
    covariates = covariates
  )
  rownames(base_df) <- base_df$PatientID
  audit_df <- sample_summary$audit[rownames(base_df), , drop = FALSE]

  taxa <- build_stratified_taxa_features(
    sample_physeq,
    patient_ids = rownames(base_df),
    patient_id_col = survival_config$patient_id_col,
    sample_strata_col = spec$sample_strata_col,
    observed_strata = sample_summary$observed_strata,
    tax_agg_level = spec$tax_agg_level,
    norm_method = spec$norm_method,
    pseudocount = spec$pseudocount,
    min_prevalence = spec$taxa_min_prevalence,
    min_mean_relative_abundance = spec$taxa_min_mean_relative_abundance
  )

  taxa_df <- taxa$matrix[rownames(base_df), , drop = FALSE]
  feature_df <- cbind(base_df, audit_df, taxa_df)
  feature_map <- taxa$map
  rownames(feature_df) <- NULL

  list(
    patient_features = feature_df,
    feature_map = feature_map,
    taxa_filter_stats = taxa$filter_stats
  )
}

empty_coda_signature <- function() {
  data.frame(
    analysis = character(0),
    sample_strata_col = character(0),
    sample_stratum = character(0),
    taxon = character(0),
    label = character(0),
    log_contrast_coefficient = numeric(0),
    risk_direction = character(0),
    abs_coefficient = numeric(0),
    stringsAsFactors = FALSE
  )
}

empty_coda_risk_scores <- function(covariates = character(0)) {
  base <- data.frame(
    PatientID = character(0),
    sample_strata_col = character(0),
    sample_stratum = character(0),
    .survival_time = numeric(0),
    .survival_status = numeric(0),
    stringsAsFactors = FALSE
  )
  for (covariate in covariates) {
    base[[covariate]] <- character(0)
  }
  base$risk_score <- numeric(0)
  base$risk_group_median <- character(0)
  base
}

empty_coda_metrics <- function() {
  data.frame(
    sample_strata_col = character(0),
    sample_stratum = character(0),
    status = character(0),
    reason = character(0),
    n_total = integer(0),
    n_with_stratum = integer(0),
    n_used = integer(0),
    events_used = integer(0),
    dropped_missing = integer(0),
    n_taxa_retained = integer(0),
    n_taxa_selected = integer(0),
    lambda = character(0),
    alpha = numeric(0),
    nfolds = integer(0),
    apparent_cindex = numeric(0),
    mean_cv_cindex = numeric(0),
    sd_cv_cindex = numeric(0),
    zero_handling = character(0),
    stringsAsFactors = FALSE
  )
}

coda_scalar_character <- function(value, default = NA_character_) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) {
    return(default)
  }
  as.character(value[[1]])
}

coda_scalar_numeric <- function(value, default = NA_real_) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) {
    return(default)
  }
  as.numeric(value[[1]])
}

coda_metric_row <- function(input,
                            status,
                            reason = NA_character_,
                            n_used = NA_integer_,
                            events_used = NA_integer_,
                            dropped_missing = NA_integer_,
                            n_taxa_selected = NA_integer_,
                            apparent_cindex = NA_real_,
                            mean_cv_cindex = NA_real_,
                            sd_cv_cindex = NA_real_,
                            options = list()) {
  data.frame(
    sample_strata_col = input$sample_strata_col,
    sample_stratum = input$sample_stratum,
    status = status,
    reason = reason,
    n_total = as.integer(input$n_total),
    n_with_stratum = as.integer(input$n_with_stratum),
    n_used = as.integer(n_used),
    events_used = as.integer(events_used),
    dropped_missing = as.integer(dropped_missing),
    n_taxa_retained = as.integer(ncol(input$taxa_matrix)),
    n_taxa_selected = as.integer(n_taxa_selected),
    lambda = coda_scalar_character(options$lambda),
    alpha = coda_scalar_numeric(options$alpha),
    nfolds = as.integer(coda_scalar_numeric(options$nfolds)),
    apparent_cindex = as.numeric(apparent_cindex),
    mean_cv_cindex = as.numeric(mean_cv_cindex),
    sd_cv_cindex = as.numeric(sd_cv_cindex),
    zero_handling = "coda4microbiome::impute_zeros",
    stringsAsFactors = FALSE
  )
}

build_coda_taxa_input <- function(patient_physeq,
                                  tax_agg_level = "Genus",
                                  min_prevalence = 0.10,
                                  min_mean_relative_abundance = 1e-5,
                                  sample_strata_col = NA_character_,
                                  sample_stratum = NA_character_) {
  rel_physeq <- relative_abundance_physeq(patient_physeq, tax_agg_level)
  filter_info <- filter_taxa_by_relative_abundance(
    rel_physeq,
    min_prevalence = min_prevalence,
    min_mean_relative_abundance = min_mean_relative_abundance
  )
  filter_stats <- filter_info$stats
  filter_stats <- cbind(
    data.frame(
      sample_strata_col = rep(as.character(sample_strata_col)[[1]], nrow(filter_stats)),
      sample_stratum = rep(as.character(sample_stratum)[[1]], nrow(filter_stats)),
      stringsAsFactors = FALSE
    ),
    filter_stats
  )

  if (length(filter_info$taxa) == 0) {
    empty <- data.frame(row.names = phyloseq::sample_names(patient_physeq))
    return(list(matrix = empty, filter_stats = filter_stats))
  }

  model_physeq <- patient_physeq |>
    tax_glom_rename(tax_agg_level) |>
    prune_empty_physeq()
  model_mat <- otu_samples_by_taxa(model_physeq)
  retained_taxa <- intersect(filter_info$taxa, colnames(model_mat))
  taxa_df <- as.data.frame(model_mat[, retained_taxa, drop = FALSE], stringsAsFactors = FALSE)
  list(matrix = taxa_df, filter_stats = filter_stats)
}

build_coda_model_inputs <- function(sample_physeq,
                                    spec,
                                    survival_config) {
  covariates <- as.character(unlist(spec$covariates %||% character(0), use.names = FALSE))
  metadata_df <- as.data.frame(phyloseq::sample_data(sample_physeq), stringsAsFactors = FALSE)
  sample_summary <- build_patient_sample_summary(
    metadata_df,
    patient_id_col = survival_config$patient_id_col,
    sample_strata_col = spec$sample_strata_col
  )
  base_df <- make_survival_base_df_from_metadata(
    sample_summary$metadata,
    time_col = survival_config$time_col,
    status_col = survival_config$status_col,
    patient_id_col = survival_config$patient_id_col,
    event_code = survival_config$event_code,
    censor_code = survival_config$censor_code,
    covariates = covariates
  )
  rownames(base_df) <- base_df$PatientID

  make_input <- function(patient_physeq, sample_strata_col, sample_stratum, n_with_stratum) {
    taxa <- build_coda_taxa_input(
      patient_physeq,
      tax_agg_level = spec$tax_agg_level,
      min_prevalence = spec$taxa_min_prevalence,
      min_mean_relative_abundance = spec$taxa_min_mean_relative_abundance,
      sample_strata_col = sample_strata_col,
      sample_stratum = sample_stratum
    )
    list(
      analysis = as.character(spec$name),
      base_df = base_df,
      taxa_matrix = taxa$matrix,
      filter_stats = taxa$filter_stats,
      sample_strata_col = sample_strata_col,
      sample_stratum = sample_stratum,
      n_total = nrow(base_df),
      n_with_stratum = n_with_stratum
    )
  }

  if (!survival_strata_enabled(spec$sample_strata_col)) {
    patient_physeq <- collapse_physeq_by_patient(sample_physeq, patient_id_col = survival_config$patient_id_col)
    input <- make_input(
      patient_physeq,
      sample_strata_col = NA_character_,
      sample_stratum = NA_character_,
      n_with_stratum = phyloseq::nsamples(patient_physeq)
    )
    return(list(inputs = list(input), filter_stats = input$filter_stats))
  }

  strata_col <- as.character(spec$sample_strata_col[[1]])
  fail_missing_columns(names(metadata_df), c(survival_config$patient_id_col, strata_col), "phyloseq sample_data")
  strata_values <- as.character(standardize_metadata_missing(metadata_df[[strata_col]]))
  strata_values[is_metadata_missing_like(strata_values)] <- NA_character_
  observed_strata <- sample_summary$observed_strata
  if (length(observed_strata) == 0) {
    return(list(inputs = list(), filter_stats = empty_taxa_filter_stats()))
  }

  inputs <- lapply(observed_strata, function(stratum) {
    keep <- !is.na(strata_values) & strata_values == stratum
    stratum_physeq <- phyloseq::prune_samples(keep, sample_physeq)
    patient_stratum_physeq <- collapse_physeq_by_patient(
      stratum_physeq,
      patient_id_col = survival_config$patient_id_col
    )
    make_input(
      patient_stratum_physeq,
      sample_strata_col = strata_col,
      sample_stratum = stratum,
      n_with_stratum = phyloseq::nsamples(patient_stratum_physeq)
    )
  })

  filter_stats <- do.call(rbind, lapply(inputs, `[[`, "filter_stats"))
  if (is.null(filter_stats)) filter_stats <- empty_taxa_filter_stats()
  list(inputs = inputs, filter_stats = filter_stats)
}

prepare_coda_covariates <- function(model_df, covariates = character(0)) {
  if (length(covariates) == 0) {
    return(NULL)
  }
  covar_df <- model_df[, 0, drop = FALSE]
  for (covariate in covariates) {
    value <- standardize_metadata_missing(model_df[[covariate]])
    if (any(is_metadata_missing_like(value))) {
      stop("internal error: missing covariate after complete-case filtering", call. = FALSE)
    }
    numeric_value <- suppressWarnings(as.numeric(as.character(value)))
    numeric_invalid <- is.na(numeric_value) | !is.finite(numeric_value)
    if (!any(numeric_invalid)) {
      if (stats::sd(numeric_value) == 0) {
        stop("covariate has fewer than two usable values: ", covariate, call. = FALSE)
      }
      covar_df[[covariate]] <- numeric_value
    } else {
      factor_value <- factor(as.character(value))
      if (length(levels(factor_value)) < 2) {
        stop("covariate has fewer than two usable levels: ", covariate, call. = FALSE)
      }
      covar_df[[covariate]] <- factor_value
    }
  }
  covar_df
}

resolve_coda_coxnet <- function(coda_coxnet_fn = NULL) {
  if (!is.null(coda_coxnet_fn)) {
    return(coda_coxnet_fn)
  }
  if (!requireNamespace("coda4microbiome", quietly = TRUE)) {
    stop("The R package 'coda4microbiome' is required for coda survival analysis.", call. = FALSE)
  }
  getExportedValue("coda4microbiome", "coda_coxnet")
}

coda_result_value <- function(result, name) {
  if (!name %in% names(result)) {
    return(NULL)
  }
  result[[name]]
}

coda_signature_table <- function(result, input) {
  coefficients <- coda_result_value(result, "log-contrast coefficients")
  taxa <- coda_result_value(result, "taxa.name")
  if (is.null(coefficients) || length(coefficients) == 0) {
    return(empty_coda_signature())
  }
  coefficients <- as.numeric(coefficients)
  if (is.null(taxa) || length(taxa) != length(coefficients)) {
    taxa <- names(coefficients)
  }
  if (is.null(taxa) || length(taxa) != length(coefficients)) {
    taxa <- paste0("taxon_", seq_along(coefficients))
  }
  taxa <- as.character(taxa)
  label <- taxa
  if (!is.na(input$sample_strata_col) && nzchar(input$sample_strata_col)) {
    label <- paste0(input$sample_strata_col, "=", input$sample_stratum, " | ", taxa)
  }
  signature <- data.frame(
    analysis = input$analysis,
    sample_strata_col = input$sample_strata_col,
    sample_stratum = input$sample_stratum,
    taxon = taxa,
    label = label,
    log_contrast_coefficient = coefficients,
    risk_direction = ifelse(coefficients >= 0, "positive", "negative"),
    abs_coefficient = abs(coefficients),
    stringsAsFactors = FALSE
  )
  signature[order(-signature$abs_coefficient, signature$taxon), , drop = FALSE]
}

coda_risk_score_table <- function(result, input, model_df, covariates = character(0)) {
  risk_score <- coda_result_value(result, "risk.score")
  if (is.null(risk_score) || length(risk_score) == 0) {
    return(empty_coda_risk_scores(covariates))
  }
  risk_score <- as.numeric(risk_score)
  if (length(risk_score) != nrow(model_df)) {
    stop(
      "coda4microbiome returned ",
      length(risk_score),
      " risk score(s) for ",
      nrow(model_df),
      " modeled patient(s).",
      call. = FALSE
    )
  }
  out <- model_df[, unique(c("PatientID", ".survival_time", ".survival_status", covariates)), drop = FALSE]
  out <- cbind(
    out[, "PatientID", drop = FALSE],
    data.frame(
      sample_strata_col = input$sample_strata_col,
      sample_stratum = input$sample_stratum,
      stringsAsFactors = FALSE
    ),
    out[, setdiff(names(out), "PatientID"), drop = FALSE]
  )
  median_score <- stats::median(risk_score, na.rm = TRUE)
  out$risk_score <- risk_score
  out$risk_group_median <- ifelse(risk_score <= median_score, "below_or_equal_median", "above_median")
  out
}

fit_coda_model_input <- function(input,
                                 covariates = character(0),
                                 min_n = 30,
                                 min_events = 10,
                                 options = list(),
                                 coda_coxnet_fn = NULL) {
  taxa_aligned <- align_taxa_matrix_to_patients(input$taxa_matrix, rownames(input$base_df))
  model_df <- cbind(input$base_df, taxa_aligned)
  taxa_cols <- colnames(taxa_aligned)
  required <- unique(c(".survival_time", ".survival_status", covariates, taxa_cols))
  model_df <- standardize_metadata_missing_df(model_df, intersect(required, names(model_df)))
  complete <- stats::complete.cases(model_df[, required, drop = FALSE])
  n_used <- sum(complete)
  dropped_missing <- nrow(model_df) - n_used
  if (n_used < as.integer(min_n)) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", "too_few_complete_cases", n_used, NA, dropped_missing, 0, options = options)
    ))
  }

  model_df <- model_df[complete, , drop = FALSE]
  model_df$.survival_time <- as.numeric(model_df$.survival_time)
  model_df$.survival_status <- as.numeric(model_df$.survival_status)
  events_used <- sum(model_df$.survival_status == 1)
  if (events_used < as.integer(min_events)) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", "too_few_events", n_used, events_used, dropped_missing, 0, options = options)
    ))
  }
  if (length(taxa_cols) < 2) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", "too_few_taxa", n_used, events_used, dropped_missing, 0, options = options)
    ))
  }
  if (n_used < as.integer(options$nfolds)) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", "too_few_cv_folds", n_used, events_used, dropped_missing, 0, options = options)
    ))
  }

  covar_df <- tryCatch(
    prepare_coda_covariates(model_df, covariates = covariates),
    error = function(e) e
  )
  if (inherits(covar_df, "error")) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", conditionMessage(covar_df), n_used, events_used, dropped_missing, 0, options = options)
    ))
  }
  if (inherits(coda_coxnet_fn, "error")) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(input, "skipped", conditionMessage(coda_coxnet_fn), n_used, events_used, dropped_missing, 0, options = options)
    ))
  }

  taxa_mat <- as.data.frame(model_df[, taxa_cols, drop = FALSE], stringsAsFactors = FALSE)
  rownames(taxa_mat) <- model_df$PatientID
  if (!is.na(options$seed)) {
    set.seed(as.integer(options$seed))
  }
  coda_fit <- tryCatch(
    coda_coxnet_fn(
      x = taxa_mat,
      time = model_df$.survival_time,
      status = model_df$.survival_status,
      covar = covar_df,
      lambda = options$lambda,
      nvar = options$nvar,
      alpha = options$alpha,
      nfolds = options$nfolds,
      showPlots = options$show_plots,
      coef_threshold = options$coef_threshold
    ),
    error = function(e) e
  )
  if (inherits(coda_fit, "error")) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(
        input,
        "failed",
        paste0("coda_coxnet failed: ", conditionMessage(coda_fit)),
        n_used,
        events_used,
        dropped_missing,
        0,
        options = options
      )
    ))
  }

  signature <- coda_signature_table(coda_fit, input)
  risk_scores <- tryCatch(
    coda_risk_score_table(coda_fit, input, model_df, covariates = covariates),
    error = function(e) e
  )
  if (inherits(risk_scores, "error")) {
    return(list(
      signature = signature,
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = coda_metric_row(
        input,
        "failed",
        paste0("risk score output failed: ", conditionMessage(risk_scores)),
        n_used,
        events_used,
        dropped_missing,
        nrow(signature),
        options = options
      )
    ))
  }

  list(
    signature = signature,
    risk_scores = risk_scores,
    metrics = coda_metric_row(
      input,
      "fit",
      NA_character_,
      n_used,
      events_used,
      dropped_missing,
      nrow(signature),
      apparent_cindex = coda_scalar_numeric(coda_result_value(coda_fit, "apparent Cindex")),
      mean_cv_cindex = coda_scalar_numeric(coda_result_value(coda_fit, "mean cv-Cindex")),
      sd_cv_cindex = coda_scalar_numeric(coda_result_value(coda_fit, "sd cv-Cindex")),
      options = options
    )
  )
}

run_survival_coda_models <- function(sample_physeq,
                                     spec,
                                     survival_config,
                                     covariates = character(0),
                                     min_n = 30,
                                     min_events = 10,
                                     coda_coxnet_fn = NULL) {
  coda_inputs <- build_coda_model_inputs(sample_physeq, spec, survival_config)
  if (length(coda_inputs$inputs) == 0) {
    return(list(
      signature = empty_coda_signature(),
      risk_scores = empty_coda_risk_scores(covariates),
      metrics = empty_coda_metrics(),
      taxa_filter_stats = coda_inputs$filter_stats
    ))
  }

  resolved_coda <- tryCatch(
    resolve_coda_coxnet(coda_coxnet_fn),
    error = function(e) e
  )
  if (inherits(resolved_coda, "error")) {
    warning(conditionMessage(resolved_coda), call. = FALSE)
  }

  fits <- lapply(coda_inputs$inputs, function(input) {
    fit_coda_model_input(
      input,
      covariates = covariates,
      min_n = min_n,
      min_events = min_events,
      options = spec$coda4microbiome,
      coda_coxnet_fn = resolved_coda
    )
  })

  signature <- do.call(rbind, lapply(fits, `[[`, "signature"))
  risk_scores <- do.call(rbind, lapply(fits, `[[`, "risk_scores"))
  metrics <- do.call(rbind, lapply(fits, `[[`, "metrics"))
  if (is.null(signature)) signature <- empty_coda_signature()
  if (is.null(risk_scores)) risk_scores <- empty_coda_risk_scores(covariates)
  if (is.null(metrics)) metrics <- empty_coda_metrics()

  list(
    signature = signature,
    risk_scores = risk_scores,
    metrics = metrics,
    taxa_filter_stats = coda_inputs$filter_stats
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
                            sample_strata_col = NA_character_,
                            sample_stratum = NA_character_,
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
      skip = survival_skip_row(feature, feature_family, feature_label, sample_strata_col, sample_stratum, "too_few_complete_cases", n_total, n_used, NA, dropped_missing)
    ))
  }

  model_df <- model_df[complete, , drop = FALSE]
  model_df$.survival_time <- as.numeric(model_df$.survival_time)
  model_df$.survival_status <- as.numeric(model_df$.survival_status)
  events_used <- sum(model_df$.survival_status == 1)
  if (events_used < as.integer(min_events)) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, sample_strata_col, sample_stratum, "too_few_events", n_total, n_used, events_used, dropped_missing)
    ))
  }

  coerced <- tryCatch(
    coerce_survival_model_columns(model_df, covariates = covariates),
    error = function(e) e
  )
  if (inherits(coerced, "error")) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, sample_strata_col, sample_stratum, conditionMessage(coerced), n_total, n_used, events_used, dropped_missing)
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
      skip = survival_skip_row(feature, feature_family, feature_label, sample_strata_col, sample_stratum, paste0("coxph failed: ", conditionMessage(fit)), n_total, n_used, events_used, dropped_missing)
    ))
  }

  fit_summary <- summary(fit)
  if (!".feature_z" %in% rownames(fit_summary$coefficients)) {
    return(list(
      result = NULL,
      skip = survival_skip_row(feature, feature_family, feature_label, sample_strata_col, sample_stratum, "coxph omitted feature coefficient", n_total, n_used, events_used, dropped_missing)
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
    sample_strata_col = sample_strata_col,
    sample_stratum = sample_stratum,
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
                              sample_strata_col,
                              sample_stratum,
                              reason,
                              n_total,
                              n_used,
                              events_used,
                              dropped_missing) {
  data.frame(
    feature = feature,
    feature_family = feature_family,
    label = feature_label,
    sample_strata_col = sample_strata_col,
    sample_stratum = sample_stratum,
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
      sample_strata_col = if ("sample_strata_col" %in% names(feature_map)) feature_map$sample_strata_col[[i]] else NA_character_,
      sample_stratum = if ("sample_stratum" %in% names(feature_map)) feature_map$sample_stratum[[i]] else NA_character_,
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
  keep_cols <- c(
    "feature", "feature_family", "label", "sample_strata_col", "sample_stratum",
    "n_total", "n_used", "events_used", "dropped_missing"
  )
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
