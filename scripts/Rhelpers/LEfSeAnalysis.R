source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source_if_needed(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))

is_missing_lefse_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

normalize_lefse_config <- function(config = list(),
                                   project_config = list()) {
  config <- config %||% list()

  if (is.null(config$source_comparisons)) config$source_comparisons <- "differential_abundance"
  if (is.null(config$abundance_scale)) config$abundance_scale <- "project_norm_to_relative_abundance"
  if (is.null(config$p_adjust_method)) config$p_adjust_method <- "BH"
  if (is.null(config$kruskal_threshold)) config$kruskal_threshold <- 0.05
  if (is.null(config$wilcox_threshold)) config$wilcox_threshold <- 0.05
  if (is.null(config$lda_threshold)) config$lda_threshold <- 2
  if (is.null(config$filter)) config$filter <- "nonzero"
  if (is.null(config$lda_safe_filter)) config$lda_safe_filter <- TRUE
  if (is.null(config$relative_abundance_scale)) config$relative_abundance_scale <- 1e6
  if (is.null(config$subclass_col)) config$subclass_col <- NULL
  if (is.null(config$seed)) config$seed <- 1234
  if (is.null(config$tax_agg_level)) config$tax_agg_level <- project_config$tax_agg_level %||% "Genus"
  if (is.null(config$tax_label_level)) config$tax_label_level <- config$tax_agg_level

  config
}

merge_lefse_overrides <- function(base_specs, override_specs) {
  if (is.null(override_specs) || length(override_specs) == 0) {
    return(base_specs)
  }

  for (override in override_specs) {
    if (is.null(override$name)) {
      stop("Each LEfSe comparison override must define name.", call. = FALSE)
    }
    match_idx <- vapply(
      base_specs,
      function(spec) identical(as.character(spec$name), as.character(override$name)),
      logical(1)
    )
    if (!any(match_idx)) {
      stop("LEfSe comparison override does not match a DA comparison: ", override$name, call. = FALSE)
    }
    idx <- which(match_idx)[[1]]
    base_specs[[idx]] <- utils::modifyList(base_specs[[idx]], override)
  }

  base_specs
}

lefse_comparison_specs <- function(lefse_config, da_config = NULL) {
  source_comparisons <- as.character(lefse_config$source_comparisons %||% "differential_abundance")

  if (identical(source_comparisons, "differential_abundance")) {
    if (is.null(da_config) || is.null(da_config$comparisons) || length(da_config$comparisons) == 0) {
      stop("LEfSe source_comparisons is differential_abundance, but no DA comparisons are available.", call. = FALSE)
    }
    return(merge_lefse_overrides(da_config$comparisons, lefse_config$comparisons))
  }

  if (is.null(lefse_config$comparisons) || length(lefse_config$comparisons) == 0) {
    stop("LEfSe config must define comparisons when source_comparisons is not differential_abundance.", call. = FALSE)
  }
  lefse_config$comparisons
}

resolve_lefse_norm_method <- function(lefse_config, project_config) {
  config_value(lefse_config, "norm_method") %||%
    config_value(project_config, "norm_method") %||%
    "noNorm"
}

resolve_lefse_pseudocount <- function(lefse_config, project_config) {
  config_value(lefse_config, "pseudocount") %||%
    config_value(project_config, "pseudocount") %||%
    1
}

comparison_class_col <- function(spec) {
  spec$class_col %||% spec$group
}

apply_lefse_class_mapping <- function(metadata_df, spec) {
  class_col <- comparison_class_col(spec)
  fail_missing_columns(names(metadata_df), class_col, "phyloseq sample_data")

  class_map <- spec$class_map %||% spec$class_mapping
  if (is.null(class_map)) {
    return(list(
      metadata = metadata_df,
      class_col = class_col,
      class_levels = spec$class_levels
    ))
  }

  if (is.null(names(class_map)) ||
      length(class_map) != 2 ||
      any(!nzchar(names(class_map)))) {
    stop("LEfSe class_map for comparison '", spec$name, "' must define exactly two named classes.", call. = FALSE)
  }

  mapped <- rep(NA_character_, nrow(metadata_df))
  source_values <- as.character(metadata_df[[class_col]])
  for (target in names(class_map)) {
    mapped[source_values %in% as.character(unlist(class_map[[target]], use.names = FALSE))] <- target
  }

  keep <- !is.na(mapped)
  if (!any(keep)) {
    stop("LEfSe class_map for comparison '", spec$name, "' selected zero samples.", call. = FALSE)
  }

  metadata_df <- metadata_df[keep, , drop = FALSE]
  metadata_df$.LEfSeClass <- factor(mapped[keep], levels = names(class_map))

  list(
    metadata = metadata_df,
    class_col = ".LEfSeClass",
    class_levels = names(class_map)
  )
}

validate_lefse_binary_class <- function(metadata_df, spec) {
  mapped <- apply_lefse_class_mapping(metadata_df, spec)
  metadata_df <- mapped$metadata
  class_col <- mapped$class_col
  configured_levels <- mapped$class_levels

  values <- as.character(metadata_df[[class_col]])
  keep <- !is.na(values) & nzchar(trimws(values))
  metadata_df <- metadata_df[keep, , drop = FALSE]
  values <- as.character(metadata_df[[class_col]])

  if (!is_missing_lefse_value(configured_levels)) {
    configured_levels <- as.character(unlist(configured_levels, use.names = FALSE))
    if (length(configured_levels) != 2) {
      stop("LEfSe class_levels for comparison '", spec$name, "' must contain exactly two levels.", call. = FALSE)
    }
    keep <- values %in% configured_levels
    metadata_df <- metadata_df[keep, , drop = FALSE]
    values <- as.character(metadata_df[[class_col]])
    missing_levels <- setdiff(configured_levels, unique(values))
    if (length(missing_levels) > 0) {
      stop(
        "LEfSe comparison '",
        spec$name,
        "' is missing configured class level(s): ",
        paste(missing_levels, collapse = ", "),
        call. = FALSE
      )
    }
    class_levels <- configured_levels
  } else {
    class_levels <- unique(values)
    if (is.factor(metadata_df[[class_col]])) {
      class_levels <- levels(metadata_df[[class_col]])[levels(metadata_df[[class_col]]) %in% class_levels]
    }
    if (length(class_levels) != 2) {
      stop(
        "LEfSe comparison '",
        spec$name,
        "' must have exactly two observed class levels in ",
        class_col,
        "; found: ",
        paste(class_levels, collapse = ", "),
        ". Configure class_levels or class_map for binary LEfSe analysis.",
        call. = FALSE
      )
    }
  }

  metadata_df[[class_col]] <- factor(as.character(metadata_df[[class_col]]), levels = class_levels)

  list(
    metadata = metadata_df,
    class_col = class_col,
    class_levels = class_levels
  )
}

filter_lefse_taxa <- function(physeq, filter = "nonzero") {
  filter <- as.character(filter %||% "nonzero")
  if (!identical(filter, "nonzero")) {
    stop("Unsupported LEfSe filter: ", filter, ". Supported filter: nonzero.", call. = FALSE)
  }

  phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
}

filter_lefse_lda_safe_taxa <- function(physeq,
                                       class_col,
                                       class_levels = NULL,
                                       min_variance = 0) {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), class_col, "phyloseq sample_data")

  class_values <- as.character(metadata_df[[class_col]])
  if (!is_missing_lefse_value(class_levels)) {
    class_levels <- as.character(unlist(class_levels, use.names = FALSE))
  } else {
    class_levels <- unique(class_values[!is.na(class_values) & nzchar(trimws(class_values))])
  }
  if (length(class_levels) != 2) {
    stop("LEfSe LDA-safe filter requires exactly two class levels.", call. = FALSE)
  }

  class_factor <- factor(class_values, levels = class_levels)
  otu_mat <- otu_samples_by_taxa(physeq)
  keep_taxa <- vapply(seq_len(ncol(otu_mat)), function(i) {
    feature_values <- otu_mat[, i]
    all(vapply(class_levels, function(level) {
      group_values <- feature_values[class_factor == level]
      group_values <- group_values[is.finite(group_values)]
      length(group_values) > 1 &&
        is.finite(stats::var(group_values)) &&
        stats::var(group_values) > min_variance
    }, logical(1)))
  }, logical(1))

  if (!any(keep_taxa)) {
    stop(
      "LEfSe LDA-safe filter removed all taxa because no taxa had nonzero within-class variance in both classes.",
      call. = FALSE
    )
  }

  list(
    physeq = phyloseq::prune_taxa(colnames(otu_mat)[keep_taxa], physeq),
    before_taxa = ncol(otu_mat),
    retained_taxa = sum(keep_taxa),
    removed_taxa = sum(!keep_taxa)
  )
}

validate_lefse_abundance_matrix <- function(otu_mat,
                                            comparison_name,
                                            abundance_scale,
                                            norm_method) {
  if (any(!is.finite(otu_mat))) {
    stop(
      "LEfSe comparison '",
      comparison_name,
      "' has non-finite values after normalization. Switch LEfSe to raw-count relative abundance.",
      call. = FALSE
    )
  }
  if (any(otu_mat < 0, na.rm = TRUE)) {
    stop(
      "LEfSe comparison '",
      comparison_name,
      "' has negative values after ",
      norm_method,
      " normalization under abundance_scale=",
      abundance_scale,
      ". Switch LEfSe to raw-count relative abundance.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

close_lefse_relative_abundance <- function(physeq,
                                           comparison_name,
                                           relative_abundance_scale = 1e6) {
  if (!is.finite(relative_abundance_scale) || relative_abundance_scale <= 0) {
    stop("LEfSe relative_abundance_scale must be positive and finite.", call. = FALSE)
  }

  otu_mat <- otu_samples_by_taxa(physeq)
  row_totals <- rowSums(otu_mat)
  keep <- is.finite(row_totals) & row_totals > 0
  if (!any(keep)) {
    stop("LEfSe comparison '", comparison_name, "' has zero samples with positive abundance.", call. = FALSE)
  }
  if (any(!keep)) {
    physeq <- phyloseq::prune_samples(keep, physeq)
    otu_mat <- otu_samples_by_taxa(physeq)
    row_totals <- rowSums(otu_mat)
  }

  rel_mat <- sweep(otu_mat, 1, row_totals, FUN = "/")
  rel_mat <- rel_mat * relative_abundance_scale
  set_otu_samples_by_taxa(physeq, rel_mat)
}

physeq_to_relative_abundance <- function(physeq,
                                         comparison_name,
                                         abundance_scale,
                                         norm_method,
                                         pseudocount,
                                         relative_abundance_scale = 1e6) {
  if (!identical(abundance_scale, "project_norm_to_relative_abundance") &&
      !identical(abundance_scale, "raw_relative_abundance")) {
    stop("Unsupported LEfSe abundance_scale: ", abundance_scale, call. = FALSE)
  }

  if (identical(abundance_scale, "project_norm_to_relative_abundance")) {
    physeq <- counts_normalization(
      physeq = physeq,
      norm_method = norm_method,
      pseudocount = pseudocount
    )
  }

  otu_mat <- otu_samples_by_taxa(physeq)
  validate_lefse_abundance_matrix(
    otu_mat,
    comparison_name = comparison_name,
    abundance_scale = abundance_scale,
    norm_method = norm_method
  )

  close_lefse_relative_abundance(
    physeq,
    comparison_name = comparison_name,
    relative_abundance_scale = relative_abundance_scale
  )
}

prepare_lefse_physeq <- function(physeq,
                                 spec,
                                 lefse_config,
                                 project_config = list()) {
  physeq <- apply_sample_filter(physeq, spec$sample_filter)
  physeq <- apply_factor_levels(physeq, spec$factor_levels)

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  class_info <- validate_lefse_binary_class(metadata_df, spec)
  physeq <- phyloseq::prune_samples(rownames(class_info$metadata), physeq)
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(class_info$metadata)

  tax_agg_level <- spec$tax_agg_level %||% lefse_config$tax_agg_level
  physeq <- tax_glom_rename(physeq, tax_agg_level)
  physeq <- prune_empty_physeq(physeq)

  abundance_scale <- as.character(lefse_config$abundance_scale)
  norm_method <- resolve_lefse_norm_method(lefse_config, project_config)
  pseudocount <- resolve_lefse_pseudocount(lefse_config, project_config)
  relative_abundance_scale <- as.numeric(lefse_config$relative_abundance_scale)

  physeq <- physeq_to_relative_abundance(
    physeq,
    comparison_name = spec$name,
    abundance_scale = abundance_scale,
    norm_method = norm_method,
    pseudocount = pseudocount,
    relative_abundance_scale = relative_abundance_scale
  )
  physeq <- filter_lefse_taxa(physeq, lefse_config$filter)
  physeq <- prune_empty_physeq(physeq)
  lda_filter <- list(
    enabled = FALSE,
    before_taxa = ntaxa(physeq),
    retained_taxa = ntaxa(physeq),
    removed_taxa = 0
  )
  if (truthy_flag(lefse_config$lda_safe_filter, default = TRUE)) {
    lda_filter <- filter_lefse_lda_safe_taxa(
      physeq,
      class_col = class_info$class_col,
      class_levels = class_info$class_levels
    )
    physeq <- lda_filter$physeq
    lda_filter$enabled <- TRUE
  }
  physeq <- close_lefse_relative_abundance(
    physeq,
    comparison_name = spec$name,
    relative_abundance_scale = relative_abundance_scale
  )

  list(
    physeq = physeq,
    class_col = class_info$class_col,
    class_levels = class_info$class_levels,
    norm_method = norm_method,
    pseudocount = pseudocount,
    abundance_scale = abundance_scale,
    relative_abundance_scale = relative_abundance_scale,
    lda_filter = lda_filter
  )
}

physeq_to_lefse_se <- function(physeq) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The R package 'SummarizedExperiment' is required for LEfSe analysis.", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("The R package 'S4Vectors' is required for LEfSe analysis.", call. = FALSE)
  }

  assay_mat <- t(otu_samples_by_taxa(physeq))
  tax_df <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(relative_abundance = assay_mat),
    colData = S4Vectors::DataFrame(metadata_df),
    rowData = S4Vectors::DataFrame(tax_df)
  )
}

validate_relative_abundance_closure <- function(physeq,
                                                target_sum = 1e6,
                                                tolerance = 1e-8) {
  otu_mat <- otu_samples_by_taxa(physeq)
  row_sums <- rowSums(otu_mat)
  allowed_error <- tolerance * max(1, abs(target_sum))
  all(is.finite(row_sums)) && all(abs(row_sums - target_sum) <= allowed_error)
}

format_lefse_results <- function(results_df, class_levels) {
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(data.frame(
      taxon = character(0),
      lda_score = numeric(0),
      lda_abs_score = numeric(0),
      enriched_class = character(0),
      rank = integer(0)
    ))
  }

  colnames(results_df)[colnames(results_df) == "features"] <- "taxon"
  colnames(results_df)[colnames(results_df) == "scores"] <- "lda_score"
  results_df$lda_score <- as.numeric(results_df$lda_score)
  results_df$lda_abs_score <- abs(results_df$lda_score)
  results_df$enriched_class <- ifelse(results_df$lda_score >= 0, class_levels[[2]], class_levels[[1]])
  results_df <- results_df[order(-results_df$lda_abs_score, results_df$taxon), , drop = FALSE]
  results_df$rank <- seq_len(nrow(results_df))
  results_df[, c("taxon", "lda_score", "lda_abs_score", "enriched_class", "rank")]
}

write_lefse_results <- function(results_df, out_dir, trial_id, comparison_name) {
  comparison_dir <- file.path(out_dir, "LEfSe", comparison_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(
    comparison_dir,
    paste0(trial_id, "_", comparison_name, "_LEfSeResults.tsv")
  )
  write.table(
    results_df,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  out_file
}

write_lefse_run_note <- function(out_dir,
                                 trial_id,
                                 comparison_name,
                                 note_lines) {
  comparison_dir <- file.path(out_dir, "LEfSe", comparison_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  note_file <- file.path(
    comparison_dir,
    paste0(trial_id, "_", comparison_name, "_LEfSeRunNote.txt")
  )
  writeLines(note_lines, note_file)
  note_file
}
