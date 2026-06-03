source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

read_da_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The R package 'yaml' is required to read analysis config files.", call. = FALSE)
  }

  config <- yaml::read_yaml(path)
  if ("differential_abundance" %in% names(config)) {
    config <- config$differential_abundance
  }
  normalize_da_config(config)
}

legacy_comparison_spec <- function(name,
                                   tax_agg_level = "Genus",
                                   tax_label_level = "Genus",
                                   alpha = 0.05,
                                   lfc_cutoff = 1) {
  common <- list(
    name = name,
    tax_agg_level = tax_agg_level,
    tax_label_level = tax_label_level,
    alpha = alpha,
    lfc_cutoff = lfc_cutoff
  )

  if (name %in% c("CellLineControltoTumor", "CellLineControltoNontumor", "NegativeControl", "AllControl")) {
    sample_values <- switch(
      name,
      CellLineControltoTumor = c("CellLineControl", "Tumor"),
      CellLineControltoNontumor = c("CellLineControl", "Nontumor", "NormalTissue"),
      NegativeControl = c("NegativeControl", "Tumor", "Nontumor", "NormalTissue"),
      AllControl = NULL
    )
    common$sample_filter <- if (is.null(sample_values)) {
      list(SampleType = "*")
    } else {
      list(SampleType = sample_values)
    }
    common$factor_levels <- list(ControlStatus = c("Control", "PatientSample"))
    common$formula <- "ControlStatus"
    common$group <- "ControlStatus"
    common$coefficient <- "ControlStatusPatientSample"
    common$structural_zero_groups <- c(
      "structural_zero (ControlStatus = Control)",
      "structural_zero (ControlStatus = PatientSample)"
    )
    return(common)
  }

  if (name == "PatientSample") {
    common$sample_filter <- list(SampleType = c("Tumor", "Nontumor", "NormalTissue"))
    common$factor_levels <- list(SampleType = c("Nontumor", "NormalTissue", "Tumor"))
    common$formula <- "SampleType + PatientID"
    common$group <- "SampleType"
    common$coefficient <- "SampleTypeTumor"
    common$structural_zero_groups <- c(
      "structural_zero (SampleType = Nontumor)",
      "structural_zero (SampleType = Tumor)"
    )
    return(common)
  }

  if (name == "TumorType") {
    common$sample_filter <- list(SampleType = "Tumor")
    common$factor_levels <- list(TumorType = c("HCC", "iCC"))
    common$formula <- "TumorType"
    common$group <- "TumorType"
    common$coefficient <- "TumorTypeiCC"
    common$structural_zero_groups <- c(
      "structural_zero (TumorType = HCC)",
      "structural_zero (TumorType = iCC)"
    )
    return(common)
  }

  stop("Unknown legacy DA comparison: ", name, call. = FALSE)
}

build_legacy_da_config <- function(trial_id,
                                   comparisons,
                                   methods = c("ANCOMBC"),
                                   out_dir,
                                   norm_method = "noNorm",
                                   pseudocount = 1,
                                   tax_agg_level = "Genus",
                                   tax_label_level = "Genus",
                                   alpha = 0.05,
                                   lfc_cutoff = 1) {
  list(
    trialID = trial_id,
    methods = methods,
    output_dir = out_dir,
    norm_method = norm_method,
    pseudocount = pseudocount,
    tax_agg_level = tax_agg_level,
    tax_label_level = tax_label_level,
    alpha = alpha,
    lfc_cutoff = lfc_cutoff,
    comparisons = lapply(
      comparisons,
      legacy_comparison_spec,
      tax_agg_level = tax_agg_level,
      tax_label_level = tax_label_level,
      alpha = alpha,
      lfc_cutoff = lfc_cutoff
    )
  )
}

normalize_da_config <- function(config) {
  if (is.null(config$methods)) config$methods <- c("ANCOMBC")
  if (is.null(config$norm_method)) config$norm_method <- "noNorm"
  if (is.null(config$pseudocount)) config$pseudocount <- 1
  if (is.null(config$tax_agg_level)) config$tax_agg_level <- "Genus"
  if (is.null(config$tax_label_level)) config$tax_label_level <- config$tax_agg_level
  if (is.null(config$alpha)) config$alpha <- 0.05
  if (is.null(config$lfc_cutoff)) config$lfc_cutoff <- 1

  if (is.null(config$comparisons) || length(config$comparisons) == 0) {
    stop("DA config must define at least one comparison.", call. = FALSE)
  }

  unsupported <- setdiff(config$methods, "ANCOMBC")
  if (length(unsupported) > 0) {
    stop(
      "Only ANCOMBC is supported by the config-driven DA interface. Unsupported method(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  config
}

apply_sample_filter <- function(physeq, sample_filter) {
  if (is.null(sample_filter) || length(sample_filter) == 0) {
    return(physeq)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  keep <- rep(TRUE, nrow(metadata_df))

  for (column in names(sample_filter)) {
    fail_missing_columns(names(metadata_df), column, "phyloseq sample_data")
    values <- unlist(sample_filter[[column]], use.names = FALSE)
    if (length(values) == 1 && values == "*") {
      keep <- keep & !is.na(metadata_df[[column]])
    } else {
      keep <- keep & as.character(metadata_df[[column]]) %in% as.character(values)
    }
  }

  if (!any(keep)) {
    stop("Sample filter selected zero samples.", call. = FALSE)
  }

  phyloseq::prune_samples(keep, physeq)
}

apply_factor_levels <- function(physeq, factor_levels) {
  if (is.null(factor_levels) || length(factor_levels) == 0) {
    return(physeq)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  for (column in names(factor_levels)) {
    levels <- as.character(unlist(factor_levels[[column]], use.names = FALSE))
    fail_missing_columns(names(metadata_df), column, "phyloseq sample_data")
    observed <- unique(as.character(metadata_df[[column]]))
    missing_levels <- setdiff(observed[!is.na(observed)], levels)
    if (length(missing_levels) > 0) {
      stop(
        "Observed value(s) in ",
        column,
        " are absent from configured factor levels: ",
        paste(missing_levels, collapse = ", "),
        call. = FALSE
      )
    }
    metadata_df[[column]] <- factor(metadata_df[[column]], levels = levels)
  }
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(metadata_df)
  physeq
}

validate_comparison_spec <- function(physeq, spec) {
  required <- c("name", "formula", "group", "coefficient", "structural_zero_groups")
  missing <- required[!vapply(required, function(x) x %in% names(spec) && !is.null(spec[[x]]), logical(1))]
  if (length(missing) > 0) {
    stop(
      "DA comparison spec is missing required field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  formula_terms <- all.vars(stats::as.formula(paste0("~ ", spec$formula)))
  fail_missing_columns(names(metadata_df), unique(c(spec$group, formula_terms)), "phyloseq sample_data")

  grouping_values <- unique(metadata_df[[spec$group]])
  grouping_values <- grouping_values[!is.na(grouping_values)]
  if (length(grouping_values) < 2) {
    stop(
      "Comparison '",
      spec$name,
      "' has fewer than two non-empty group levels in ",
      spec$group,
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

prepare_da_physeq <- function(physeq, spec, global_config, select_taxa = NULL) {
  physeq <- apply_sample_filter(physeq, spec$sample_filter)
  physeq <- apply_factor_levels(physeq, spec$factor_levels)
  validate_comparison_spec(physeq, spec)

  tax_agg_level <- if (!is.null(spec$tax_agg_level)) spec$tax_agg_level else global_config$tax_agg_level
  physeq <- tax_glom_rename(physeq, tax_agg_level)
  physeq <- prune_empty_physeq(physeq)

  if (!is.null(select_taxa)) {
    keep_taxa <- intersect(select_taxa, phyloseq::taxa_names(physeq))
    if (length(keep_taxa) == 0) {
      stop("Selected taxa list has no overlap with comparison taxa.", call. = FALSE)
    }
    physeq <- phyloseq::prune_taxa(keep_taxa, physeq)
  }

  physeq
}

extract_ancombc_table <- function(ancombc_output, table_name, coefficient) {
  table <- ancombc_output[["res"]][[table_name]]
  required <- c("taxon", coefficient)
  fail_missing_columns(colnames(table), required, paste0("ANCOMBC ", table_name, " result"))
  table[, required, drop = FALSE]
}

format_structural_zero <- function(ancombc_output, structural_zero_groups) {
  zero_ind <- ancombc_output[["zero_ind"]]
  if (is.null(zero_ind)) {
    return(data.frame(taxon = character(0), struc0 = character(0)))
  }

  required <- c("taxon", structural_zero_groups)
  fail_missing_columns(colnames(zero_ind), required, "ANCOMBC structural-zero result")
  struc0_df <- zero_ind[, required, drop = FALSE]
  colnames(struc0_df) <- c("taxon", "struc0_group1", "struc0_group2")

  struc0_df$struc0_group1 <- as.logical(struc0_df$struc0_group1)
  struc0_df$struc0_group2 <- as.logical(struc0_df$struc0_group2)
  struc0_df$struc0 <- ifelse(
    struc0_df$struc0_group1 & !struc0_df$struc0_group2,
    "group1",
    ifelse(!struc0_df$struc0_group1 & struc0_df$struc0_group2, "group2", NA)
  )
  struc0_df[, c("taxon", "struc0"), drop = FALSE]
}

standardize_ancombc_results <- function(ancombc_output, spec, alpha, lfc_cutoff) {
  lfc_df <- extract_ancombc_table(ancombc_output, "lfc", spec$coefficient)
  p_df <- extract_ancombc_table(ancombc_output, "p_val", spec$coefficient)
  q_df <- extract_ancombc_table(ancombc_output, "q_val", spec$coefficient)
  se_df <- extract_ancombc_table(ancombc_output, "se", spec$coefficient)

  colnames(lfc_df) <- c("taxon", "log2FoldChange")
  colnames(p_df) <- c("taxon", "p")
  colnames(q_df) <- c("taxon", "padj")
  colnames(se_df) <- c("taxon", "se")

  results <- Reduce(
    function(x, y) merge(x, y, by = "taxon", all.x = TRUE),
    list(lfc_df, p_df, q_df, se_df)
  )

  struc0_df <- format_structural_zero(ancombc_output, spec$structural_zero_groups)
  if (nrow(struc0_df) > 0) {
    results <- merge(results, struc0_df, by = "taxon", all.x = TRUE)
  } else {
    results$struc0 <- NA_character_
  }

  results$log2FoldChange <- as.numeric(results$log2FoldChange)
  results$p <- as.numeric(results$p)
  results$padj <- as.numeric(results$padj)
  results$se <- as.numeric(results$se)
  results$significance <- ifelse(
    !is.na(results$padj) &
      results$padj < alpha &
      (abs(results$log2FoldChange) > lfc_cutoff | !is.na(results$struc0)),
    "Sig",
    "NotSig"
  )
  results$direction <- ifelse(
    results$log2FoldChange > 0 | results$struc0 == "group1",
    "pos",
    ifelse(results$log2FoldChange < 0 | results$struc0 == "group2", "neg", "none")
  )

  results[, c("taxon", "log2FoldChange", "p", "padj", "struc0", "se", "significance", "direction")]
}

run_ancombc_comparison <- function(physeq, spec, global_config) {
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    stop("The R package 'ANCOMBC' is required for ANCOMBC differential abundance.", call. = FALSE)
  }

  tax_agg_level <- if (!is.null(spec$tax_agg_level)) spec$tax_agg_level else global_config$tax_agg_level
  alpha <- if (!is.null(spec$alpha)) spec$alpha else global_config$alpha
  lfc_cutoff <- if (!is.null(spec$lfc_cutoff)) spec$lfc_cutoff else global_config$lfc_cutoff

  ancombc_output <- ANCOMBC::ancombc(
    data = physeq,
    tax_level = tax_agg_level,
    formula = spec$formula,
    p_adj_method = "holm",
    group = spec$group,
    struc_zero = TRUE,
    alpha = alpha
  )

  standardize_ancombc_results(ancombc_output, spec, alpha, lfc_cutoff)
}

write_da_results <- function(results_df, out_dir, trial_id, method, comparison_name) {
  comparison_dir <- file.path(out_dir, method, comparison_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(
    comparison_dir,
    paste0(trial_id, "_", comparison_name, "_", method, "Results.tsv")
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

