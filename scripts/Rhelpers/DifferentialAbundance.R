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
    common$sample_filter <- list(SampleType = c("Tumor", "Nontumor"))
    common$factor_levels <- list(SampleType = c("Nontumor", "Tumor"))
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

normalize_da_method <- function(method) {
  method <- toupper(as.character(method))
  method[method == "ANCOMBC"] <- "ANCOMBC2"
  method
}

normalize_ancombc2_tests <- function(tests = NULL) {
  if (is.null(tests) || length(tests) == 0) {
    return("primary")
  }
  tests <- tolower(as.character(unlist(tests, use.names = FALSE)))
  aliases <- c(
    reference = "primary",
    contrast = "primary",
    contrasts = "primary",
    global_test = "global",
    pair = "pairwise",
    pairs = "pairwise",
    dunnett = "dunnet",
    dunnett_type = "dunnet",
    pattern = "trend",
    patterns = "trend"
  )
  tests <- ifelse(tests %in% names(aliases), aliases[tests], tests)
  supported <- c("primary", "global", "pairwise", "dunnet", "trend")
  unsupported <- setdiff(tests, supported)
  if (length(unsupported) > 0) {
    stop(
      "Unsupported ANCOM-BC2 test(s): ",
      paste(unsupported, collapse = ", "),
      ". Supported tests: ",
      paste(supported, collapse = ", "),
      call. = FALSE
    )
  }
  unique(tests)
}

normalize_da_config <- function(config) {
  if (is.null(config$methods)) config$methods <- c("ANCOMBC2")
  config$methods <- normalize_da_method(config$methods)
  if (is.null(config$norm_method)) config$norm_method <- "noNorm"
  if (is.null(config$pseudocount)) config$pseudocount <- 1
  if (is.null(config$tax_agg_level)) config$tax_agg_level <- "Genus"
  if (is.null(config$tax_label_level)) config$tax_label_level <- config$tax_agg_level
  if (is.null(config$alpha)) config$alpha <- 0.05
  if (is.null(config$lfc_cutoff)) config$lfc_cutoff <- 1
  if (is.null(config$p_adj_method)) config$p_adj_method <- "holm"
  if (is.null(config$pseudo_sens)) config$pseudo_sens <- TRUE
  if (is.null(config$prv_cut)) config$prv_cut <- 0.10
  if (is.null(config$lib_cut)) config$lib_cut <- 1000
  if (is.null(config$s0_perc)) config$s0_perc <- 0.05
  if (is.null(config$struc_zero)) config$struc_zero <- TRUE
  if (is.null(config$neg_lb)) config$neg_lb <- TRUE
  if (is.null(config$n_cl)) config$n_cl <- 1
  if (is.null(config$verbose)) config$verbose <- TRUE
  if (is.null(config$mdfdr_control)) {
    config$mdfdr_control <- list(fwer_ctrl_method = "holm", B = 100)
  }
  config$patient_duplicate_policy <- normalize_patient_duplicate_policy(
    config$patient_duplicate_policy,
    default_action = "keep",
    allowed_actions = c("keep", "drop", "error"),
    context = "differential_abundance.patient_duplicate_policy"
  )

  if (is.null(config$comparisons) || length(config$comparisons) == 0) {
    stop("DA config must define at least one comparison.", call. = FALSE)
  }

  unsupported <- setdiff(config$methods, "ANCOMBC2")
  if (length(unsupported) > 0) {
    stop(
      "Only ANCOMBC2 is supported by the config-driven DA interface. Unsupported method(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  config$comparisons <- lapply(config$comparisons, function(spec) {
    spec <- spec %||% list()
    spec_name <- if (!is.null(spec$name) && length(spec$name) > 0) as.character(spec$name)[[1]] else "comparison"
    spec$method <- normalize_da_method(spec$method %||% config$methods[[1]])
    if (!identical(spec$method, "ANCOMBC2")) {
      stop(
        "Only ANCOMBC2 is supported for DA comparison '",
        spec_name,
        "'. Unsupported method: ",
        spec$method,
        call. = FALSE
      )
    }
    if (is.null(spec$fix_formula) && !is.null(spec$formula)) {
      spec$fix_formula <- spec$formula
    }
    spec$tests <- normalize_ancombc2_tests(spec$tests)
    spec$patient_duplicate_policy <- normalize_patient_duplicate_policy(
      spec$patient_duplicate_policy %||% config$patient_duplicate_policy,
      default_action = config$patient_duplicate_policy$action,
      allowed_actions = c("keep", "drop", "error"),
      context = paste0("differential_abundance.comparisons['", spec_name, "'].patient_duplicate_policy")
    )
    spec
  })

  config
}

is_missing_da_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
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

da_fix_formula <- function(spec) {
  spec$fix_formula %||% spec$formula
}

da_formula_uses_patient_id <- function(spec, patient_id_col = "PatientID") {
  formula_value <- da_fix_formula(spec)
  if (is_missing_da_value(formula_value)) {
    return(FALSE)
  }
  formula_text <- as.character(formula_value)[[1]]
  formula_terms <- tryCatch(
    all.vars(stats::as.formula(paste0("~ ", formula_text))),
    error = function(e) character(0)
  )
  patient_id_col %in% formula_terms
}

apply_da_patient_duplicate_policy <- function(physeq,
                                              spec,
                                              global_config,
                                              patient_id_col = "PatientID") {
  policy <- spec$patient_duplicate_policy %||% global_config$patient_duplicate_policy
  policy <- normalize_patient_duplicate_policy(
    policy,
    default_action = "keep",
    allowed_actions = c("keep", "drop", "error"),
    context = "differential_abundance.patient_duplicate_policy"
  )

  if (!da_formula_uses_patient_id(spec, patient_id_col = patient_id_col) ||
      is_missing_da_value(spec$group)) {
    return(list(
      physeq = physeq,
      audit = empty_patient_duplicate_policy_audit(),
      policy = policy
    ))
  }

  comparison_name <- if (!is.null(spec$name) && length(spec$name) > 0) as.character(spec$name)[[1]] else "comparison"
  group_col <- as.character(spec$group)[[1]]
  apply_patient_duplicate_policy_physeq(
    physeq,
    policy = policy,
    patient_id_col = patient_id_col,
    unit_cols = group_col,
    default_action = "keep",
    allowed_actions = c("keep", "drop", "error"),
    context = paste0("DA comparison '", comparison_name, "'")
  )
}

inferred_structural_zero_groups <- function(group, levels) {
  paste0("structural_zero (", group, " = ", levels, ")")
}

ancombc_direction_labels <- function(spec) {
  default_labels <- list(
    x = "Effect size: log2(Fold Change)",
    caption = NULL,
    legend_title = "Direction",
    legend_labels = c(neg = "Negative", pos = "Positive", none = "No direction")
  )

  if (is.null(spec$group) ||
      is.null(spec$factor_levels) ||
      is.null(names(spec$factor_levels)) ||
      !spec$group %in% names(spec$factor_levels)) {
    return(default_labels)
  }

  group <- as.character(spec$group)[[1]]
  levels <- as.character(unlist(spec$factor_levels[[group]], use.names = FALSE))
  if (length(levels) != 2) {
    return(default_labels)
  }

  list(
    x = paste0("Effect size: log2(", levels[[2]], " / ", levels[[1]], ")"),
    caption = paste0(
      "Negative/left = ",
      group,
      " ",
      levels[[1]],
      "; positive/right = ",
      group,
      " ",
      levels[[2]],
      "."
    ),
    legend_title = group,
    legend_labels = c(
      neg = paste0("Negative: ", levels[[1]]),
      pos = paste0("Positive: ", levels[[2]]),
      none = "No direction"
    )
  )
}

stop_da_inferred_field_conflict <- function(spec_name,
                                            field,
                                            inferred_label,
                                            configured_label) {
  stop(
    "DA comparison '",
    spec_name,
    "' has configured ",
    field,
    " that does not match the inferred value after sample filtering. ",
    "Inferred ",
    field,
    ": ",
    inferred_label,
    "; configured ",
    field,
    ": ",
    configured_label,
    ". Remove ",
    field,
    " from the config or correct it.",
    call. = FALSE
  )
}

resolve_configured_group_levels <- function(spec, observed_levels) {
  group <- as.character(spec$group)[[1]]
  levels <- NULL

  if (!is.null(spec$ordered_levels) && length(spec$ordered_levels) > 0) {
    levels <- as.character(unlist(spec$ordered_levels, use.names = FALSE))
  } else if (!is.null(spec$factor_levels) &&
             length(spec$factor_levels) > 0 &&
             !is.null(names(spec$factor_levels)) &&
             group %in% names(spec$factor_levels)) {
    levels <- as.character(unlist(spec$factor_levels[[group]], use.names = FALSE))
  }

  if (is.null(levels)) {
    levels <- sort(observed_levels)
  }

  levels <- levels[!is.na(levels) & nzchar(trimws(levels))]
  missing_observed <- setdiff(observed_levels, levels)
  if (length(missing_observed) > 0) {
    stop(
      "Observed value(s) in ",
      group,
      " are absent from configured ANCOM-BC2 group levels: ",
      paste(missing_observed, collapse = ", "),
      call. = FALSE
    )
  }
  missing_configured <- setdiff(levels, observed_levels)
  if (length(missing_configured) > 0) {
    stop(
      "Configured ANCOM-BC2 group level(s) in comparison '",
      spec$name,
      "' were not observed after sample filtering: ",
      paste(missing_configured, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is_missing_da_value(spec$reference_level)) {
    reference_level <- as.character(spec$reference_level)[[1]]
    if (!reference_level %in% levels) {
      stop(
        "ANCOM-BC2 reference_level for comparison '",
        spec$name,
        "' is not among the resolved group levels: ",
        reference_level,
        call. = FALSE
      )
    }
    levels <- c(reference_level, levels[levels != reference_level])
  }

  levels
}

validate_ancombc2_trend_spec <- function(spec, levels) {
  if (!"trend" %in% spec$tests) {
    return(invisible(TRUE))
  }
  if (is.null(spec$ordered_levels) || length(spec$ordered_levels) == 0) {
    stop(
      "ANCOM-BC2 trend comparison '",
      spec$name,
      "' must define ordered_levels. Alphabetical ordering is not inferred for trend tests.",
      call. = FALSE
    )
  }
  ordered_levels <- as.character(unlist(spec$ordered_levels, use.names = FALSE))
  if (!identical(ordered_levels, levels)) {
    stop(
      "ANCOM-BC2 trend comparison '",
      spec$name,
      "' must use ordered_levels as the resolved factor order. Resolved levels: ",
      paste(levels, collapse = ", "),
      "; ordered_levels: ",
      paste(ordered_levels, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_legacy_factor_levels <- function(spec, inferred_levels) {
  if (is.null(spec$factor_levels) || length(spec$factor_levels) == 0) {
    return(invisible(TRUE))
  }

  group <- as.character(spec$group)[[1]]
  factor_level_names <- names(spec$factor_levels)
  if (is.null(factor_level_names) || !group %in% factor_level_names) {
    stop_da_inferred_field_conflict(
      as.character(spec$name)[[1]],
      "factor_levels",
      paste0(group, " = [", paste(inferred_levels, collapse = ", "), "]"),
      if (is.null(factor_level_names)) "<unnamed>" else paste(factor_level_names, collapse = ", ")
    )
  }

  configured_levels <- as.character(unlist(spec$factor_levels[[group]], use.names = FALSE))
  if (!setequal(configured_levels, inferred_levels)) {
    stop_da_inferred_field_conflict(
      as.character(spec$name)[[1]],
      "factor_levels",
      paste0(group, " = [", paste(inferred_levels, collapse = ", "), "]"),
      paste0(group, " = [", paste(configured_levels, collapse = ", "), "]")
    )
  }

  invisible(TRUE)
}

validate_legacy_coefficient <- function(spec, inferred_coefficient) {
  if (is.null(spec$coefficient) || length(spec$coefficient) == 0) {
    return(invisible(TRUE))
  }

  configured_coefficient <- as.character(spec$coefficient)
  if (!identical(configured_coefficient, inferred_coefficient)) {
    stop_da_inferred_field_conflict(
      as.character(spec$name)[[1]],
      "coefficient",
      inferred_coefficient,
      paste(configured_coefficient, collapse = ", ")
    )
  }

  invisible(TRUE)
}

validate_legacy_structural_zero_groups <- function(spec, inferred_groups) {
  if (is.null(spec$structural_zero_groups) || length(spec$structural_zero_groups) == 0) {
    return(invisible(TRUE))
  }

  configured_groups <- as.character(unlist(spec$structural_zero_groups, use.names = FALSE))
  if (!identical(configured_groups, inferred_groups)) {
    stop_da_inferred_field_conflict(
      as.character(spec$name)[[1]],
      "structural_zero_groups",
      paste(inferred_groups, collapse = "; "),
      paste(configured_groups, collapse = "; ")
    )
  }

  invisible(TRUE)
}

resolve_ancombc_comparison_spec <- function(physeq, spec) {
  if (is.null(spec$fix_formula) && !is.null(spec$formula)) {
    spec$fix_formula <- spec$formula
  }
  if (is.null(spec$tests)) {
    spec$tests <- normalize_ancombc2_tests(spec$tests)
  }

  required <- c("name", "fix_formula", "group")
  missing <- required[!vapply(required, function(x) x %in% names(spec) && !is_missing_da_value(spec[[x]]), logical(1))]
  if (length(missing) > 0) {
    stop(
      "DA comparison spec is missing required field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  spec$name <- as.character(spec$name)[[1]]
  spec$fix_formula <- as.character(spec$fix_formula)[[1]]
  spec$formula <- spec$fix_formula
  spec$group <- as.character(spec$group)[[1]]
  spec$tests <- normalize_ancombc2_tests(spec$tests)

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  formula_terms <- all.vars(stats::as.formula(paste0("~ ", spec$fix_formula)))
  fail_missing_columns(names(metadata_df), unique(c(spec$group, formula_terms)), "phyloseq sample_data")

  grouping_values <- as.character(metadata_df[[spec$group]])
  grouping_values <- grouping_values[!is.na(grouping_values) & nzchar(trimws(grouping_values))]
  observed_levels <- sort(unique(grouping_values))
  if (length(observed_levels) < 2) {
    stop(
      "DA comparison '",
      spec$name,
      "' must have at least two observed non-empty group levels in ",
      spec$group,
      " after sample filtering; found ",
      length(observed_levels),
      if (length(observed_levels) > 0) paste0(": ", paste(observed_levels, collapse = ", ")) else "",
      ". Update sample_filter or factor/ordered levels.",
      call. = FALSE
    )
  }

  inferred_levels <- resolve_configured_group_levels(spec, observed_levels)
  inferred_coefficient <- paste0(spec$group, inferred_levels[[2]])
  inferred_struc0_groups <- inferred_structural_zero_groups(spec$group, inferred_levels)

  validate_legacy_factor_levels(spec, inferred_levels)
  if (length(inferred_levels) == 2) {
    validate_legacy_coefficient(spec, inferred_coefficient)
    validate_legacy_structural_zero_groups(spec, inferred_struc0_groups)
  }
  validate_ancombc2_trend_spec(spec, inferred_levels)

  spec$factor_levels <- stats::setNames(list(inferred_levels), spec$group)
  spec$coefficient <- if (length(inferred_levels) == 2) inferred_coefficient else NULL
  spec$structural_zero_groups <- inferred_struc0_groups
  spec$reference_level <- inferred_levels[[1]]

  metadata_df[[spec$group]] <- factor(as.character(metadata_df[[spec$group]]), levels = inferred_levels)
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(metadata_df)

  list(physeq = physeq, spec = spec)
}

validate_comparison_spec <- function(physeq, spec) {
  resolve_ancombc_comparison_spec(physeq, spec)
  invisible(TRUE)
}

prepare_da_physeq <- function(physeq, spec, global_config, select_taxa = NULL) {
  physeq <- apply_sample_filter(physeq, spec$sample_filter)
  duplicate_policy <- apply_da_patient_duplicate_policy(
    physeq,
    spec = spec,
    global_config = global_config,
    patient_id_col = "PatientID"
  )
  physeq <- duplicate_policy$physeq
  resolved <- resolve_ancombc_comparison_spec(physeq, spec)
  physeq <- resolved$physeq

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

  attr(physeq, "da_comparison_spec") <- resolved$spec
  attr(physeq, "patient_duplicate_policy_audit") <- duplicate_policy$audit
  attr(physeq, "patient_duplicate_policy") <- duplicate_policy$policy
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

assign_da_direction <- function(log2_fold_change, struc0) {
  direction <- rep("none", length(log2_fold_change))
  struc0 <- as.character(struc0)
  pos <- (!is.na(log2_fold_change) & log2_fold_change > 0) |
    (!is.na(struc0) & struc0 == "group1")
  neg <- (!is.na(log2_fold_change) & log2_fold_change < 0) |
    (!is.na(struc0) & struc0 == "group2")
  direction[pos] <- "pos"
  direction[neg] <- "neg"
  direction
}

logical_or_na <- function(x) {
  if (is.null(x)) {
    return(NA)
  }
  as.logical(x)
}

numeric_or_na <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}

ancombc2_group_levels <- function(spec) {
  if (!is.null(spec$factor_levels) &&
      !is.null(names(spec$factor_levels)) &&
      spec$group %in% names(spec$factor_levels)) {
    return(as.character(unlist(spec$factor_levels[[spec$group]], use.names = FALSE)))
  }
  character(0)
}

ancombc2_contrast_specs <- function(spec, test = "primary") {
  levels <- ancombc2_group_levels(spec)
  if (length(levels) < 2) {
    return(data.frame(
      coefficient = character(0),
      contrast = character(0),
      reference_level = character(0),
      target_level = character(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- list()
  group <- spec$group
  if (test %in% c("primary", "dunnet")) {
    for (i in 2:length(levels)) {
      rows[[length(rows) + 1]] <- data.frame(
        coefficient = paste0(group, levels[[i]]),
        contrast = paste0(levels[[i]], " vs ", levels[[1]]),
        reference_level = levels[[1]],
        target_level = levels[[i]],
        stringsAsFactors = FALSE
      )
    }
  } else if (identical(test, "pairwise")) {
    for (i in seq_len(length(levels) - 1)) {
      for (j in (i + 1):length(levels)) {
        rows[[length(rows) + 1]] <- data.frame(
          coefficient = if (i == 1) {
            paste0(group, levels[[j]])
          } else {
            paste0(group, levels[[j]], "_", group, levels[[i]])
          },
          contrast = paste0(levels[[j]], " vs ", levels[[i]]),
          reference_level = levels[[i]],
          target_level = levels[[j]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      coefficient = character(0),
      contrast = character(0),
      reference_level = character(0),
      target_level = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

ancombc2_col <- function(table, prefix, coefficient, fallback = NULL) {
  candidates <- unique(c(
    paste0(prefix, "_", coefficient),
    paste0(prefix, coefficient),
    fallback
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  match <- intersect(candidates, names(table))
  if (length(match) == 0) NULL else match[[1]]
}

ancombc2_vector <- function(table, column, default = NA_real_) {
  if (is.null(column) || !column %in% names(table)) {
    return(rep(default, nrow(table)))
  }
  table[[column]]
}

format_structural_zero_levels <- function(ancombc_output, spec) {
  zero_ind <- ancombc_output[["zero_ind"]]
  levels <- ancombc2_group_levels(spec)
  structural_zero_groups <- inferred_structural_zero_groups(spec$group, levels)
  if (is.null(zero_ind) || length(levels) == 0) {
    return(data.frame(taxon = character(0), struc0 = character(0), stringsAsFactors = FALSE))
  }

  required <- c("taxon", structural_zero_groups)
  available <- intersect(required, colnames(zero_ind))
  if (!"taxon" %in% available || length(available) == 1) {
    return(data.frame(taxon = character(0), struc0 = character(0), stringsAsFactors = FALSE))
  }

  zero_df <- zero_ind[, available, drop = FALSE]
  group_cols <- setdiff(names(zero_df), "taxon")
  struc0 <- vapply(seq_len(nrow(zero_df)), function(i) {
    flags <- as.logical(zero_df[i, group_cols, drop = TRUE])
    flagged <- levels[match(group_cols[flags], structural_zero_groups)]
    flagged <- flagged[!is.na(flagged)]
    if (length(flagged) == 0) {
      return(NA_character_)
    }
    if (length(levels) == 2 && length(flagged) == 1) {
      return(if (identical(flagged[[1]], levels[[1]])) "group1" else "group2")
    }
    paste(flagged, collapse = ";")
  }, character(1))

  data.frame(taxon = zero_df$taxon, struc0 = struc0, stringsAsFactors = FALSE)
}

empty_ancombc2_results <- function() {
  data.frame(
    taxon = character(0),
    test = character(0),
    contrast = character(0),
    coefficient = character(0),
    reference_level = character(0),
    target_level = character(0),
    log2FoldChange = numeric(0),
    p = numeric(0),
    padj = numeric(0),
    struc0 = character(0),
    se = numeric(0),
    W = numeric(0),
    diff_abn = logical(0),
    passed_ss = logical(0),
    diff_robust = logical(0),
    significance = character(0),
    direction = character(0),
    stringsAsFactors = FALSE
  )
}

ancombc2_significance <- function(log2_fold_change,
                                  p,
                                  padj,
                                  diff_abn,
                                  passed_ss,
                                  diff_robust,
                                  struc0,
                                  alpha,
                                  lfc_cutoff) {
  base_sig <- ifelse(
    !is.na(diff_robust),
    diff_robust,
    ifelse(
      !is.na(diff_abn),
      diff_abn & (is.na(passed_ss) | passed_ss),
      !is.na(padj) & padj < alpha
    )
  )
  has_effect <- is.na(log2_fold_change) |
    abs(log2_fold_change) > lfc_cutoff |
    (!is.na(struc0) & nzchar(struc0))
  ifelse(base_sig & has_effect, "Sig", "NotSig")
}

standardize_ancombc2_contrast_results <- function(table,
                                                  spec,
                                                  test,
                                                  alpha,
                                                  lfc_cutoff,
                                                  structural_zero_df = NULL) {
  if (is.null(table) || nrow(table) == 0) {
    return(empty_ancombc2_results())
  }
  contrasts <- ancombc2_contrast_specs(spec, test)
  if (nrow(contrasts) == 0) {
    return(empty_ancombc2_results())
  }

  rows <- lapply(seq_len(nrow(contrasts)), function(i) {
    contrast <- contrasts[i, , drop = FALSE]
    coefficient <- contrast$coefficient[[1]]
    lfc_col <- ancombc2_col(table, "lfc", coefficient)
    se_col <- ancombc2_col(table, "se", coefficient)
    w_col <- ancombc2_col(table, "W", coefficient)
    p_col <- ancombc2_col(table, "p", coefficient)
    q_col <- ancombc2_col(table, "q", coefficient)
    diff_col <- ancombc2_col(table, "diff", coefficient)
    passed_col <- ancombc2_col(table, "passed_ss", coefficient)
    robust_col <- ancombc2_col(table, "diff_robust", coefficient)

    out <- data.frame(
      taxon = table$taxon,
      test = test,
      contrast = contrast$contrast[[1]],
      coefficient = coefficient,
      reference_level = contrast$reference_level[[1]],
      target_level = contrast$target_level[[1]],
      log2FoldChange = numeric_or_na(ancombc2_vector(table, lfc_col)) / log(2),
      p = numeric_or_na(ancombc2_vector(table, p_col)),
      padj = numeric_or_na(ancombc2_vector(table, q_col)),
      struc0 = NA_character_,
      se = numeric_or_na(ancombc2_vector(table, se_col)) / log(2),
      W = numeric_or_na(ancombc2_vector(table, w_col)),
      diff_abn = logical_or_na(ancombc2_vector(table, diff_col, default = NA)),
      passed_ss = logical_or_na(ancombc2_vector(table, passed_col, default = NA)),
      diff_robust = logical_or_na(ancombc2_vector(table, robust_col, default = NA)),
      stringsAsFactors = FALSE
    )
    if (!is.null(structural_zero_df) && nrow(structural_zero_df) > 0) {
      out <- merge(out, structural_zero_df, by = "taxon", all.x = TRUE, suffixes = c("", ".zero"))
      out$struc0 <- out$struc0.zero %||% out$struc0
      out$struc0.zero <- NULL
    }
    out$significance <- ancombc2_significance(
      out$log2FoldChange,
      out$p,
      out$padj,
      out$diff_abn,
      out$passed_ss,
      out$diff_robust,
      out$struc0,
      alpha,
      lfc_cutoff
    )
    out$direction <- assign_da_direction(out$log2FoldChange, out$struc0)
    out
  })

  do.call(rbind, rows)
}

standardize_ancombc2_omnibus_results <- function(table,
                                                 spec,
                                                 test,
                                                 alpha,
                                                 lfc_cutoff,
                                                 contrast = test) {
  if (is.null(table) || nrow(table) == 0) {
    return(empty_ancombc2_results())
  }
  p_col <- if ("p_val" %in% names(table)) "p_val" else if ("p" %in% names(table)) "p" else NULL
  q_col <- if ("q_val" %in% names(table)) "q_val" else if ("q" %in% names(table)) "q" else NULL
  diff_col <- if ("diff_robust_abn" %in% names(table)) "diff_abn" else if ("diff_abn" %in% names(table)) "diff_abn" else NULL
  robust_col <- if ("diff_robust_abn" %in% names(table)) "diff_robust_abn" else NULL
  passed_col <- if ("passed_ss" %in% names(table)) "passed_ss" else NULL
  out <- data.frame(
    taxon = table$taxon,
    test = test,
    contrast = contrast,
    coefficient = NA_character_,
    reference_level = spec$reference_level %||% NA_character_,
    target_level = NA_character_,
    log2FoldChange = NA_real_,
    p = numeric_or_na(ancombc2_vector(table, p_col)),
    padj = numeric_or_na(ancombc2_vector(table, q_col)),
    struc0 = NA_character_,
    se = NA_real_,
    W = numeric_or_na(ancombc2_vector(table, if ("W" %in% names(table)) "W" else NULL)),
    diff_abn = logical_or_na(ancombc2_vector(table, diff_col, default = NA)),
    passed_ss = logical_or_na(ancombc2_vector(table, passed_col, default = NA)),
    diff_robust = logical_or_na(ancombc2_vector(table, robust_col, default = NA)),
    stringsAsFactors = FALSE
  )
  out$significance <- ancombc2_significance(
    out$log2FoldChange,
    out$p,
    out$padj,
    out$diff_abn,
    out$passed_ss,
    out$diff_robust,
    out$struc0,
    alpha,
    lfc_cutoff
  )
  out$direction <- "none"
  out
}

standardize_ancombc2_trend_results <- function(table,
                                               spec,
                                               alpha,
                                               lfc_cutoff) {
  if (is.null(table) || nrow(table) == 0) {
    return(empty_ancombc2_results())
  }
  contrasts <- ancombc2_contrast_specs(spec, "primary")
  if (nrow(contrasts) == 0) {
    return(standardize_ancombc2_omnibus_results(table, spec, "trend", alpha, lfc_cutoff, "trend"))
  }
  omnibus <- standardize_ancombc2_omnibus_results(table, spec, "trend", alpha, lfc_cutoff, "trend")
  rows <- lapply(seq_len(nrow(contrasts)), function(i) {
    contrast <- contrasts[i, , drop = FALSE]
    coefficient <- contrast$coefficient[[1]]
    lfc_col <- ancombc2_col(table, "lfc", coefficient)
    se_col <- ancombc2_col(table, "se", coefficient)
    row <- omnibus
    row$contrast <- contrast$contrast[[1]]
    row$coefficient <- coefficient
    row$reference_level <- contrast$reference_level[[1]]
    row$target_level <- contrast$target_level[[1]]
    row$log2FoldChange <- numeric_or_na(ancombc2_vector(table, lfc_col)) / log(2)
    row$se <- numeric_or_na(ancombc2_vector(table, se_col)) / log(2)
    row$direction <- assign_da_direction(row$log2FoldChange, row$struc0)
    row$significance <- ancombc2_significance(
      row$log2FoldChange,
      row$p,
      row$padj,
      row$diff_abn,
      row$passed_ss,
      row$diff_robust,
      row$struc0,
      alpha,
      lfc_cutoff
    )
    row
  })
  do.call(rbind, rows)
}

standardize_ancombc2_results <- function(ancombc_output, spec, alpha, lfc_cutoff) {
  structural_zero_df <- format_structural_zero_levels(ancombc_output, spec)
  parts <- list()
  if (!is.null(ancombc_output$res)) {
    parts[[length(parts) + 1]] <- standardize_ancombc2_contrast_results(
      ancombc_output$res,
      spec,
      "primary",
      alpha,
      lfc_cutoff,
      structural_zero_df
    )
  }
  if ("global" %in% spec$tests && !is.null(ancombc_output$res_global)) {
    parts[[length(parts) + 1]] <- standardize_ancombc2_omnibus_results(
      ancombc_output$res_global,
      spec,
      "global",
      alpha,
      lfc_cutoff,
      "global"
    )
  }
  if ("pairwise" %in% spec$tests && !is.null(ancombc_output$res_pair)) {
    parts[[length(parts) + 1]] <- standardize_ancombc2_contrast_results(
      ancombc_output$res_pair,
      spec,
      "pairwise",
      alpha,
      lfc_cutoff,
      structural_zero_df
    )
  }
  if ("dunnet" %in% spec$tests && !is.null(ancombc_output$res_dunn)) {
    parts[[length(parts) + 1]] <- standardize_ancombc2_contrast_results(
      ancombc_output$res_dunn,
      spec,
      "dunnet",
      alpha,
      lfc_cutoff,
      structural_zero_df
    )
  }
  if ("trend" %in% spec$tests && !is.null(ancombc_output$res_trend)) {
    parts[[length(parts) + 1]] <- standardize_ancombc2_trend_results(
      ancombc_output$res_trend,
      spec,
      alpha,
      lfc_cutoff
    )
  }

  parts <- parts[vapply(parts, function(x) !is.null(x) && nrow(x) > 0, logical(1))]
  if (length(parts) == 0) {
    return(empty_ancombc2_results())
  }
  results <- do.call(rbind, parts)
  rownames(results) <- NULL
  results
}

standardize_ancombc_results <- function(ancombc_output, spec, alpha, lfc_cutoff) {
  if (is.data.frame(ancombc_output[["res"]])) {
    return(standardize_ancombc2_results(ancombc_output, spec, alpha, lfc_cutoff))
  }

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
  results$direction <- assign_da_direction(results$log2FoldChange, results$struc0)

  results[, c("taxon", "log2FoldChange", "p", "padj", "struc0", "se", "significance", "direction")]
}

default_ancombc2_trend_matrix <- function(pattern, n_coef) {
  mat <- matrix(0, nrow = n_coef, ncol = n_coef)
  if (identical(pattern, "increasing")) {
    mat[1, 1] <- 1
    if (n_coef > 1) {
      for (i in 2:n_coef) {
        mat[i, i - 1] <- -1
        mat[i, i] <- 1
      }
    }
  } else if (identical(pattern, "decreasing")) {
    mat[1, 1] <- -1
    if (n_coef > 1) {
      for (i in 2:n_coef) {
        mat[i, i - 1] <- 1
        mat[i, i] <- -1
      }
    }
  } else {
    stop("Unsupported default ANCOM-BC2 trend pattern: ", pattern, call. = FALSE)
  }
  mat
}

as_trend_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }
  values <- as.numeric(unlist(x, use.names = FALSE))
  n <- sqrt(length(values))
  if (n != floor(n)) {
    stop("Explicit ANCOM-BC2 trend contrast must be a square matrix or square-length numeric vector.", call. = FALSE)
  }
  matrix(values, nrow = n, byrow = TRUE)
}

normalize_ancombc2_trend_patterns <- function(spec) {
  levels <- ancombc2_group_levels(spec)
  n_coef <- length(levels) - 1
  if (n_coef < 1) {
    return(list())
  }
  patterns <- spec$trend_patterns %||% spec$trend_pattern
  if (is.null(patterns) || length(patterns) == 0) {
    patterns <- c("increasing", "decreasing")
  }

  if (is.character(patterns)) {
    return(lapply(patterns, function(pattern) {
      pattern <- tolower(pattern)
      list(
        name = pattern,
        contrast = default_ancombc2_trend_matrix(pattern, n_coef),
        node = n_coef
      )
    }))
  }

  lapply(patterns, function(pattern) {
    if (is.character(pattern)) {
      pattern <- tolower(pattern[[1]])
      return(list(
        name = pattern,
        contrast = default_ancombc2_trend_matrix(pattern, n_coef),
        node = n_coef
      ))
    }
    pattern_name <- as.character(pattern$name %||% pattern$type %||% "custom")
    pattern_type <- tolower(as.character(pattern$type %||% pattern_name)[[1]])
    contrast <- if (!is.null(pattern$contrast)) {
      as_trend_matrix(pattern$contrast)
    } else {
      default_ancombc2_trend_matrix(pattern_type, n_coef)
    }
    node <- as.integer(pattern$node %||% if (pattern_type %in% c("increasing", "decreasing")) n_coef else 1L)
    list(name = pattern_name, contrast = contrast, node = node)
  })
}

build_ancombc2_trend_control <- function(spec, global_config) {
  trend_control <- spec$trend_control %||% global_config$trend_control %||% list()
  if (!"trend" %in% spec$tests) {
    return(trend_control)
  }
  if (!is.null(trend_control$contrast) && !is.null(trend_control$node)) {
    return(trend_control)
  }
  patterns <- normalize_ancombc2_trend_patterns(spec)
  trend_control$contrast <- lapply(patterns, `[[`, "contrast")
  trend_control$node <- lapply(patterns, `[[`, "node")
  if (is.null(trend_control$B)) {
    trend_control$B <- 100
  }
  trend_control
}

ancombc2_arg <- function(spec, global_config, name, default = NULL) {
  spec[[name, exact = TRUE]] %||% global_config[[name, exact = TRUE]] %||% default
}

run_ancombc_comparison <- function(physeq, spec, global_config) {
  ancombc_load_error <- tryCatch({
    loadNamespace("ANCOMBC")
    NULL
  }, error = function(e) e)
  if (!is.null(ancombc_load_error)) {
    stop(
      "The R package 'ANCOMBC' is required for ANCOM-BC2 differential abundance ",
      "but could not be loaded: ",
      conditionMessage(ancombc_load_error),
      call. = FALSE
    )
  }
  if (!exists("ancombc2", envir = asNamespace("ANCOMBC"), inherits = FALSE)) {
    stop(
      "The installed R package 'ANCOMBC' does not export ancombc2(). ",
      "Install a current ANCOMBC release before running differential abundance.",
      call. = FALSE
    )
  }

  tax_agg_level <- if (!is.null(spec$tax_agg_level)) spec$tax_agg_level else global_config$tax_agg_level
  alpha <- if (!is.null(spec$alpha)) spec$alpha else global_config$alpha
  lfc_cutoff <- if (!is.null(spec$lfc_cutoff)) spec$lfc_cutoff else global_config$lfc_cutoff
  resolved <- resolve_ancombc_comparison_spec(physeq, spec)
  physeq <- resolved$physeq
  spec <- resolved$spec

  ancombc_output <- ANCOMBC::ancombc2(
    data = physeq,
    tax_level = tax_agg_level,
    fix_formula = spec$fix_formula,
    rand_formula = spec$rand_formula %||% NULL,
    p_adj_method = ancombc2_arg(spec, global_config, "p_adj_method", "holm"),
    pseudo_sens = as.logical(ancombc2_arg(spec, global_config, "pseudo_sens", TRUE)),
    prv_cut = as.numeric(ancombc2_arg(spec, global_config, "prv_cut", 0.10)),
    lib_cut = as.numeric(ancombc2_arg(spec, global_config, "lib_cut", 1000)),
    s0_perc = as.numeric(ancombc2_arg(spec, global_config, "s0_perc", 0.05)),
    group = spec$group,
    struc_zero = as.logical(ancombc2_arg(spec, global_config, "struc_zero", TRUE)),
    neg_lb = as.logical(ancombc2_arg(spec, global_config, "neg_lb", TRUE)),
    alpha = alpha,
    n_cl = as.integer(ancombc2_arg(spec, global_config, "n_cl", 1)),
    verbose = as.logical(ancombc2_arg(spec, global_config, "verbose", TRUE)),
    global = "global" %in% spec$tests,
    pairwise = "pairwise" %in% spec$tests,
    dunnet = "dunnet" %in% spec$tests,
    trend = "trend" %in% spec$tests,
    iter_control = ancombc2_arg(spec, global_config, "iter_control", list(tol = 1e-2, max_iter = 20, verbose = FALSE)),
    em_control = ancombc2_arg(spec, global_config, "em_control", list(tol = 1e-5, max_iter = 100)),
    lme_control = ancombc2_arg(spec, global_config, "lme_control", NULL),
    mdfdr_control = ancombc2_arg(spec, global_config, "mdfdr_control", list(fwer_ctrl_method = "holm", B = 100)),
    trend_control = build_ancombc2_trend_control(spec, global_config)
  )

  standardize_ancombc2_results(ancombc_output, spec, alpha, lfc_cutoff)
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

write_ancombc2_results <- function(results_df, out_dir, trial_id, comparison_name) {
  method <- "ANCOMBC2"
  main_file <- write_da_results(
    results_df,
    out_dir = out_dir,
    trial_id = trial_id,
    method = method,
    comparison_name = comparison_name
  )

  comparison_dir <- file.path(out_dir, method, comparison_name)
  test_suffix <- c(
    primary = "Primary",
    global = "Global",
    pairwise = "Pairwise",
    dunnet = "Dunnett",
    trend = "Trend"
  )
  family_files <- list()
  if ("test" %in% names(results_df)) {
    for (test in intersect(names(test_suffix), unique(as.character(results_df$test)))) {
      family_df <- results_df[as.character(results_df$test) == test, , drop = FALSE]
      out_file <- file.path(
        comparison_dir,
        paste0(trial_id, "_", comparison_name, "_ANCOMBC2", test_suffix[[test]], ".tsv")
      )
      write.table(
        family_df,
        file = out_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      family_files[[test]] <- out_file
    }
  }

  c(main = main_file, unlist(family_files, use.names = TRUE))
}
