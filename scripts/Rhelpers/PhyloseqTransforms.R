source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "MetadataSchema.R"))

otu_samples_by_taxa <- function(physeq) {
  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
  if (phyloseq::taxa_are_rows(physeq)) {
    otu_mat <- t(otu_mat)
  }
  as.matrix(otu_mat)
}

set_otu_samples_by_taxa <- function(physeq, otu_mat) {
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(as.matrix(otu_mat), taxa_are_rows = FALSE)
  physeq
}

prune_empty_physeq <- function(physeq) {
  physeq <- phyloseq::prune_samples(phyloseq::sample_sums(physeq) > 0, physeq)
  phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
}

tax_glom_rename <- function(physeq, tax_agg_level = NULL) {
  if (is.null(tax_agg_level) ||
      is.na(tax_agg_level) ||
      tax_agg_level %in% c("", "none", "None", "ASV")) {
    return(physeq)
  }

  tax_df <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  if (!tax_agg_level %in% colnames(tax_df)) {
    stop("Taxonomic rank not found in tax_table: ", tax_agg_level, call. = FALSE)
  }

  physeq <- phyloseq::tax_glom(physeq, taxrank = tax_agg_level, NArm = FALSE)
  labels <- as.character(phyloseq::tax_table(physeq)[, tax_agg_level])
  labels[is.na(labels) | labels == ""] <- phyloseq::taxa_names(physeq)[is.na(labels) | labels == ""]
  phyloseq::taxa_names(physeq) <- make.unique(labels)
  physeq
}

average_by_techrep <- function(physeq, sample_id_col = "SampleID") {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), sample_id_col, "phyloseq sample_data")

  otu_mat <- otu_samples_by_taxa(physeq)
  group <- as.character(metadata_df[[sample_id_col]])
  group_factor <- factor(group, levels = unique(group))
  group_counts <- table(group_factor)

  avg_otu_mat <- rowsum(otu_mat, group = group_factor, reorder = FALSE)
  avg_otu_mat <- sweep(
    avg_otu_mat,
    1,
    as.numeric(group_counts[rownames(avg_otu_mat)]),
    FUN = "/"
  )

  avg_metadata_df <- metadata_df[match(rownames(avg_otu_mat), group), , drop = FALSE]
  rownames(avg_metadata_df) <- rownames(avg_otu_mat)

  phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(avg_otu_mat), taxa_are_rows = FALSE),
    phyloseq::sample_data(avg_metadata_df),
    phyloseq::tax_table(physeq)
  )
}

divide_by_sample_factor <- function(physeq,
                                    factor_column,
                                    scale = 1e6,
                                    drop_invalid = TRUE) {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), factor_column, "phyloseq sample_data")

  factors <- suppressWarnings(as.numeric(metadata_df[[factor_column]]))
  bad <- is.na(factors) | !is.finite(factors) | factors <= 0
  if (any(bad)) {
    bad_samples <- rownames(metadata_df)[bad]
    if (!drop_invalid) {
      stop(
        "Cannot normalize by ",
        factor_column,
        ": all values must be positive and finite.",
        call. = FALSE
      )
    }
    if (all(bad)) {
      stop(
        "Cannot normalize by ",
        factor_column,
        ": all samples have missing, non-finite, zero, or negative values.",
        call. = FALSE
      )
    }

    warning(
      "Dropping ",
      length(bad_samples),
      " sample(s) before ",
      factor_column,
      " normalization because their divisor is missing, non-finite, zero, or negative: ",
      paste(bad_samples, collapse = ", "),
      call. = FALSE
    )
    physeq <- phyloseq::prune_samples(!bad, physeq)
    metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
    factors <- suppressWarnings(as.numeric(metadata_df[[factor_column]]))
  }

  otu_mat <- otu_samples_by_taxa(physeq)
  otu_mat <- sweep(otu_mat, 1, factors, FUN = "/") * scale
  set_otu_samples_by_taxa(physeq, otu_mat)
}

counts_normalization <- function(physeq, norm_method = "noNorm", pseudocount = 1) {
  if (is.null(norm_method) || norm_method %in% c("", "noNorm")) {
    return(physeq)
  }

  if (norm_method == "log2") {
    return(phyloseq::transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
  }
  if (norm_method == "RelAbund") {
    return(phyloseq::transform_sample_counts(physeq, function(x) {
      denom <- sum(x)
      if (denom == 0) x else x / denom
    }))
  }
  if (norm_method == "RawTSS") {
    return(divide_by_sample_factor(physeq, "Raw_reads"))
  }
  if (norm_method == "HostMapped") {
    return(divide_by_sample_factor(physeq, "Host_mapped_reads"))
  }
  if (norm_method == "log2HostMapped") {
    physeq <- divide_by_sample_factor(physeq, "Host_mapped_reads")
    return(phyloseq::transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
  }
  if (norm_method == "log2RelAbund") {
    physeq <- phyloseq::transform_sample_counts(physeq, function(x) {
      y <- x + pseudocount
      y / sum(y)
    })
    return(phyloseq::transform_sample_counts(physeq, log2))
  }

  stop("Unknown normalization method: ", norm_method, call. = FALSE)
}

batch_adjustment <- function(physeq,
                             batch_column = NULL,
                             method = c("removeBatchEffect", "ComBat"),
                             design_formula = "~ SampleType") {
  method <- match.arg(method)
  if (is.null(batch_column) || is.na(batch_column) || batch_column == "") {
    message("Skipping batch adjustment: no batch column provided.")
    return(physeq)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  if (!batch_column %in% names(metadata_df)) {
    message("Skipping batch adjustment: metadata column not found: ", batch_column)
    return(physeq)
  }

  batch <- factor(ifelse(is.na(metadata_df[[batch_column]]), "Unknown", metadata_df[[batch_column]]))
  if (length(levels(batch)) < 2) {
    message("Skipping batch adjustment: fewer than two batch levels.")
    return(physeq)
  }

  otu_mat <- t(otu_samples_by_taxa(physeq))
  design <- stats::model.matrix(stats::as.formula(design_formula), data = metadata_df)

  if (method == "removeBatchEffect") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("The R package 'limma' is required for removeBatchEffect.", call. = FALSE)
    }
    adjusted <- limma::removeBatchEffect(otu_mat, batch = batch, design = design)
  } else {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("The R package 'sva' is required for ComBat.", call. = FALSE)
    }
    adjusted <- sva::ComBat(dat = otu_mat, batch = batch, mod = NULL, par.prior = TRUE)
  }

  set_otu_samples_by_taxa(physeq, t(adjusted))
}

limma_voom_normalization <- function(physeq, design_formula = "~ SampleType + PatientID") {
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE)) {
    stop("The R packages 'limma' and 'edgeR' are required for voom normalization.", call. = FALSE)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  counts <- t(otu_samples_by_taxa(physeq))
  design <- stats::model.matrix(stats::as.formula(design_formula), data = metadata_df)

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  voom_out <- limma::voom(dge, design, plot = FALSE)

  set_otu_samples_by_taxa(physeq, t(as.matrix(voom_out$E)))
}
