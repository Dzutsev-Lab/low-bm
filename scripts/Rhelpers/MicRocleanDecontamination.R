source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

build_microclean_metadata <- function(physeq,
                                      sample_type_column = "SampleType",
                                      control_sample_types = "NegativeControl",
                                      batch_column = "ProcessingBatch") {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)

  if (!sample_type_column %in% names(metadata_df)) {
    stop("micRoclean sample type column not found in phyloseq sample_data: ", sample_type_column, call. = FALSE)
  }

  if (batch_column %in% names(metadata_df)) {
    batch <- metadata_df[[batch_column]]
  } else {
    batch <- rep("NoProcessingBatches", nrow(metadata_df))
  }

  control_sample_types <- as.character(unlist(control_sample_types, use.names = FALSE))
  control_sample_types <- control_sample_types[!is.na(control_sample_types) & nzchar(trimws(control_sample_types))]
  if (length(control_sample_types) == 0) {
    stop("At least one micRoclean control sample type must be configured.", call. = FALSE)
  }

  sample_type <- as.character(metadata_df[[sample_type_column]])
  microclean_meta <- data.frame(
    is_control = sample_type %in% control_sample_types,
    sample_type = factor(sample_type),
    batch = factor(batch),
    row.names = rownames(metadata_df),
    stringsAsFactors = FALSE
  )

  if (!any(microclean_meta$is_control)) {
    stop(
      "No micRoclean control samples were found using ",
      sample_type_column,
      " value(s): ",
      paste(control_sample_types, collapse = ", "),
      call. = FALSE
    )
  }

  microclean_meta
}

apply_decontaminated_taxa <- function(physeq, decontaminated_count) {
  decontaminated_count <- as.matrix(decontaminated_count)
  retained_taxa <- colnames(decontaminated_count)
  retained_taxa <- retained_taxa[retained_taxa %in% phyloseq::taxa_names(physeq)]

  if (length(retained_taxa) == 0) {
    stop("micRoclean retained zero ASVs present in the input phyloseq object.", call. = FALSE)
  }

  sample_ids <- phyloseq::sample_names(physeq)
  missing_samples <- setdiff(sample_ids, rownames(decontaminated_count))
  if (length(missing_samples) > 0) {
    stop(
      "micRoclean output is missing sample(s) from the input phyloseq object: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  cleaned_physeq <- phyloseq::prune_taxa(retained_taxa, physeq)
  cleaned_counts <- decontaminated_count[sample_ids, retained_taxa, drop = FALSE]
  set_otu_samples_by_taxa(cleaned_physeq, cleaned_counts)
}

annotate_microclean_filter_report <- function(contaminant_id, physeq) {
  report_df <- as.data.frame(contaminant_id, stringsAsFactors = FALSE)
  asv_ids <- rownames(report_df)
  if (is.null(asv_ids)) {
    asv_ids <- rep(NA_character_, nrow(report_df))
  }

  report_df <- data.frame(
    ASV = asv_ids,
    report_df,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tax_table <- phyloseq::tax_table(physeq, errorIfNULL = FALSE)
  if (is.null(tax_table)) {
    return(report_df)
  }

  tax_df <- as.data.frame(as(tax_table, "matrix"), stringsAsFactors = FALSE)
  tax_df$ASV <- rownames(tax_df)
  tax_df <- tax_df[match(report_df$ASV, tax_df$ASV), , drop = FALSE]
  rownames(tax_df) <- NULL

  cbind(report_df, tax_df[, setdiff(names(tax_df), "ASV"), drop = FALSE])
}
