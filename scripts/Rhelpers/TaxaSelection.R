source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

rank_taxa_by_abundance <- function(physeq) {
  otu_mat <- otu_samples_by_taxa(physeq)
  taxa_totals <- colSums(otu_mat, na.rm = TRUE)
  taxa_order <- names(sort(taxa_totals, decreasing = TRUE))

  out <- data.frame(
    taxon = taxa_order,
    total_abundance = as.numeric(taxa_totals[taxa_order]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$taxon) & nzchar(trimws(out$taxon)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

select_top_taxa_by_abundance <- function(physeq,
                                         top_n = 20,
                                         positive_only = FALSE) {
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }

  ranked_taxa <- rank_taxa_by_abundance(physeq)
  if (isTRUE(positive_only)) {
    ranked_taxa <- ranked_taxa[ranked_taxa$total_abundance > 0, , drop = FALSE]
  }

  head(ranked_taxa$taxon, top_n)
}

taxa_selection_columns <- function() {
  c(
    "source",
    "comparison",
    "taxa_level",
    "taxon",
    "selection_metric",
    "selection_reason"
  )
}

is_missing_taxa_selection_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

taxa_selection_from_list <- function(taxa_by_comparison,
                                     source,
                                     taxa_level = "Genus",
                                     selection_metric = NA_character_,
                                     selection_reason = NA_character_) {
  if (is.null(taxa_by_comparison) || length(taxa_by_comparison) == 0) {
    stop("No taxa were supplied for taxa-selection output.", call. = FALSE)
  }

  rows <- lapply(names(taxa_by_comparison), function(comparison) {
    taxa <- unique(as.character(unlist(taxa_by_comparison[[comparison]], use.names = FALSE)))
    taxa <- taxa[!is.na(taxa) & nzchar(trimws(taxa))]
    if (length(taxa) == 0) {
      return(NULL)
    }
    data.frame(
      source = source,
      comparison = comparison,
      taxa_level = taxa_level,
      taxon = taxa,
      selection_metric = as.character(selection_metric),
      selection_reason = as.character(selection_reason),
      stringsAsFactors = FALSE
    )
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    stop("No non-empty taxa were supplied for taxa-selection output.", call. = FALSE)
  }
  do.call(rbind, rows)
}

normalize_taxa_selection_table <- function(selection,
                                           source = "manual",
                                           comparison = "SelectedTaxa",
                                           taxa_level = "Genus",
                                           selection_metric = NA_character_,
                                           selection_reason = NA_character_) {
  if (is.null(selection)) {
    stop("No taxa-selection data supplied.", call. = FALSE)
  }

  if (is.vector(selection) && !is.list(selection)) {
    selection <- data.frame(taxon = as.character(selection), stringsAsFactors = FALSE)
  } else {
    selection <- as.data.frame(selection, stringsAsFactors = FALSE)
  }

  if (!"taxon" %in% names(selection)) {
    stop("Taxa-selection data must include a taxon column.", call. = FALSE)
  }

  defaults <- list(
    source = source,
    comparison = comparison,
    taxa_level = taxa_level,
    selection_metric = selection_metric,
    selection_reason = selection_reason
  )
  for (column in names(defaults)) {
    if (!column %in% names(selection)) {
      selection[[column]] <- defaults[[column]]
    }
  }

  selection <- selection[, taxa_selection_columns(), drop = FALSE]
  for (column in taxa_selection_columns()) {
    selection[[column]] <- as.character(selection[[column]])
  }
  selection$taxon <- trimws(selection$taxon)
  selection <- selection[!is.na(selection$taxon) & nzchar(selection$taxon), , drop = FALSE]
  selection <- unique(selection)
  rownames(selection) <- NULL

  if (nrow(selection) == 0) {
    stop("Taxa-selection data did not contain any non-empty taxa.", call. = FALSE)
  }

  selection
}

write_taxa_selection_table <- function(selection,
                                       out_file,
                                       source = "manual",
                                       comparison = "SelectedTaxa",
                                       taxa_level = "Genus",
                                       selection_metric = NA_character_,
                                       selection_reason = NA_character_) {
  if (is.list(selection) && !is.data.frame(selection) && !is.null(names(selection))) {
    selection <- taxa_selection_from_list(
      selection,
      source = source,
      taxa_level = taxa_level,
      selection_metric = selection_metric,
      selection_reason = selection_reason
    )
  }

  selection <- normalize_taxa_selection_table(
    selection,
    source = source,
    comparison = comparison,
    taxa_level = taxa_level,
    selection_metric = selection_metric,
    selection_reason = selection_reason
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  write.table(selection, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  out_file
}

read_taxa_selection_table <- function(path,
                                      comparisons = NULL,
                                      taxa_level = NULL) {
  if (!file.exists(path)) {
    stop("Missing taxa-selection file: ", path, call. = FALSE)
  }

  selection <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  selection <- normalize_taxa_selection_table(selection)

  if (!is_missing_taxa_selection_value(taxa_level)) {
    taxa_level <- as.character(taxa_level[[1]])
    selection <- selection[selection$taxa_level == taxa_level, , drop = FALSE]
  }
  if (!is_missing_taxa_selection_value(comparisons)) {
    comparisons <- as.character(unlist(comparisons, use.names = FALSE))
    comparisons <- comparisons[!is.na(comparisons) & nzchar(trimws(comparisons))]
    selection <- selection[selection$comparison %in% comparisons, , drop = FALSE]
  }

  if (nrow(selection) == 0) {
    stop("No taxa remained after filtering taxa-selection file: ", path, call. = FALSE)
  }

  split(selection$taxon, selection$comparison)
}
