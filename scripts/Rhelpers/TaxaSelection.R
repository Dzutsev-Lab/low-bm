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
