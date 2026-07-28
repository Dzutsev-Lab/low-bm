desired_kraken_taxa <- function() {
  c(
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
}

kraken_tax_rank_names <- function() {
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

taxa_level_prefix_addition <- function(tax_matrix) {
  tax_prefixes <- c(
    "Domain" = "d",
    "Phylum" = "p",
    "Class" = "c",
    "Order" = "o",
    "Family" = "f",
    "Genus" = "g",
    "Species" = "s",
    "Strain" = "st",
    "Substrain" = "sst"
  )

  for (i in seq_len(ncol(tax_matrix))) {
    rank <- colnames(tax_matrix)[i]
    for (j in seq_len(nrow(tax_matrix))) {
      if (!is.na(tax_matrix[j, i]) && tax_matrix[j, i] != "") {
        tax_matrix[j, i] <- paste0(tax_prefixes[[rank]], "__", tax_matrix[j, i])
      } else {
        tax_matrix[j, i] <- NA_character_
      }
    }
  }

  tax_matrix
}

unclassified_label_propagation <- function(tax_matrix) {
  is_unassigned <- function(x) is.na(x) || x == "" || x == "Unclassified"

  for (i in seq_len(nrow(tax_matrix))) {
    row <- tax_matrix[i, , drop = TRUE]
    assigned_index <- which(!vapply(row, is_unassigned, logical(1)))
    if (length(assigned_index) == 0) next

    lowest_assigned_index <- max(assigned_index)
    fill_value <- paste0("UC_", row[lowest_assigned_index])

    if (lowest_assigned_index < ncol(tax_matrix)) {
      for (j in (lowest_assigned_index + 1):ncol(tax_matrix)) {
        if (is_unassigned(row[j])) {
          tax_matrix[i, j] <- fill_value
        }
      }
    }
  }

  tax_matrix
}

default_kraken_taxonomy_fetcher <- function(ids, sqlFile, desiredTaxa) {
  taxonomizr::getTaxonomy(ids = ids, sqlFile = sqlFile, desiredTaxa = desiredTaxa)
}

build_kraken_tax_matrix <- function(kraken_info,
                                    asv_ids,
                                    sql_db,
                                    add_unclassified_prefix = FALSE,
                                    taxonomy_fetcher = default_kraken_taxonomy_fetcher) {
  asv_ids <- as.character(asv_ids)
  ranks <- kraken_tax_rank_names()
  tax_matrix <- matrix(
    NA_character_,
    nrow = length(asv_ids),
    ncol = length(ranks),
    dimnames = list(asv_ids, ranks)
  )

  if (length(asv_ids) == 0) {
    return(tax_matrix)
  }

  required_cols <- c("status", "ASVid", "taxid")
  missing_cols <- setdiff(required_cols, names(kraken_info))
  if (length(missing_cols) > 0) {
    stop(
      "Kraken classification table is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  kraken_info$ASVid <- as.character(kraken_info$ASVid)
  kraken_info$taxid <- suppressWarnings(as.integer(kraken_info$taxid))

  classified <- kraken_info[
    kraken_info$status == "C" &
      !is.na(kraken_info$taxid) &
      kraken_info$taxid != 0 &
      kraken_info$ASVid %in% asv_ids,
    c("ASVid", "taxid"),
    drop = FALSE
  ]

  duplicated_asvs <- classified$ASVid[duplicated(classified$ASVid)]
  if (length(duplicated_asvs) > 0) {
    stop(
      "Kraken classification table has duplicated classified ASV IDs: ",
      paste(unique(duplicated_asvs), collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(classified) > 0) {
    desired_taxa <- desired_kraken_taxa()
    tax_df <- as.data.frame(
      taxonomy_fetcher(
        ids = classified$taxid,
        sqlFile = sql_db,
        desiredTaxa = desired_taxa
      ),
      stringsAsFactors = FALSE
    )

    if (nrow(tax_df) != nrow(classified)) {
      stop("Taxonomy lookup returned a different row count than Kraken classifications.", call. = FALSE)
    }

    for (rank in desired_taxa) {
      if (!rank %in% names(tax_df)) {
        tax_df[[rank]] <- NA_character_
      }
    }

    tax_df <- tax_df[, desired_taxa, drop = FALSE]
    colnames(tax_df) <- ranks
    rownames(tax_df) <- classified$ASVid
    tax_matrix[rownames(tax_df), ] <- as.matrix(tax_df)
  }

  tax_matrix <- taxa_level_prefix_addition(tax_matrix)
  if (isTRUE(add_unclassified_prefix)) {
    tax_matrix <- unclassified_label_propagation(tax_matrix)
  }

  tax_matrix
}
