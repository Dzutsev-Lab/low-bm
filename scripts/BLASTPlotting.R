library(argparse)
library(data.table)
library(ggplot2)
library(ggalluvial)
library(taxonomizr)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with blast_confirmation settings")
parser$add_argument("--trial-dir",
                    type = "character",
                    help = "BlastAnalysis directory containing comparison directories")
parser$add_argument("--comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Comparison directories to process")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "Legacy argument retained for compatibility; plots are written beside BLAST hits")
parser$add_argument("--tax-db-sql",
                    type = "character",
                    help = "taxonomizr SQLite database for resolving BLAST taxids to names")
parser$add_argument("--taxa-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level used for candidate taxa")
parser$add_argument("--plot-ranks",
                    type = "character",
                    nargs = "+",
                    default = c("genus", "species"),
                    help = "Two BLAST taxonomy ranks to show after ASV in alluvial plots")
parser$add_argument("--min-pident",
                    type = "double",
                    default = 90,
                    help = "Minimum percent identity for BLAST hit inclusion")
parser$add_argument("--min-qcovs",
                    type = "double",
                    default = 90,
                    help = "Minimum query coverage for BLAST hit inclusion")
parser$add_argument("--max-evalue",
                    type = "double",
                    default = 1e-20,
                    help = "Maximum e-value for BLAST hit inclusion")
parser$add_argument("--top-n-asvs",
                    type = "integer",
                    default = 20,
                    help = "Maximum number of ASVs to show in each taxon plot")
parser$add_argument("--min-rel-abund",
                    type = "double",
                    default = 0.05,
                    help = "Minimum relative abundance within taxon for ASV plotting")

args <- parser$parse_args()

is_missing_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

as_config_number <- function(value, default) {
  if (is_missing_value(value)) {
    return(default)
  }
  as.numeric(value)
}

as_config_integer <- function(value, default) {
  if (is_missing_value(value)) {
    return(default)
  }
  as.integer(value)
}

fail_missing_input <- function(required_input, context) {
  missing_input <- names(required_input)[vapply(required_input, is_missing_value, logical(1))]
  if (length(missing_input) > 0) {
    stop(
      "Missing required ",
      context,
      ": ",
      paste(missing_input, collapse = ", "),
      call. = FALSE
    )
  }
}

canonical_rank <- function(rank) {
  rank <- tolower(trimws(as.character(rank)))
  aliases <- c(
    domain = "superkingdom",
    kingdom = "superkingdom",
    superkingdom = "superkingdom",
    phylum = "phylum",
    class = "class",
    order = "order",
    family = "family",
    genus = "genus",
    species = "species"
  )
  out <- aliases[rank]
  out[is.na(out)] <- rank[is.na(out)]
  unname(out)
}

validate_plot_ranks <- function(plot_ranks) {
  plot_ranks <- canonical_rank(unlist(plot_ranks, use.names = FALSE))
  plot_ranks <- plot_ranks[!is.na(plot_ranks) & nzchar(plot_ranks)]
  allowed <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  invalid <- setdiff(plot_ranks, allowed)

  if (length(invalid) > 0) {
    stop("Unsupported BLAST plot rank(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (length(plot_ranks) != 2) {
    stop("Provide exactly two BLAST plot ranks, such as genus and species.", call. = FALSE)
  }

  plot_ranks
}

rank_prefix <- function(rank) {
  prefixes <- c(
    superkingdom = "d",
    phylum = "p",
    class = "c",
    order = "o",
    family = "f",
    genus = "g",
    species = "s"
  )
  prefixes[[rank]]
}

rank_col <- function(rank) {
  paste0("Blast", paste0(toupper(substr(rank, 1, 1)), substr(rank, 2, nchar(rank))))
}

blast_rank_cols <- function(plot_ranks) {
  unname(vapply(plot_ranks, rank_col, character(1)))
}

rank_label <- function(rank) {
  paste0(toupper(substr(rank, 1, 1)), substr(rank, 2, nchar(rank)))
}

read_taxon_manifest <- function(manifest_file, taxa_level, taxon_path_name) {
  m <- fread(manifest_file)
  if (!"ASVid" %in% names(m)) {
    stop("Manifest is missing ASVid column: ", manifest_file, call. = FALSE)
  }

  m[, ASVid := as.character(ASVid)]
  if ("TotalCount" %in% names(m)) {
    m[, TotalCount := as.numeric(TotalCount)]
  }

  if ("RelativeAbundanceWithinTaxon" %in% names(m)) {
    m[, RelativeAbundanceWithinTaxon := as.numeric(RelativeAbundanceWithinTaxon)]
  } else if ("RelativeAbundanceWithinGenus" %in% names(m)) {
    m[, RelativeAbundanceWithinTaxon := as.numeric(RelativeAbundanceWithinGenus)]
  } else if ("TotalCount" %in% names(m)) {
    total_count <- sum(m$TotalCount, na.rm = TRUE)
    if (total_count <= 0) {
      stop("Cannot derive relative abundance from non-positive TotalCount values: ", manifest_file, call. = FALSE)
    }
    m[, RelativeAbundanceWithinTaxon := TotalCount / total_count]
  } else {
    stop(
      "Manifest is missing RelativeAbundanceWithinTaxon, RelativeAbundanceWithinGenus, or TotalCount: ",
      manifest_file,
      call. = FALSE
    )
  }

  if (!"TaxaLevel" %in% names(m)) {
    m[, TaxaLevel := taxa_level]
  }
  if (!"Taxon" %in% names(m)) {
    if (taxa_level %in% names(m)) {
      m[, Taxon := as.character(.SD[[1]]), .SDcols = taxa_level]
    } else {
      m[, Taxon := taxon_path_name]
    }
  }
  if (!"TaxonPathName" %in% names(m)) {
    m[, TaxonPathName := taxon_path_name]
  }

  m
}

empty_blast_hits <- function() {
  data.table(
    qseqid = character(0),
    sseqid = character(0),
    pident = numeric(0),
    align_length = numeric(0),
    qcovs = numeric(0),
    evalue = numeric(0),
    bitscore = numeric(0),
    staxids = integer(0)
  )
}

parse_first_taxid <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(trimws(x)) || x == "0") {
    return(NA_integer_)
  }

  vals <- unlist(strsplit(x, ";", fixed = TRUE))
  vals <- suppressWarnings(as.integer(trimws(vals)))
  vals <- vals[!is.na(vals) & vals > 0L]

  if (length(vals) == 0) NA_integer_ else vals[1]
}

extract_taxid_from_sseqid <- function(x) {
  x <- as.character(x)
  m <- regexec("kraken:taxid\\|([0-9]+)$", x)
  parts <- regmatches(x, m)[[1]]
  if (length(parts) == 2) as.integer(parts[2]) else NA_integer_
}

read_blast_hits <- function(blast_file) {
  if (!file.exists(blast_file)) {
    stop("Missing BLAST hits file: ", blast_file, call. = FALSE)
  }
  if (file.info(blast_file)$size == 0) {
    return(empty_blast_hits())
  }

  hits <- fread(
    blast_file,
    header = FALSE,
    sep = "\t",
    col.names = c(
      "qseqid",
      "sseqid",
      "pident",
      "align_length",
      "qcovs",
      "evalue",
      "bitscore",
      "staxids"
    )
  )

  if (nrow(hits) == 0) {
    return(empty_blast_hits())
  }

  hits[, `:=`(
    qseqid = as.character(qseqid),
    sseqid = as.character(sseqid),
    pident = as.numeric(pident),
    align_length = as.numeric(align_length),
    qcovs = as.numeric(qcovs),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    staxids = vapply(staxids, parse_first_taxid, integer(1))
  )]

  missing_taxids <- is.na(hits$staxids)
  if (any(missing_taxids)) {
    hits[missing_taxids, staxids := vapply(sseqid, extract_taxid_from_sseqid, integer(1))]
  }

  hits
}

filter_blast_hits <- function(hits,
                              min_pident = 90,
                              min_qcovs = 90,
                              max_evalue = 1e-20) {
  hits[pident >= min_pident &
         qcovs >= min_qcovs &
         evalue <= max_evalue]
}

resolve_taxid_labels <- function(taxids, sql_db, plot_ranks) {
  taxids <- unique(as.integer(taxids))
  taxids <- taxids[!is.na(taxids) & taxids > 0L]
  rank_cols <- blast_rank_cols(plot_ranks)

  if (length(taxids) == 0) {
    out <- data.table(taxid = integer(0))
    for (col in rank_cols) {
      out[, (col) := character(0)]
    }
    return(out)
  }

  desired_taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_df <- getTaxonomy(ids = taxids, sqlFile = sql_db, desiredTaxa = desired_taxa)
  tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
  tax_df$taxid <- taxids

  for (rank in plot_ranks) {
    col <- rank_col(rank)
    prefix <- rank_prefix(rank)
    labels <- as.character(tax_df[[rank]])
    labels[is.na(labels) | !nzchar(labels)] <- "Unresolved"
    tax_df[[col]] <- paste0(prefix, "__", labels)
  }

  as.data.table(tax_df[, c("taxid", rank_cols), drop = FALSE])
}

collapse_hits_to_weights <- function(hits, plot_ranks) {
  rank_cols <- blast_rank_cols(plot_ranks)

  if (nrow(hits) == 0) {
    out <- data.table(qseqid = character(0))
    for (col in rank_cols) {
      out[, (col) := character(0)]
    }
    out[, hit_share := numeric(0)]
    return(out)
  }

  hits <- hits[, .(
    bitscore = sum(bitscore, na.rm = TRUE),
    pident = max(pident, na.rm = TRUE),
    qcovs = max(qcovs, na.rm = TRUE),
    evalue = min(evalue, na.rm = TRUE)
  ), by = c("qseqid", rank_cols)]

  hits[, hit_share := bitscore / sum(bitscore, na.rm = TRUE), by = qseqid]
  hits[, c("qseqid", rank_cols, "hit_share"), with = FALSE]
}

no_acceptable_hit_input <- function(manifest, plot_ranks) {
  manifest_cols <- c(
    "ASVid",
    "TaxaLevel",
    "Taxon",
    "TaxonPathName",
    "RelativeAbundanceWithinTaxon"
  )
  manifest_cols <- manifest_cols[manifest_cols %in% names(manifest)]
  out <- manifest[, manifest_cols, with = FALSE]

  for (rank in plot_ranks) {
    col <- rank_col(rank)
    out[, (col) := paste0(rank_prefix(rank), "__No acceptable hit")]
  }

  out[, hit_share := 1]
  out[, flow_value := RelativeAbundanceWithinTaxon]
  out
}

build_alluvial_input <- function(manifest, hits_weighted, plot_ranks) {
  rank_cols <- blast_rank_cols(plot_ranks)
  manifest_cols <- c(
    "ASVid",
    "TaxaLevel",
    "Taxon",
    "TaxonPathName",
    "RelativeAbundanceWithinTaxon"
  )
  manifest_cols <- manifest_cols[manifest_cols %in% names(manifest)]
  hits_cols <- c("qseqid", rank_cols, "hit_share")

  out <- merge(
    manifest[, manifest_cols, with = FALSE],
    hits_weighted[, hits_cols, with = FALSE],
    by.x = "ASVid",
    by.y = "qseqid",
    all.x = TRUE
  )

  for (rank in plot_ranks) {
    col <- rank_col(rank)
    out[is.na(get(col)) | !nzchar(get(col)), (col) := paste0(rank_prefix(rank), "__No acceptable hit")]
  }
  out[is.na(hit_share), hit_share := 1]
  out[, flow_value := RelativeAbundanceWithinTaxon * hit_share]
  out
}

plot_taxon_alluvial <- function(alluvial_df,
                                taxon_label,
                                taxa_level,
                                plot_ranks,
                                out_file = NULL,
                                top_n_asvs = 20,
                                min_rel_abund = 0.05) {
  if (nrow(alluvial_df) == 0) {
    warning("No rows available for alluvial plot: ", taxon_label, call. = FALSE)
    return(invisible(NULL))
  }

  rank_cols <- blast_rank_cols(plot_ranks)
  plot_df <- alluvial_df[RelativeAbundanceWithinTaxon >= min_rel_abund]
  if (nrow(plot_df) == 0) {
    warning(
      "No ASVs passed min_rel_abund for ",
      taxon_label,
      "; plotting available ASVs instead.",
      call. = FALSE
    )
    plot_df <- copy(alluvial_df)
  }

  asv_rank <- plot_df[, .(
    ASV_abund = sum(RelativeAbundanceWithinTaxon, na.rm = TRUE)
  ), by = ASVid][order(-ASV_abund)]

  top_asvs <- head(asv_rank$ASVid, top_n_asvs)
  if (length(top_asvs) == 0) {
    warning("No ASVs available for alluvial plot: ", taxon_label, call. = FALSE)
    return(invisible(NULL))
  }

  plot_df <- plot_df[ASVid %in% top_asvs]
  asv_order <- asv_rank[ASVid %in% top_asvs]$ASVid
  plot_df[, ASVid := factor(ASVid, levels = rev(asv_order))]
  plot_df[, ASVid_short := substr(as.character(ASVid), 1, 10)]

  for (col in rank_cols) {
    plot_df[, (col) := factor(get(col))]
  }

  fill_col <- rank_cols[length(rank_cols)]
  p <- ggplot(
    plot_df,
    aes(
      axis1 = ASVid_short,
      axis2 = .data[[rank_cols[1]]],
      axis3 = .data[[rank_cols[2]]],
      y = flow_value
    )
  ) +
    geom_alluvium(
      aes(fill = .data[[fill_col]], color = .data[[fill_col]]),
      alpha = 0.85,
      width = 0.08
    ) +
    geom_stratum(width = 0.08, color = "grey45", fill = "grey92") +
    scale_x_discrete(
      limits = c(
        "ASVs",
        paste0("BLAST ", rank_label(plot_ranks[1])),
        paste0("BLAST ", rank_label(plot_ranks[2]))
      ),
      expand = c(0.05, 0.05)
    ) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    labs(
      title = paste0("BLAST support for Kraken ", taxa_level, ": ", taxon_label),
      x = NULL,
      y = paste0("Relative abundance within Kraken ", taxa_level)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 12, height = 8, dpi = 300)
  }

  p
}

process_taxon <- function(manifest_file,
                          blast_file,
                          out_plot_file,
                          out_table_file,
                          sql_db,
                          taxa_level,
                          taxon_path_name,
                          plot_ranks,
                          min_pident = 90,
                          min_qcovs = 90,
                          max_evalue = 1e-20,
                          top_n_asvs = 20,
                          min_rel_abund = 0.05) {
  manifest <- read_taxon_manifest(manifest_file, taxa_level, taxon_path_name)
  if (nrow(manifest) == 0) {
    return(invisible(NULL))
  }

  taxon_label <- unique(as.character(manifest$Taxon))
  taxon_label <- taxon_label[!is.na(taxon_label) & nzchar(taxon_label)]
  taxon_label <- if (length(taxon_label) == 0) taxon_path_name else paste(taxon_label, collapse = ", ")

  hits <- read_blast_hits(blast_file)
  if (nrow(hits) > 0) {
    hits <- filter_blast_hits(
      hits,
      min_pident = min_pident,
      min_qcovs = min_qcovs,
      max_evalue = max_evalue
    )
  }

  if (nrow(hits) == 0) {
    alluvial_df <- no_acceptable_hit_input(manifest, plot_ranks)
  } else {
    tax_map <- resolve_taxid_labels(hits$staxids, sql_db, plot_ranks)
    hits <- merge(hits, tax_map, by.x = "staxids", by.y = "taxid", all.x = TRUE)

    for (rank in plot_ranks) {
      col <- rank_col(rank)
      hits[is.na(get(col)) | !nzchar(get(col)), (col) := paste0(rank_prefix(rank), "__Unresolved")]
    }

    hits_weighted <- collapse_hits_to_weights(hits, plot_ranks)
    alluvial_df <- build_alluvial_input(manifest, hits_weighted, plot_ranks)
  }

  fwrite(alluvial_df, out_table_file, sep = "\t", quote = FALSE)
  plot_taxon_alluvial(
    alluvial_df = alluvial_df,
    taxon_label = taxon_label,
    taxa_level = taxa_level,
    plot_ranks = plot_ranks,
    out_file = out_plot_file,
    top_n_asvs = top_n_asvs,
    min_rel_abund = min_rel_abund
  )

  invisible(alluvial_df)
}

if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  blast_config <- cfg$blast_confirmation %||% list()

  trial_dir <- if (!is_missing_value(blast_config$io_dir)) {
    file.path(blast_config$io_dir, "BlastAnalysis")
  } else {
    NULL
  }
  comparisons <- blast_config$DA_comparisons
  tax_db_sql <- blast_config$tax_db_sql
  taxa_level <- blast_config$taxa_level %||% "Genus"
  plot_ranks <- blast_config$plot_ranks %||% args$plot_ranks
  min_pident <- as_config_number(blast_config$min_pident, args$min_pident)
  min_qcovs <- as_config_number(blast_config$min_qcovs, args$min_qcovs)
  max_evalue <- as_config_number(blast_config$max_evalue, args$max_evalue)
  top_n_asvs <- as_config_integer(blast_config$top_n_asvs, args$top_n_asvs)
  min_rel_abund <- as_config_number(blast_config$min_rel_abund, args$min_rel_abund)
} else {
  trial_dir <- args$trial_dir
  comparisons <- args$comparisons
  tax_db_sql <- args$tax_db_sql
  taxa_level <- args$taxa_level
  plot_ranks <- args$plot_ranks
  min_pident <- args$min_pident
  min_qcovs <- args$min_qcovs
  max_evalue <- args$max_evalue
  top_n_asvs <- args$top_n_asvs
  min_rel_abund <- args$min_rel_abund
}

plot_ranks <- validate_plot_ranks(plot_ranks)

required_input <- list(
  trial_dir = trial_dir,
  comparisons = comparisons,
  tax_db_sql = tax_db_sql,
  taxa_level = taxa_level
)
fail_missing_input(required_input, "BLAST plotting input")

if (!dir.exists(trial_dir)) {
  stop("Missing BlastAnalysis directory: ", trial_dir, call. = FALSE)
}
if (!file.exists(tax_db_sql)) {
  stop("Missing taxonomizr SQLite database: ", tax_db_sql, call. = FALSE)
}

for (comparison in comparisons) {
  comp_dir <- file.path(trial_dir, comparison)

  if (!dir.exists(comp_dir)) {
    warning("Skipping missing comparison directory: ", comp_dir, call. = FALSE)
    next
  }

  taxon_dirs <- list.dirs(comp_dir, recursive = FALSE, full.names = TRUE)

  if (length(taxon_dirs) == 0) {
    warning("No taxon directories found in: ", comp_dir, call. = FALSE)
    next
  }

  for (taxon_dir in taxon_dirs) {
    taxon_path_name <- basename(taxon_dir)
    manifest_file <- file.path(taxon_dir, paste0(taxon_path_name, "_ASV_manifest.tsv"))
    blast_file <- file.path(taxon_dir, paste0(taxon_path_name, "_blast_hits.tsv"))

    if (!file.exists(manifest_file)) {
      warning("Missing manifest: ", manifest_file, call. = FALSE)
      next
    }
    if (!file.exists(blast_file)) {
      warning("Missing BLAST hits: ", blast_file, call. = FALSE)
      next
    }

    out_plot_file <- file.path(taxon_dir, paste0(taxon_path_name, "_alluvial.png"))
    out_table_file <- file.path(taxon_dir, paste0(taxon_path_name, "_alluvial_input.tsv"))

    message("Processing ", comparison, " / ", taxon_path_name)

    process_taxon(
      manifest_file = manifest_file,
      blast_file = blast_file,
      out_plot_file = out_plot_file,
      out_table_file = out_table_file,
      sql_db = tax_db_sql,
      taxa_level = taxa_level,
      taxon_path_name = taxon_path_name,
      plot_ranks = plot_ranks,
      min_pident = min_pident,
      min_qcovs = min_qcovs,
      max_evalue = max_evalue,
      top_n_asvs = top_n_asvs,
      min_rel_abund = min_rel_abund
    )
  }
}
