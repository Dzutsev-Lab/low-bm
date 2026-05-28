library(argparse)
library(data.table)
library(Biostrings)
library(phyloseq)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(taxonomizr)
library(tibble)
library(stringr)

parser <- ArgumentParser()

parser$add_argument("--trial-dir",
                    type = "character",
                    help = "Trial directory containing comparison directories")

parser$add_argument("--comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Comparison directories to process")

parser$add_argument("--out-dir",
                    type = "character",
                    help = "Base output directory for plots and intermediate tables")

parser$add_argument("--tax-db-sql",
                    type = "character",
                    help = "taxonomizr SQLite database for resolving BLAST taxids to names")

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
                    help = "Maximum number of ASVs to show in each genus plot")

args <- parser$parse_args()



#-----------------------------
# Helpers
#-----------------------------
normalize_taxon <- function(x) {
  x <- as.character(x)
  x <- sub("^[a-zA-Z]__", "", x)
  x <- sub("^UC_", "", x)
  x[x == "" | x == "Unclassified"] <- NA_character_
  x
}

get_asv_abundance <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    rowSums(otu)
  } else {
    colSums(otu)
  }
}

read_genus_manifest <- function(manifest_file) {
  m <- fread(manifest_file)
  m[, ASVid := as.character(ASVid)]
  if ("TotalCount" %in% names(m)) {
    m[, TotalCount := as.numeric(TotalCount)]
  }
  if ("RelativeAbundanceWithinGenus" %in% names(m)) {
    m[, RelativeAbundanceWithinGenus := as.numeric(RelativeAbundanceWithinGenus)]
  }
  m
}

extract_taxid_from_sseqid <- function(x) {
  x <- as.character(x)

  m <- regexec("kraken:taxid\\|([0-9]+)$", x)
  parts <- regmatches(x, m)[[1]]

  if (length(parts) == 2) {
    return(as.integer(parts[2]))
  }

  NA_integer_
}

read_blast_hits <- function(blast_file) {
  hits <- fread(
    blast_file,
    header = FALSE,
    sep = "\t",
    col.names = c("qseqid", "sseqid", "pident", "align_length", "qcovs",
                  "evalue", "bitscore", "staxids")
  )

  hits[, `:=`(
    qseqid = as.character(qseqid),
    sseqid = as.character(sseqid),
    pident = as.numeric(pident),
    align_length = as.numeric(align_length),
    qcovs = as.numeric(qcovs),
    evalue = as.numeric(evalue),
    bitscore = as.numeric(bitscore),
    staxids = as.character(staxids)
  )]
  hits[, staxids := vapply(sseqid, extract_taxid_from_sseqid, integer(1))]

  hits
}

filter_blast_hits <- function(hits,
                              min_pident = 97,
                              min_qcovs = 90,
                              max_evalue = 1e-20) {
  hits[pident >= min_pident &
         qcovs >= min_qcovs &
         evalue <= max_evalue]
}

# TODO: determine function
parse_first_taxid <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(trimws(x)) || x == "0") return(NA_integer_)

  vals <- unlist(strsplit(x, ";", fixed = TRUE))
  vals <- suppressWarnings(as.integer(trimws(vals)))
  vals <- vals[!is.na(vals) & vals > 0L]

  if (length(vals) == 0) return(NA_integer_)
  vals[1]
}




resolve_taxid_labels <- function(taxids, sql_db) {
    taxids <- unique(as.integer(taxids))
    taxids <- taxids[!is.na(taxids) & taxids > 0L]

    if (length(taxids) == 0) {
        return(data.table(
            taxid = integer(0), 
            BlastGenus = character(0),
            BlastSpecies = character(0)
        ))
    }

    desired_taxa <- c("superkingdom", "phylum", "class", "order",
                        "family", "genus", "species")

    tax_df <- getTaxonomy(ids = taxids, sqlFile = sql_db, desiredTaxa = desired_taxa)
    tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
    tax_df$taxid <- taxids

    # choose the deepest non-NA rank available
    rank_order <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")
    rank_prefix <- c(
        superkingdom = "d",
        phylum = "p",
        class = "c",
        order = "o",
        family = "f",
        genus = "g",
        species = "s"
    )

    tax_df$BlastGenus <- ifelse(
        !is.na(tax_df$genus) & nzchar(tax_df$genus),
        paste0("g__", tax_df$genus),
        "g__Unresolved"
    )

    tax_df$BlastSpecies <- ifelse(
        !is.na(tax_df$species) & nzchar(tax_df$species),
        paste0("s__", tax_df$species),
        "s__Unresolved"
    )

    as.data.table(tax_df[, c("taxid", "BlastGenus", "BlastSpecies")])
}

collapse_hits_to_weights <- function(hits) {
    if (nrow(hits) == 0) {
        return(data.table(
            qseqid = character(0),
            BlastGenus = character(0),
            BlastSpecies = character(0),
            hit_share = numeric(0)
        ))
    }

    hits <- hits[, .(
        bitscore = sum(bitscore),
        pident = max(pident, na.rm = TRUE),
        qcovs = max(qcovs, na.rm = TRUE),
        evalue = min(evalue, na.rm = TRUE)
    ), by = .(qseqid, BlastGenus, BlastSpecies)]

    hits[, hit_share := bitscore / sum(bitscore), by = qseqid]
    hits
}

build_alluvial_input <- function(manifest, hits_weighted) {
    out <- merge(
        manifest[, .(ASVid, RelativeAbundanceWithinGenus)],
        hits_weighted[, .(qseqid, BlastGenus, BlastSpecies, hit_share)],
        by.x = "ASVid",
        by.y = "qseqid",
        all.x = TRUE
    )

    out[is.na(BlastGenus), BlastGenus := "g__No acceptable hit"]
    out[is.na(BlastSpecies), BlastSpecies := "s__No acceptable hit"]
    out[is.na(hit_share), hit_share := 1]

    out[, flow_value := RelativeAbundanceWithinGenus * hit_share]
    out
}


plot_genus_alluvial <- function(alluvial_df, 
                                genus_name, 
                                out_file = NULL, 
                                top_n_asvs = 20,
                                min_rel_abund = 0.05) {
    # Keep only ASVs above the abundance threshold
    plot_df <- alluvial_df[RelativeAbundanceWithinGenus >= min_rel_abund]

    # Keep only top ASVs by abundance for readability (after filtering)
    asv_rank <- plot_df[, .(
        ASV_abund = sum(RelativeAbundanceWithinGenus, na.rm = TRUE)
    ), by = ASVid][order(-ASV_abund)]

    top_asvs <- head(asv_rank$ASVid, top_n_asvs)
    plot_df <- plot_df[ASVid %in% top_asvs]

    # Order ASVs by abundance
    asv_order <- asv_rank[ASVid %in% top_asvs]$ASVid
    plot_df[, ASVid := factor(ASVid, levels = rev(asv_order))]

    # Short display labels
    plot_df[, ASVid_short := substr(ASVid, 1, 10)]
    plot_df[, BlastGenus := factor(BlastGenus)]
    plot_df[, BlastSpecies := factor(BlastSpecies)]

    p <- ggplot(plot_df,
                aes(axis1 = ASVid_short, 
                    axis2 = BlastGenus, 
                    axis3 = BlastSpecies, 
                    y = flow_value)) +
        geom_alluvium(aes(fill = BlastSpecies, color = BlastSpecies), 
                      alpha = 0.85, 
                      width = 0.08) +
        geom_stratum(width = 0.08, 
                     color = "grey45", 
                     fill = "grey92") +
        scale_x_discrete(
            limits = c("ASVs", "BLAST Genus", "BLAST Species"),
            expand = c(0.05, 0.05)
        ) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
        labs(
            title = paste0("BLAST support Kraken Assignment: ", genus_name),
            x = NULL,
            y = "Relative abundance within Kraken genus"
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

process_genus <- function(manifest_file,
                          blast_file,
                          out_plot_file,
                          out_table_file,
                          sql_db,
                          min_pident = 97,
                          min_qcovs = 90,
                          max_evalue = 1e-20,
                          top_n_asvs = 20,
                          min_rel_abund = 0.05) {
    manifest <- read_genus_manifest(manifest_file)
    if (nrow(manifest) == 0) return(invisible(NULL))

    hits <- read_blast_hits(blast_file)
    if (nrow(hits) == 0) {
        out_dt <- manifest[, .(ASVid, RelativeAbundanceWithinGenus)]
        out_dt[, `:=`(BlastGenus = "No acceptable genus hit", BlastSpecies = "No acceptable species hit", hit_share = 1, flow_value = RelativeAbundanceWithinGenus)]
        fwrite(out_dt, out_table_file, sep = "\t", quote = FALSE)
        plot_genus_alluvial(out_dt, unique(manifest$Genus), out_plot_file, top_n_asvs)
        return(invisible(out_dt))
    }

    hits <- filter_blast_hits(hits,
                                min_pident = min_pident,
                                min_qcovs = min_qcovs,
                                max_evalue = max_evalue)

    if (nrow(hits) == 0) {
        out_dt <- manifest[, .(ASVid, RelativeAbundanceWithinGenus)]
        out_dt[, `:=`(
            BlastGenus = "g__No acceptable hit",
            BlastSpecies = "s__No acceptable hit",
            hit_share = 1,
            flow_value = RelativeAbundanceWithinGenus
        )]
        fwrite(out_dt, out_table_file, sep = "\t", quote = FALSE)
        plot_genus_alluvial(
            alluvial_df = out_dt,
            genus_name = unique(manifest$Genus),
            out_file = out_plot_file,
            top_n_asvs = top_n_asvs,
            min_rel_abund = min_rel_abund
        )
        return(invisible(out_dt))
    }

    tax_map <- resolve_taxid_labels(hits$staxids, sql_db)
    hits <- merge(hits, tax_map, by.x = "staxids", by.y = "taxid", all.x = TRUE)
    
    hits_weighted <- collapse_hits_to_weights(hits)
    alluvial_df <- build_alluvial_input(manifest, hits_weighted)

    fwrite(alluvial_df, out_table_file, sep = "\t", quote = FALSE)
    plot_genus_alluvial(
        alluvial_df = alluvial_df,
        genus_name = unique(manifest$Genus),
        out_file = out_plot_file,
        top_n_asvs = top_n_asvs,
        min_rel_abund = min_rel_abund
    )

    invisible(alluvial_df)
}


#-----------------------------
# Main directory walk
#-----------------------------
for (comparison in args$comparisons) {
  comp_dir <- file.path(args$trial_dir, comparison)

  if (!dir.exists(comp_dir)) {
    warning("Skipping missing comparison directory: ", comp_dir)
    next
  }

  genus_dirs <- list.dirs(comp_dir, recursive = FALSE, full.names = TRUE)

  if (length(genus_dirs) == 0) {
    warning("No genus directories found in: ", comp_dir)
    next
  }

  for (genus_dir in genus_dirs) {
    genus_name <- basename(genus_dir)

    manifest_file <- file.path(genus_dir, paste0(genus_name, "_ASV_manifest.tsv"))
    blast_file <- file.path(genus_dir, paste0(genus_name, "_blast_hits.tsv"))

    if (!file.exists(manifest_file)) {
      warning("Missing manifest: ", manifest_file)
      next
    }
    if (!file.exists(blast_file)) {
      warning("Missing BLAST hits: ", blast_file)
      next
    }

    out_plot_file <- file.path(genus_dir, paste0(genus_name, "_alluvial.png"))
    out_table_file <- file.path(genus_dir, paste0(genus_name, "_alluvial_input.tsv"))

    message("Processing ", comparison, " / ", genus_name)

    process_genus(
      manifest_file = manifest_file,
      blast_file = blast_file,
      out_plot_file = out_plot_file,
      out_table_file = out_table_file,
      sql_db = args$tax_db_sql,
      min_pident = args$min_pident,
      min_qcovs = args$min_qcovs,
      max_evalue = args$max_evalue,
      top_n_asvs = args$top_n_asvs
    )
  }
}