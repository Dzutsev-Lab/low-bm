library(data.table)
library(argparse)

source("scripts/Rhelpers/SecondaryClassificationHelpers.R")

parser <- ArgumentParser()

parser$add_argument("--kraken-file",
                    type = "character",
                    help = "Kraken output file with per-ASV classifications")
parser$add_argument("--blast-hits",
                    type = "character",
                    help = "blast_hits.tsv containing best BLAST hit per ASV")
parser$add_argument("--tax-db-nodes",
                    type = "character",
                    help = "nodes.dmp file from the Kraken database")
parser$add_argument("--reconciled-out",
                    type = "character",
                    help = "Output Kraken-format file with reconciled taxonomy")
parser$add_argument("--summary-out",
                    type = "character",
                    help = "Output summary table of Kraken vs BLAST reconciliation")
parser$add_argument("--conflict-policy",
                    type = "character",
                    default = "lca",
                    choices = c("lca", "keep_kraken", "keep_blast"),
                    help = "How to resolve Kraken/BLAST conflicts when neither is ancestor of the other")

args <- parser$parse_args()

#---------------------------
# Read nodes.dmp
#---------------------------
nodes_df <- fread(
  args$tax_db_nodes,
  sep = "|",
  header = FALSE,
  fill = TRUE,
  quote = "",
  data.table = TRUE
)

nodes_df <- nodes_df[, .(
  TaxID = as.integer(trimws(V1)),
  Parent = as.integer(trimws(V2)),
  rank = trimws(V3)
)]

nodes_df <- nodes_df[!is.na(TaxID) & !is.na(Parent)]

rank_lookup <- setNames(nodes_df$rank, as.character(nodes_df$TaxID))
parent_lookup <- setNames(nodes_df$Parent, as.character(nodes_df$TaxID))

#---------------------------
# Read Kraken output
#---------------------------
kraken_df <- fread(
  args$kraken_file,
  header = FALSE,
  sep = "\t",
  col.names = c("status", "ASVid", "TaxID", "length", "kmer_trace")
)

kraken_df[, TaxID := as.integer(TaxID)]
kraken_df[, kraken_taxid := TaxID]

#---------------------------
# Read BLAST best hits
#---------------------------
blast_df <- fread(
  args$blast_hits,
  header = FALSE,
  sep = "\t",
  #TODO: fix to match actual blast hits file
  col.names = c("qseqid", "sseqid", "pident", "align_length", "qcovs",
                "evalue", "bitscore", "staxids", "sscinames")
)

setnames(blast_df, "qseqid", "ASVid")

# Make sure character fields stay character
blast_df[, `:=`(
  sseqid = as.character(sseqid)
)]

# Keep only the best hit per ASV
setorder(blast_df, ASVid, -bitscore, -pident, -qcovs, evalue)
blast_df <- blast_df[, .SD[1], by = ASVid]

# Safety check
stopifnot(!anyDuplicated(blast_df$ASVid))

#-----------------------------
# Formatting Helper Functions
#-----------------------------
extract_taxid_from_sseqid <- function(sseqid) {
  sseqid <- as.character(sseqid)

  m <- regexec("kraken:taxid\\|([0-9]+)$", sseqid)
  parts <- regmatches(sseqid, m)[[1]]

  if (length(parts) == 2) {
    return(as.integer(parts[2]))
  }

  0L
}

parse_staxids <- function(x) {
  x <- as.character(x)

  if (is.na(x) || !nzchar(trimws(x))) {
    return(integer(0))
  }

  vals <- unlist(strsplit(x, ";", fixed = TRUE))
  vals <- suppressWarnings(as.integer(trimws(vals)))
  vals <- vals[!is.na(vals) & vals > 0L]
  unique(vals)
}

collapse_blast_taxid <- function(staxids_string, parent_lookup) {
  taxids <- parse_staxids(staxids_string)
  if (length(taxids) == 0) return(0L)
  if (length(taxids) == 1) return(taxids[1])

  cur <- taxids[1]
  for (tx in taxids[-1]) {
    cur <- lowest_common_ancestor(cur, tx, parent_lookup)
    if (cur == 0L) break
  }
  cur
}


#---------------------------------
# Primary Reconciliation Function
#---------------------------------
reconcile_taxids <- function(kraken_taxid, blast_taxid, parent_lookup, conflict_policy = "lca") {
  kraken_taxid <- as.integer(kraken_taxid)
  blast_taxid <- as.integer(blast_taxid)

  if (is.na(kraken_taxid)) kraken_taxid <- 0L
  if (is.na(blast_taxid)) blast_taxid <- 0L

  if (blast_taxid == 0L) {
    return(list(final_taxid = kraken_taxid, action = "keep_kraken"))
  }

  if (kraken_taxid == 0L) {
    return(list(final_taxid = blast_taxid, action = "adopt_blast"))
  }

  if (kraken_taxid == blast_taxid) {
    return(list(final_taxid = kraken_taxid, action = "agree"))
  }

  if (is_ancestor_of(kraken_taxid, blast_taxid, parent_lookup)) {
    return(list(final_taxid = blast_taxid, action = "blast_more_specific"))
  }

  if (is_ancestor_of(blast_taxid, kraken_taxid, parent_lookup)) {
    return(list(final_taxid = kraken_taxid, action = "kraken_more_specific"))
  }

  if (conflict_policy == "keep_kraken") {
    return(list(final_taxid = kraken_taxid, action = "conflict_keep_kraken"))
  }

  if (conflict_policy == "keep_blast") {
    return(list(final_taxid = blast_taxid, action = "conflict_keep_blast"))
  }

  anc <- lowest_common_ancestor(kraken_taxid, blast_taxid, parent_lookup)
  if (is.na(anc) || anc %in% c(0L, 1L)) {
    return(list(final_taxid = 0L, action = "conflict_unclassified"))
  }

  list(final_taxid = anc, action = "conflict_lca")
}

#---------------------------
# Prepare BLAST taxonomy IDs
#---------------------------
blast_df[, blast_taxid := as.integer(mapply(
  extract_taxid_from_sseqid,
  sseqid = sseqid
))]


#---------------------------
# Merge Kraken and BLAST
#---------------------------
kraken_blast_comp_df <- merge(
  kraken_df,
  blast_df,
  by = "ASVid",
  all.x = TRUE
)

#---------------------------
# Reconcile taxonomy
#---------------------------
kraken_blast_reconciliation <- mapply(
  FUN = reconcile_taxids,
  kraken_taxid = kraken_blast_comp_df$kraken_taxid,
  blast_taxid = kraken_blast_comp_df$blast_taxid,
  MoreArgs = list(parent_lookup = parent_lookup, conflict_policy = args$conflict_policy),
  SIMPLIFY = FALSE
)

kraken_blast_comp_df[, final_taxid := as.integer(vapply(kraken_blast_reconciliation, `[[`, integer(1), "final_taxid"))]
kraken_blast_comp_df[, reconcile_action := vapply(kraken_blast_reconciliation, `[[`, character(1), "action")]

# Final Kraken-style status
kraken_blast_comp_df[, final_status := fifelse(final_taxid == 0L, "U", "C")]

cat("Kraken duplicates: ", anyDuplicated(kraken_df$ASVid), "\n")
cat("BLAST duplicates: ", anyDuplicated(blast_df$ASVid), "\n")
cat("Merged duplicates: ", anyDuplicated(kraken_blast_comp_df$ASVid), "\n")
cat("Kraken row count: ", nrow(kraken_df), "\n")
cat("BLAST row count: ", nrow(blast_df), "\n")
cat("Merged row count: ", nrow(kraken_blast_comp_df), "\n")


#---------------------------
# Write final Kraken-format output
#---------------------------
reconciled_kraken <- kraken_blast_comp_df[, .(
  status = final_status,
  ASVid,
  TaxID = final_taxid,
  length,
  kmer_trace
)]

fwrite(
  reconciled_kraken,
  file = args$reconciled_out,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

#---------------------------
# Write reconciliation summary
#---------------------------
summary_tbl <- kraken_blast_comp_df[, .(
  ASVid,
  kraken_taxid,
  blast_taxid,
  final_taxid,
  final_status,
  reconcile_action,
  pident,
  qcovs,
  bitscore,
  sscinames
)]

fwrite(
  summary_tbl,
  file = args$summary_out,
  sep = "\t",
  quote = FALSE
)