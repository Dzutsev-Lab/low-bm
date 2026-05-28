library(Biostrings)
library(data.table)
library(argparse)

source("scripts/Rhelpers/SecondaryClassificationHelpers.R")

parser <- ArgumentParser()

parser$add_argument("--kraken-file",
                    type = "character",
                    help = "Path to kraken output file (.kraken2) with per ASV classification")
parser$add_argument("--tax-db-nodes",
                    type = "character",
                    help = "Path to nodes.dmp file created as a part of kraken database construction")
parser$add_argument("--conf-cutoff",
                    type = "double",
                    default = 0.2,
                    help = "Confidence threshold used to flag genus and species labels for confirmatory BLAST")
parser$add_argument("--asv-fasta",
                    type = "character",
                    help = "Path to fasta file with ASV sequences paired with ASV ids")
parser$add_argument("--candidate-fasta",
                    type = "character",
                    help = "Path to fasta file with candidate ASV sequences for BLAST classification")

args <- parser$parse_args()

#---------------------------
# Import Kraken Results
#---------------------------
kraken_df <- fread(
    args$kraken_file,
    header = FALSE,
    sep = "\t",
    col.names = c("status", "ASVid", "TaxID", "length", "kmer_trace")
)
kraken_df[, TaxID := as.integer(TaxID)]

#---------------------------
# Import Database Nodes
#---------------------------
nodes_df <- read_kraken_nodes(args$tax_db_nodes)

#-------------------------------------
# Merge Official Tax ID Ranks to ASVs
#-------------------------------------
kraken_df <- merge(
    kraken_df,
    nodes_df[, .(TaxID, rank)],
    by = "TaxID",
    all.x = TRUE
)

#---------------------------------------------------
# Calculate Confidence Scores for Kraken Taxa Labels
#---------------------------------------------------
parent_lookup <- build_parent_lookup(nodes_df)
kraken_df[status == "C" & rank %chin% c("genus", "species", "subspecies"), conf_score := mapply(
    assigned_label_score,
    kmer_trace = kmer_trace,
    assigned_TaxID = TaxID,
    MoreArgs = list(parent_lookup = parent_lookup)
    )]


#---------------------------
# Identify BLAST Candidates
#---------------------------
# - all unclassfied
# - anything classified above genus
# - genus/species/subspecies calls below confidence threshold
unclassified_candidates <- kraken_df[status == "U" | TaxID == 0]

above_genus_candidates <- kraken_df[
  status == "C" &
    !is.na(rank) &
    !(rank %chin% c("genus", "species", "subspecies"))
]

low_conf_candidates <- kraken_df[
  status == "C" &
    rank %chin% c("genus", "species", "subspecies") &
    !is.na(conf_score) &
    conf_score < args$conf_cutoff
]

blast_candidates <- unique(rbind(
  unclassified_candidates,
  above_genus_candidates,
  low_conf_candidates,
  fill = TRUE
))

blast_candidate_IDs <- blast_candidates$ASVid

#---------------------------
# Filter ASV FASTA
#---------------------------
asv_fasta <- readDNAStringSet(args$asv_fasta)
keep <- names(asv_fasta) %in% blast_candidate_IDs
candidate_fasta <- asv_fasta[keep]

# Report any IDs not found
missing_IDs <- setdiff(blast_candidate_IDs, names(asv_fasta))
if (length(missing_IDs) > 0) {
  warning("These candidate IDs were not found in the ASV FASTA: ",
          paste(missing_IDs, collapse = ", "))
}

# Write filtered FASTA
writeXStringSet(candidate_fasta, filepath = args$candidate_fasta)