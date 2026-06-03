library(phyloseq)
library(Biostrings)
library(data.table)
library(tibble)
library(stringr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--fasta-files",
                    type = "character",
                    nargs = "+",
                    help = "List of fasta files to be merged")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "Output directory for merged FASTA file")

args <- parser$parse_args()

seqs_list <- lapply(args$fasta_files, readDNAStringSet)
merged_seqs <- do.call(c, seqs_list)

if (length(merged_seqs) == 0) {
    stop("No sequences found in the input FASTA files")
}

asv_dt <- data.table(
  ASV_ID = names(merged_seqs),
  Sequence = as.character(merged_seqs)
)

#---------------------------------
# Deduplicate by ASV ID
# - keep identical duplicates
# - stop if same ASV ID has different sequences
#---------------------------------
dup_check <- asv_dt[, .(
  n_rows = .N,
  n_unique_seqs = uniqueN(Sequence)
), by = ASV_ID][n_rows > 1]

if (nrow(dup_check) > 0) {
  bad_ids <- dup_check[n_unique_seqs > 1, ASV_ID]

  if (length(bad_ids) > 0) {
    stop(
      "Found duplicated ASV IDs with different sequences: ",
      paste(bad_ids, collapse = ", ")
    )
  }

  message("Found ", nrow(dup_check), " duplicated ASV IDs with identical sequences; keeping one copy.")
}

dt_unique <- asv_dt[!duplicated(ASV_ID)]

#---------------------------------
# Write merged FASTA
#---------------------------------
merged_unique <- DNAStringSet(dt_unique$Sequence)
names(merged_unique) <- dt_unique$ASV_ID

writeXStringSet(merged_unique, filepath = file.path(args$out_dir, "MergedASV.fasta"))

message("Wrote merged FASTA with ", length(merged_unique), " unique ASVs to: ", args$out_dir, "/MergedASV.fasta")
