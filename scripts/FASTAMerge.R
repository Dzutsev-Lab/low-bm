library(phyloseq)
library(Biostrings)
library(data.table)
library(tibble)
library(stringr)
library(argparse)

source(file.path("scripts", "Rhelpers", "ASVFasta.R"))

parser <- ArgumentParser()

parser$add_argument("--fasta-files",
                    type = "character",
                    nargs = "+",
                    help = "List of fasta files to be merged")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "Output directory for merged FASTA file")

args <- parser$parse_args()

#---------------------------------
# Write merged FASTA
#---------------------------------
merged_unique <- read_merged_asv_fasta(args$fasta_files)
writeXStringSet(merged_unique, filepath = file.path(args$out_dir, "MergedASV.fasta"))

message("Wrote merged FASTA with ", length(merged_unique), " unique ASVs to: ", args$out_dir, "/MergedASV.fasta")
