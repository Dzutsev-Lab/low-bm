library(Biostrings)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--fasta-in",
                    type = "character",
                    help = "ID number to attach to output files")
parser$add_argument("--map-out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()

asv_fasta <- readDNAStringSet(args$fasta_in)

asv_map <- data.frame(
  ASV_ID = names(asv_fasta),
  Sequence = as.character(asv_fasta),
  SeqLength = width(asv_fasta),
  stringsAsFactors = FALSE
)

write.table(
  asv_map,
  file = args$map_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)