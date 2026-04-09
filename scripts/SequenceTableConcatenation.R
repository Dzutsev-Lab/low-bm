library(dplyr)
library(tibble)
library(data.table)
library(Biostrings)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--seq-tables",
                    type = "character",
                    nargs='+',
                    help = "List of run-specific sequence tables to be concatenated")
parser$add_argument("--out",
                    type = "character",
                    help = "Desired output directory")

args <- parser$parse_args()

# Read in sequence tables and prepare them for concatenation
read_seq_table <- function(seq_table_path) {
    # Using fread over read.table for speed, but it returns a data.table.
    # Convert it to a data.frame and set the rownames to fit needed data structure.
    seq_table_dt <- fread(seq_table_path)
    sampleID <- seq_table_dt[[1]]
    seq_table_dt <- seq_table_dt[, -1, with = FALSE]
    seq_table_df <- as.data.frame(seq_table_dt)
    # Remove *_S## suffix from sample IDs to ensure consistent sample nameing with metadata.
    sampleID <- sub("_S\\d+$", "", sampleID)
    rownames(seq_table_df) <- sampleID
    seq_table_df
}

seq_tables <- lapply(args$seq_tables, read_seq_table)

# Collect all unique ASVs across all sequence tables in order and create stable ASV ID mapping.
# Done to increase speed of table merging and to ensure consistent ASV IDs between sequences.
all_ASVs <- sort(unique(unlist(lapply(seq_tables, colnames), use.names = FALSE)))

ASV_map <- setNames(sprintf("ASV%d", seq_along(all_ASVs)), all_ASVs)

# Update column names of each sequence table to use the stable ASV IDs.
seq_tables <- lapply(seq_tables, function(seq_table_df) {
    colnames(seq_table_df) <- unname(ASV_map[colnames(seq_table_df)])
    seq_table_df <- rownames_to_column(seq_table_df, var = "SampleID")
    seq_table_df
})

# Bind rows of all sequence tables together, filling in missing ASVs with zeros.
combined_seq_table <- rbindlist(seq_tables, 
                                use.names = TRUE, 
                                fill = TRUE)
combined_seq_table <- column_to_rownames(combined_seq_table, var = "SampleID")
combined_seq_table[is.na(combined_seq_table)] <- 0


#-------------------------
# Construct FASTA 
#-------------------------
# names: ASV IDs 
# sequences: nucleotide strings
dna_strings <- DNAStringSet(names(ASV_map))
names(dna_strings) <- unname(ASV_map)
writeXStringSet(dna_strings, 
                filepath = paste0(args$out, "/ASV.fasta"), 
                format = "fasta")


write.table(combined_seq_table,
            file = paste0(args$out, "/SeqTable.tsv"), 
            sep = '\t', 
            quote = FALSE, 
            col.names = NA)

