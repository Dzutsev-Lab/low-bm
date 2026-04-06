library(dplyr)
library(tibble)
library(Biostrings)

parser <- ArgumentParser()

parser$add_argument("--seq-tables",
                    type = "character",
                    nargs='+',
                    help = "List of run-specific sequence tables to be concatenated")
parser$add_argument("--out",
                    type = "character",
                    help = "Desired output directory")

args <- parser$parse_args()

seq_tables <- lapply(args$seq_tables, function(f) {
    read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
})

combine_seq_tables <- function(seq_table_list) {
    seq_table_list <- lapply(seq_table_list, function(x) {
        rownames_to_column(x, var = "sampleID")
    })

    out <- bind_rows(seq_table_list)
    column_to_rownames(out, var = "sampleID")
    out[is.na(out)] <- 0

    out
}

combined_seq_table <- combine_seq_tables(seq_tables)
# Constructing stable ASV IDs
ASV_IDs <- paste0("ASV", seq_len(ncol(combined_seq_table)))

#-------------------------
# Construct FASTA 
#-------------------------
dna_strings <- DNAStringSet(colnames(combined_seq_table))
names(dna_strings) <- ASV_IDs
writeXStringSet(dna_strings, 
                filepath = paste0(args$out, "/ASV.fasta"), 
                format = "fasta")

colnames(combined_seq_table) <- ASV_IDs
write.table(combined_seq_table,
            file = paste0(args$out, "/SeqTable.tsv"), 
            sep = '\t', 
            quote = FALSE, 
            col.names = NA)

