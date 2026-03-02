library(ShortRead)
library(dplyr)
library(tibble)
library(argparse)

#---------------------------------
# Execution Arguments
#---------------------------------
parser <- ArgumentParser()

#
parser$add_argument("--sample-name-file", type="character", help="Text file containing list of sample names, one per line")
parser$add_argument("--raw-dir", type="character", help="File path to directory with raw reads")
parser$add_argument("--selected-dir", type="character", help="File path to directory with selected authentic reads")

parser$add_argument("--dada-filter-counts", type="character", help="File path to tsv with Dada filtering stage counts")
parser$add_argument("--seq-table", type="character", help="File path to tsv of denoised sequence table")

parser$add_argument("--host-names", type="character", help="File path to list of ASV's unmapped to host reference")
parser$add_argument("--viral-names", type="character", help="File path to list of ASV's unmapped to viral reference")
parser$add_argument("--bacterial-names", type="character", help="File path to list of ASV's mapped to bacterial reference")

parser$add_argument("--combined-counts", type="character", help="Output combined read count tsv file path")

args <- parser$parse_args()

sample_names <- readLines(args$sample_name_file)
sample_names <- trimws(sample_names)
sample_names <- sample_names[sample_names != ""]

#----------------------------------------
# Read Counts from UMI Selection Summary
#----------------------------------------
# - Raw reads
# - Selected authentic reads

# Function to read FastQ read counts from UMI selection summary files
read_umi_selection_summary <- function(sample_name, selected_dir) {
  summary_file <- file.path(selected_dir, paste0("CountSummary.", sample_name, ".tsv"))
  summary_df <- read.delim(summary_file, header = TRUE, stringsAsFactors = FALSE)
  row <- summary_df[1, ]
  tibble(
    SampleID = as.character(row$SampleID),
    Raw_reads = row$Raw_reads,
    Selected_reads = row$Selected_reads
  )
}

# Read and combine UMI selection summaries for all samples
umi_selection_summary_counts_list <- lapply(sample_names, read_umi_selection_summary, selected_dir = args$selected_dir)
umi_selection_summary_counts <- bind_rows(umi_selection_summary_counts_list)
str(umi_selection_summary_counts)


#---------------------------------
# Read Counts from Dada Denoising
#---------------------------------
dada_counts <- read.delim(args$dada_filter_counts, header =TRUE)
str(dada_counts)

#---------------------------------
# Filter Read Counts from ASV
#---------------------------------
seq_table <- read.delim(args$seq_table, header = TRUE)
seq_table <- rename(seq_table, 'SampleID' = 'X')
str(seq_table)

make_ASV_summary_df <- function(tag, name_file, seq_table) {
  ASV_names <- readLines(name_file)
  ASV_cols <- intersect(colnames(seq_table), ASV_names)
  
  tibble(
    SampleID = seq_table$SampleID,
    !!tag := if (length(ASV_cols) == 0) {
      0 
    } else {
      rowSums(seq_table[, ASV_cols, drop = FALSE])
    }
  )
}
host_unmapped_summary_df    <- make_ASV_summary_df('HostUnmapped_reads', args$host_names, seq_table)
viral_unmapped_summary_df   <- make_ASV_summary_df('ViralUnmapped_reads', args$viral_names, seq_table)
bacterial_mapped_summary_df <- make_ASV_summary_df('BacterialMapped_reads', args$bacterial_names, seq_table)

str(bacterial_mapped_summary_df)

#---------------------------------
# Combine All Read Counts
#---------------------------------
all_summary_tables <- list(
  umi_selection_summary_counts,
  dada_counts, 
  host_unmapped_summary_df, 
  viral_unmapped_summary_df,
  bacterial_mapped_summary_df 
)

combined_counts <- Reduce(
  function(x, y) full_join(x, y, by = "SampleID"), all_summary_tables
)

write.table(combined_counts, 
            file = args$combined_counts, 
            sep = '\t', 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

