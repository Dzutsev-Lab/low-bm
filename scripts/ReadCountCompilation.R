# #----------------------------
# # Log Construction
# #----------------------------  
# logfile <- snakemake@log[[1]]
# if (!is.null(logfile) && nzchar(logfile)) {
#   logcon <- file(logfile, open = "wt")
#   sink(logcon, type = "output")
#   sink(logcon, type = "message")
#   cat("----- R script started -----\n")
#   cat("Working dir:", getwd(), "\n")
#   cat("Sys.getenv PATH:", Sys.getenv("PATH"), "\n")
#   flush.console()
# }


library(ShortRead)
library(dplyr)
library(tibble)
library(argparse)

#---------------------------------
# Execution Arguments
#---------------------------------
parser <- ArgumentParser()

#
parser$add_argument("--sample-names", type="character", nargs='+', help="List of sample names")
parser$add_argument("--raw-dir", type="character", help="File path to directory with raw reads")
parser$add_argument("--selected-dir", type="character", help="File path to directory with selected authentic reads")

parser$add_argument("--dada-filter-counts", type="character", help="File path to tsv with Dada filtering stage counts")
parser$add_argument("--seq-table", type="character", help="File path to tsv of denoised sequence table")

parser$add_argument("--host-names", type="character", help="File path to list of ASV's unmapped to host reference")
parser$add_argument("--viral-names", type="character", help="File path to list of ASV's unmapped to viral reference")
parser$add_argument("--bacterial-names", type="character", help="File path to list of ASV's mapped to bacterial reference")

parser$add_argument("--combined-counts", type="character", help="Output combined read count tsv file path")

args <- parser$parse_args()

sample_names <- args$sample_names


fastq_read_counter <- function(fq_path) {
  wc_out <- system2("wc", c("-l", fq_path), stdout = TRUE, stderr = "")
  line_count <- as.numeric(strsplit(wc_out[1], "\\s+")[[1]][1])
  read_count <- line_count / 4
  read_count
}

fastq_counts <- lapply(sample_names, function(s){
  normalized_path <- file.path(args$raw_dir, paste0(s, "_R1_001.fastq"))
  selected_path   <- file.path(args$selected_dir, paste0("Selected.", s, ".UMI_R1.fastq"))
  tibble(
    SampleID = s,
    Raw_reads      = fastq_read_counter(normalized_path),
    Selected_reads = fastq_read_counter(selected_path)
    )
}) |> bind_rows()
str(fastq_counts)

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
  fastq_counts,
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

# #----------------------------
# # Log Close
# #----------------------------
# if (!is.null(logfile) && nzchar(logfile)) {
#   sink(type = "message"); sink(type = "output")
#   close(logcon)
# }
