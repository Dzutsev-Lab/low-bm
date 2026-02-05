#----------------------------
# Log Construction
#----------------------------  
logfile <- snakemake@log[[1]]
if (!is.null(logfile) && nzchar(logfile)) {
  logcon <- file(logfile, open = "wt")
  sink(logcon, type = "output")
  sink(logcon, type = "message")
  cat("----- R script started -----\n")
  cat("Working dir:", getwd(), "\n")
  cat("Sys.getenv PATH:", Sys.getenv("PATH"), "\n")
  flush.console()
}


library(ShortRead)
library(dplyr)
library(tibble)

#---------------------------
# Read Counts from fastq's
#---------------------------
sample_names <- snakemake@params$samples

fastq_read_counter <- function(fq_path) {
  wc_out <- system2("wc", c("-l", fq_path), stdout = TRUE, stderr = "")
  line_count <- as.numeric(strsplit(wc_out[1], "\\s+")[[1]][1])
  read_count <- line_count / 4
  read_count
}

fastq_counts <- lapply(sample_names, function(s){
  normalized_path <- file.path(snakemake@params$normalized, paste0(s, "_R1_001.fastq"))
  selected_path   <- file.path(snakemake@params$selected, paste0("Selected.", s, ".UMI_R1.fastq"))
  dedpued_path    <- file.path(snakemake@params$deduped, paste0("Deduped.", s, ".fastq"))
  tibble(
    SampleID = s,
    Raw_reads      = fastq_read_counter(normalized_path),
    Selected_reads = fastq_read_counter(selected_path),
    Deduped_reads  = fastq_read_counter(dedpued_path)
    )
}) |> bind_rows()
str(fastq_counts)

#---------------------------------
# Read Counts from Dada Denoising
#---------------------------------
dada_counts <- read.delim(snakemake@input$dada_read_counts, header =TRUE)
str(dada_counts)

#---------------------------------
# Filter Read Counts from ASV
#---------------------------------
seq_table <- read.delim(snakemake@input$seq_table, header = TRUE)
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
host_unmapped_summary_df    <- make_ASV_summary_df('HostUnmapped_reads', snakemake@input$host_unmapped_names, seq_table)
viral_unmapped_summary_df   <- make_ASV_summary_df('ViralUnmapped_reads', snakemake@input$viral_unmapped_names, seq_table)
bacterial_mapped_summary_df <- make_ASV_summary_df('BacterialMapped_reads', snakemake@input$bacterial_names, seq_table)

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
            file = snakemake@output$combined_read_counts, 
            sep = '\t', 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}
