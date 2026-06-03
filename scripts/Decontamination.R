library(micRoclean)
library(readxl)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(argparse)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))

parser <- ArgumentParser()

parser$add_argument("--seq-table",
                    type = "character",
                    help = "sequence table from Dada2 (Samples x ASV)")
parser$add_argument("--bacterial-names",
                    type = "character",
                    help = "ASV IDs associated with bacterial reads according to mapping filters")
parser$add_argument("--sample-names", 
                    type="character", 
                    help="Text file containing list of sample names, one per line")
parser$add_argument("--metadata",
                    type = "character",
                    help = "standardized metadata sheet as .xslx file")
parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number for the trial")
parser$add_argument("--out",
                    type = "character",
                    help = "output directory")

args <- parser$parse_args()

sample_names <- readLines(args$sample_names)
sample_names <- trimws(sample_names)
sample_names <- sample_names[sample_names != ""]

#---------------------------
# Construct count matrix
#---------------------------
raw_seq_table <- read.delim(args$seq_table,
                            header = TRUE,
                            row.names = 1)

bacterial_ASVids <- readLines(args$bacterial_names)
bacterial_ASVids <- bacterial_ASVids[nzchar(bacterial_ASVids)]

counts_df <- raw_seq_table[, bacterial_ASVids, drop = FALSE]

#---------------------------
# Construct sample metadata
#---------------------------
metadata_df <- read_excel(args$metadata)
metadata_df <- metadata_df |> 
  filter(SampleName %in% sample_names)

if (length(setdiff(sample_names, metadata_df$SampleName)) > 0) {
  stop(
    paste("No metadata information found for:", paste(unlist(setdiff(sample_names, metadata_df$SampleName)), collapse = ", ")),
    call. = FALSE
  )
}

metadata_df <- validate_metadata_df(metadata_df, context = args$metadata)

  # Use ProcessingBatch as batch for Step 1 of micRoclean if cases where it exists in metadata
  #   If it is missing from meta data, set to default value with one level, will automatically skip
if ("ProcessingBatch" %in% names(metadata_df)) {
  metadata_df <- metadata_df |> 
    rename(batch = ProcessingBatch)
} else {
  metadata_df <- metadata_df |> 
    mutate(batch = "NoProcessingBatches")
}

metadata_df <- metadata_df |> 
  rename(sample_type = "SampleType") |>
  mutate(
    # Only counting negative controls as controls for purposes of decontam (cell line controls don't fix with their framework)
    is_control = sample_type == "NegativeControl",
    is_control = as.logical(is_control),
    sample_type = as.factor(sample_type),
    batch = as.factor(batch)) |>
  column_to_rownames(var = "SampleName") |>
  select(is_control, sample_type, batch)

#------------------------------------------
# Construct Technical Replicate Data Frame
#------------------------------------------
tr <- data.frame(
  Batch_1 = character(),
  Batch_2 = character()
)

#---------------------------
# Run micRoclean
#---------------------------
bl <- c("FakeBlockedTaxa")
biomarkerID_results <- micRoclean(counts = counts_df,
                                  meta = metadata_df,
                                  research_goal = "biomarker",
                                  control_name = "Control",
                                  blocklist = bl,
                                  technical_replicates = tr)

write.table(
  biomarkerID_results$decontaminated_count,
  file = paste0(args$out, "/DecontamSeqTable.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

writeLines(colnames(biomarkerID_results$decontaminated_count),
           con = paste0(args$out, "/decontaminated.ASV.names"))

write.table(
  biomarkerID_results$contaminant_id,
  file = paste0(args$out, "/DecontamFilterReport.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
