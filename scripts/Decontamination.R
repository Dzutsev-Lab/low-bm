library(micRoclean)
library(readxl)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--seq-table",
                    type = "character",
                    help = "sequence table from Dada2 (Samples x ASV)")
parser$add_argument("--bacterial-names",
                    type = "character",
                    help = "ASV IDs associated with bacterial reads according to mapping filters")
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

#---------------------------
# Construct count matrix
#---------------------------
raw_seq_table <- read.delim(args$seq_table,
                            header = TRUE,
                            row.names = 1)

bacterial_ASVids <- readLines(args$bacterial_names)
bacterial_ASVids <- bacterial_ASVids[nzchar(bacterial_ASVids)]

counts_df <- raw_seq_table[, bacterial_ASVids]
# Remove S## label from sample names in the sequence table to match with metadata sample names
counts_df <- counts_df |>
  rownames_to_column(var = "Sample_ID") |>
  mutate(Sample_ID = sub("_S\\d+$", "", Sample_ID)) |>
  column_to_rownames(var = "Sample_ID")


#---------------------------
# Construct sample metadata
#---------------------------
metadata_df <- read_excel(args$metadata)
metadata_df <- metadata_df |> 
  rename(batch = "ProcessingBatch",
         sample_type = "SampleType") |>
  mutate(
    is_control = str_detect(sample_type, "Control"),
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