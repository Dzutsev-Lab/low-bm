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


#---------------------------
# Construct sample metadata
#---------------------------
# sample_names <- rownames(counts_df)
# sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)

# # Sample Type
# sample_type <- sub("\\d+[A-Za-z]*$", "", sample_info)

# # Technical Rep (Will also denote patient ID for patient samples)
# tech_rep <- sub(".*?(\\d+[A-Za-z]*)$", "\\1", sample_info)

metadata_sheet_df <- read_excel(args$metadata)
metadata_df <- metadata_sheet_df |> 
  rename(batch = "ProcessingBatch",
         sample_type = "SampleType",
         Sample_ID = "SampleName") |>
  mutate(
    is_control = ifelse(str_detect(sample_type, "Control"),
                        TRUE,
                        FALSE),
    batch = factor(batch),
    sample_type = factor(sample_type)) |>
  column_to_rownames(var = "Sample_ID")

#------------------------------------------
# Construct Technical Replicate Data Frame
#------------------------------------------
technical_replicate_df <- metadata_sheet_df[!str_detect(metadata_sheet_df$SampleType, "Control"), 
                                            c("SampleID", "SampleName", "SequencingBatch"), drop = FALSE] |>
  pivot_wider(names_from = "SequencingBatch",
              values_from = "SampleName") |> 
  na.omit() |>
  select(-SampleID)


#---------------------------
# Run micRoclean
#---------------------------
bl <- c("FakeBlockedTaxa")

biomarkerID_results <- micRoclean(counts = counts_df,
                                  meta = metadata_df,
                                  research_goal = "biomarker",
                                  control_name = "Control",
                                  blocklist = bl,
                                  technical_replicates = technical_replicate_df)

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