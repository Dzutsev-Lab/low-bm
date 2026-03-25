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
counts_df <- counts_df |>
  rownames_to_column(var = "Sample_ID") |>
  mutate(Sample_ID = sub("_[^_]*$", "", Sample_ID)) |>
  column_to_rownames(var = "Sample_ID")


#---------------------------
# Construct sample metadata
#---------------------------
sample_names <- rownames(counts_df)
sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)

# Sample Type
sample_type <- sub("\\d+[A-Za-z]*$", "", sample_info)

# Technical Rep (Will also denote patient ID for patient samples)
tech_rep <- sub(".*?(\\d+[A-Za-z]*)$", "\\1", sample_info)

metadata_sheet_df <- read_excel(args$metadata)
metadata_sheet_df <- metadata_sheet_df |> 
  rename(batch = "Batch_type") |>
  mutate(
    batch = gsub("\\D", "", batch),
    batch = factor(batch)
  )

metadata_df <- data.frame(
  sample_type = factor(sample_type),
  replicate = as.character(tech_rep),
  Sample_ID = rownames(counts_df),
  stringsAsFactors = FALSE
)

metadata_df <- metadata_df |>
  mutate(
    sample_type = case_when(
                            (is.na(sample_type) | sample_type == "") &
                              (grepl("NT$", replicate) | grepl("N$", replicate))
                            ~ "NormalTissue", #check before tumor otherwise all would be labeled tumor (ending with T)
                            (is.na(sample_type) | sample_type == "") &
                              grepl("T$", replicate) ~ "Tumor",
                            TRUE ~ sample_type),
    sample_type = case_when(
                            sample_type %in% c("expcontrol", "NEGATIVECONTROL")  ~ "Control",
                            TRUE ~ as.character(sample_type)),
    sample_type = factor(sample_type),
    is_control = case_when(
                           sample_type == "Control" ~ TRUE,
                           TRUE ~ FALSE)
  )
metadata_sheet_df <- metadata_sheet_df[, c("Sample_ID", "batch")]
metadata_df <- metadata_df |>
  full_join(metadata_sheet_df, by = "Sample_ID") |>
  column_to_rownames(var = "Sample_ID")

metadata_df <- select(metadata_df, -replicate)

#------------------------------------------
# Construct Technical Replicate Data Frame
#------------------------------------------
# tr <- metadata_df[, "batch", drop = FALSE] %>%
#   rownames_to_column(var = "Sample_ID") %>%
#   mutate(
#     batch_1 = case_when(batch == 1 ~ Sample_ID,
#                         TRUE ~ NA),
#     batch_2 = case_when(batch ==2 ~ Sample_ID,
#                         TRUE ~ NA))
# tr <- select(tr, batch_1, batch_2)
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