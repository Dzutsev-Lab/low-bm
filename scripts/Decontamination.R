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
metadata_df <- read_excel(args$metadata)
metadata_df <- metadata_df |> 
  rename(batch = "ProcessingBatch",
         sample_type = "SampleType") |>
  mutate(
    is_control = ifelse(str_detect(sample_type, "Control"),
                        TRUE,
                        FALSE),
    sample_type = as.factor(sample_type),
    is_control = as.logical(is_control),
    batch = as.factor(batch),
    SequencingBatch = as.factor(SequencingBatch),
    SampleID = as.factor(SampleID),
    SampleName = as.character(SampleName)) |>
    as.data.frame()


#------------------------------------------
# Construct Technical Replicate Data Frame
#------------------------------------------
technical_replicate_df <- metadata_df[!metadata_df$is_control, c("SampleID", "SampleName", "SequencingBatch"), drop = FALSE] |>
  pivot_wider(id_cols = SampleID,
              names_from = SequencingBatch,
              values_from = SampleName) |> 
  mutate(n_non_na = rowSums(!is.na(across(-SampleID)))) |>
  filter(n_non_na == 2) |>
  select(-SampleID, -n_non_na) |>
  as.data.frame()

# Set sample file names to row names to fit micRoclean input format
metadata_df <- metadata_df |> 
  column_to_rownames(var = "SampleName") |>
  select(is_control, sample_type, batch) |>
  as.data.frame()

summary(metadata_df)
summary(technical_replicate_df)
print(technical_replicate_df)

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