library(argparse)
library(readxl)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(purrr)

library(phyloseq)
library(Biostrings)
library(taxonomizr)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source(file.path("scripts", "Rhelpers", "KrakenTaxonomy.R"))



# Keep ggplot from producing Rplots.pdf
if (!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--kraken-file",
                    type = "character",
                    help = "kraken taxonomy classification file path")
parser$add_argument("--raw-seq-table",
                    type = "character",
                    help = "un-normalized seq table tsv file path")
parser$add_argument("--bacterial-names",
                    type = "character",
                    help = "names of bacterial ASVs")
parser$add_argument("--sample-names", 
                    type="character", 
                    help="Text file containing list of sample names, one per line")
parser$add_argument("--metadata",
                    type = "character",
                    help = "standardized metadata sheet as .xlsx file")
parser$add_argument("--library-counts",
                    type = "character",
                    help = "tsv file with read counts at various pipeline stages to use as normalization denominators if needed (e.g. raw read counts, host read counts, etc.)")
parser$add_argument("--dump-dir",
                    type = "character",
                    help = "directory containing nodes.dmp and names.dmp 
                            to be used in database construction")
parser$add_argument("--add-unclassified-prefix",
                    help = "Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action = "store_true",
                    default = FALSE)
parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number to identify phyloseq output")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)



#----------------------------------------
# Positive Bacterial ASVID Extraction
#----------------------------------------
bacterial_IDs <- readLines(args$bacterial_names) # split IDs by line # nolint
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines # nolint


#-----------------------------------------
# Sequence Table Construction
#-----------------------------------------
# Need to remove S## label from sample names in the sequence table to match with metadata sample names
# UN-NORMALIZED
raw_seq_table <- read.delim(args$raw_seq_table,
                            header = TRUE,
                            row.names = 1)


#--------------------------------
# Sample Data Table Construction
#--------------------------------
sample_names <- readLines(args$sample_names)
sample_names <- trimws(sample_names)
sample_names <- sample_names[sample_names != ""]

sample_meta_data_df <- read_excel(args$metadata)
sample_meta_data_df <- sample_meta_data_df |>
  filter(SampleName %in% sample_names)

missing_metadata <- setdiff(sample_names, sample_meta_data_df$SampleName)
missing_metadata_report <- data.frame(
  SampleName = missing_metadata,
  Reason = rep("missing_metadata", length(missing_metadata)),
  stringsAsFactors = FALSE
)
write.table(
  missing_metadata_report,
  file = file.path(args$out, "DroppedSamplesMissingMetadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(missing_metadata) > 0) {
  warning(
    "Dropping ",
    length(missing_metadata),
    " sample(s) before phyloseq construction because no metadata was found: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(sample_meta_data_df) == 0) {
  stop(
    "No samples with metadata remain after dropping samples missing metadata.",
    call. = FALSE
  )
}

sample_meta_data_df <- validate_metadata_df(sample_meta_data_df, context = args$metadata)

library_counts_df <- read.delim(args$library_counts, sep = "\t", header = TRUE)
library_counts_df <- library_counts_df |>
  mutate(
      SampleName = SampleID,
      Raw_reads = as.numeric(Raw_reads),
      HostMappedReads = chimera.filtered - HostUnmapped_reads,
      HostMappedReads = as.numeric(HostMappedReads) 
  ) |>
  select(SampleName, Raw_reads = Raw_reads, Host_mapped_reads = HostMappedReads)  

sample_meta_data_df <- sample_meta_data_df |>
  left_join(library_counts_df, by = "SampleName") |>
  column_to_rownames(var = "SampleName")
raw_seq_table <- raw_seq_table[rownames(sample_meta_data_df), , drop = FALSE]

missing_bacterial_IDs <- setdiff(bacterial_IDs, colnames(raw_seq_table))
if (length(missing_bacterial_IDs) > 0) {
  warning(
    "Dropping ",
    length(missing_bacterial_IDs),
    " bacterial ASV ID(s) because they are absent from the sequence table.",
    call. = FALSE
  )
}

bacterial_IDs <- bacterial_IDs[bacterial_IDs %in% colnames(raw_seq_table)]
if (length(bacterial_IDs) == 0) {
  stop("No bacterial ASVs from --bacterial-names are present in the sequence table.", call. = FALSE)
}

raw_seq_table <- raw_seq_table[, bacterial_IDs, drop = FALSE]

#------------------------------------
# Kraken Taxonomy Table Construction
#------------------------------------
# Taxonomy Database File Paths
names_dmp <- file.path(args$dump_dir, "names.dmp")
nodes_dmp <- file.path(args$dump_dir, "nodes.dmp")
sql_db <- file.path(args$dump_dir, "tax_db.sqlite")

# Taxonomy Database Construction
if (!file.exists(sql_db)) {
  read.names.sql(names_dmp, sql_db, overwrite = TRUE)
  read.nodes.sql(nodes_dmp, sql_db, overwrite = TRUE)
}


# Import ASVid -> TaxID info from .kraken2 file. Kraken annotates the bacterial
# ASV universe; unclassified ASVs are retained with missing taxonomy.
kraken_info <- read.delim(args$kraken_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(kraken_info)[1:3] <- c("status", "ASVid", "taxid")

kraken_tax_matrix <- build_kraken_tax_matrix(
  kraken_info = kraken_info,
  asv_ids = colnames(raw_seq_table),
  sql_db = sql_db,
  add_unclassified_prefix = isTRUE(args$add_unclassified_prefix)
)


#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
anyDuplicated(colnames(raw_seq_table))
anyDuplicated(rownames(kraken_tax_matrix))
physeq <- phyloseq(otu_table(raw_seq_table, taxa_are_rows = FALSE),
                     sample_data(sample_meta_data_df),
                     tax_table(kraken_tax_matrix))

save(
    physeq,
    file = paste0(args$out, "/", args$trialID, "_physeq.RData")
)
