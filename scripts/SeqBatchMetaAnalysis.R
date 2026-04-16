library(tibble)
library(dplyr)
library(ggplot2)
library(argparse)
library(phyloseq)
library(stringr)

parser <- ArgumentParser()

parser$add_argument("--batch1-Name",
                    type = "character",
                    help = "Name for batch 1 used to locate output files")
parser$add_argument("--batch2-Name",
                    type = "character",
                    help = "Name for batch 2 used to locate output files")
parser$add_argument("--norm-methods",
                    type = "character",
                    nargs = '+',
                    help = "normalized table to compare between batches (e.g. noNorm, RelAbund, RawTSS, HostMapped, log2)")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "pseudocount value to replace zero counts (default: 1.0)")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = NULL,
                    help = "taxonomic aggregation level for output sequence tables (e.g. Genus, Species, etc.) (default: NULL, no aggregation, working at ASV level)")
parser$add_argument("--out",
                    type = "character",
                    help = "desired output directory")

args <- parser$parse_args()



Compute_LFC <- function(group1_df, group2_df, pseudocount) {
    log_group1 <- log2(group1_df + pseudocount)
    log_group2 <- log2(group2_df + pseudocount)

    mean_group1 <- colMeans(log_group1, na.rm = TRUE)
    mean_group2 <- colMeans(log_group2, na.rm = TRUE)

    LFC <- mean_group2 - mean_group1
    return(LFC)
}

Compute_Paired_LFC <- function(group1_df, group2_df, pseudocount) {
    #TODO: implement this function to compute paired LFCs for samples with matched tumor-normal samples across batches
}

LFCvLFC_plot <- function(seqTableB1, seqTableB2, physeqB1, physeqB2, expB1, expB2, method, common_taxa, out_dir) {
    # Subsetting seqtables to common taxa and samples of interest
    TumorSeqTableB1 <- seqTableB1[TumorSamplesB1, common_taxa]
    TumorSeqTableB2 <- seqTableB2[TumorSamplesB2, common_taxa]
    NormalSeqTableB1 <- seqTableB1[NormalTissueSamplesB1, common_taxa]
    NormalSeqTableB2 <- seqTableB2[NormalTissueSamplesB2, common_taxa]

    # Calculating log fold changes for each batch
    BatchCompResults <- data.frame(
        LFC_B1 = Compute_LFC(NormalSeqTableB1, TumorSeqTableB1, args$pseudocount),
        LFC_B2 = Compute_LFC(NormalSeqTableB2, TumorSeqTableB2, args$pseudocount)
    )

    # LFC vs LFC Plotting
    ggplot(BatchCompResults, aes(x = LFC_B1, y = LFC_B2)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = paste0("Log Fold Change Comparison Between ", expB1, " and ", expB2, "(", method, " Normalization)"),
        x = paste0("LFC (Batch 1: ", expB1, ")"),
        y = paste0("LFC (Batch 2: ", expB2, ")")
      )
    ggsave(paste0(out_dir, "/", expB1, "v", expB2, "_", method, "_LFCvLFC.png"), width = 14, height = 12, units = "in")

  # # Correlation Calculations
  # cor_pearson <- cor(BatchCompResults$LFC_B1,
  #                   BatchCompResults$LFC_B2,
  #                   method = "pearson")
  # cor_spearman <- cor(BatchCompResults$LFC_B1,
  #                   BatchCompResults$LFC_B2,
  #                   method = "spearman")

  # # Concordance of Direction Calculations
  # prop_same_direction <- sum(
  #   sign(BatchCompResults$LFC_B1) == sign(BatchCompResults$LFC_B2)
  # ) / nrow(BatchCompResults)
    
}
batch1_name_components <- str_split(args$batch1_Name, "_")[[1]]
batch1_ID <- batch1_name_components[1]
batch1_ExpID <- batch1_name_components[2]
batch1_TrialDescription <- batch1_name_components[3]
batch2_name_components <- str_split(args$batch2_Name, "_")[[1]]
batch2_ID <- batch2_name_components[1]
batch2_ExpID <- batch2_name_components[2]
batch2_TrialDescription <- batch2_name_components[3]

load(paste0("Exp_Output/", args$batch1_Name, "/", batch1_ID, "_raw_kraken_phyloseq.RData"))
physeqB1 <- raw_kraken_phyloseq
load(paste0("Exp_Output/", args$batch2_Name, "/", batch2_ID, "_raw_kraken_phyloseq.RData"))
physeqB2 <- raw_kraken_phyloseq

# Glomming to desired taxonomic level if specified
# TODO: shift this to phyloseq analysis given redo in both this script and seqtable normalization
if (!is.null(args$tax_agg_level)) {
    physeqB1 <- tax_glom(physeqB1, taxrank = args$tax_agg_level)
    taxa_names(physeqB1) <- as.character(tax_table(physeqB1)[, args$tax_agg_level])
    physeqB2 <- tax_glom(physeqB2, taxrank = args$tax_agg_level)
    taxa_names(physeqB2) <- as.character(tax_table(physeqB2)[, args$tax_agg_level])
}

# Getting sets of samples grouped by tumor vs normal samples for each batch
NormalTissueSamplesB1 <- sample_names(subset_samples(physeqB1, SampleType == "NormalTissue"))
NormalTissueSamplesB2 <- sample_names(subset_samples(physeqB2, SampleType == "NormalTissue"))
TumorSamplesB1 <- sample_names(subset_samples(physeqB1, SampleType == "Tumor"))
TumorSamplesB2 <- sample_names(subset_samples(physeqB2, SampleType == "Tumor"))

# Getting set of common taxa between batches to compare
common_taxa <- intersect(taxa_names(physeqB1), taxa_names(physeqB2))

for (method in args$norm_methods) {

  seqTableB1_path <- paste0("Exp_Output/", args$batch1_Name, "/CountNormalization/", batch1_ID, "_", method, "_SeqTable.tsv")
  seqTableB2_path <- paste0("Exp_Output/", args$batch2_Name, "/CountNormalization/", batch2_ID, "_", method, "_SeqTable.tsv")



  if (!file.exists(seqTableB1_path) || !file.exists(seqTableB2_path)) {
    warning(paste("Normalized sequence table for method", method, "not found in one or both batches. Skipping this method."))
    next
  }

  print(paste("Processing normalization method:", method))
  seqTableB1 <- read.delim(seqTableB1_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  seqTableB2 <- read.delim(seqTableB2_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)


  LFCvLFC_plot(seqTableB1 = seqTableB1, seqTableB2 = seqTableB2, 
               physeqB1 = physeqB1, physeqB2 = physeqB2,
                expB1 = batch1_ExpID, expB2 = batch2_ExpID, 
               method = method, common_taxa = common_taxa, out_dir = args$out)
}
