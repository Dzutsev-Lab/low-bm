library(phyloseq)
library(argparse)
library(tidyr)
library(tibble)
library(dplyr)

parser <- ArgumentParser()

parser$add_argument("--phyloseq-data",
                    type = "character",
                    help = "phyloseq object .RData file path")
parser$add_argument("--norm-methods",
                    type = "character",
                    nargs = '+',
                    help = "normalization methods to apply to sequence count table")
parser$add_argument("--read-counts",
                    type = "character",
                    help = "tsv file with read counts to use as normalization denominators if needed (e.g. raw read counts, host read counts, etc.)")
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
parser$add_argument("--trialID",
                    type = "character",
                    help = "unique identifier for this trial (used in output file names)")

args <- parser$parse_args()




normalize <- function(physeq, method, pseudocount = 1.0) {
    if (method == "noNorm") {
        return(physeq)
    }

    if (method == "log2") {
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    }

    if (method == "RelAbund") {
        return(transform_sample_counts(physeq, function(x) x / sum(x)))
    }

    otu_divide_by_sample_factor <- function(physeq, factor_column) {
        sample_factors <- sample_data(physeq)[[factor_column]]
   
        if (any(sample_factors == 0, na.rm = TRUE)) {
            warning("Samples with zero denominator detected, ; setting to NA")
            sample_factors[sample_factors == 0] <- NA
        }

        otu_mat <- as(otu_table(physeq), "matrix")

        if (taxa_are_rows(physeq)) {
            otu_mat <- sweep(otu_mat, 2, sample_factors, FUN = "/")
            otu_mat <- otu_mat * 1e6 # scaling factor to bring values back to a more interpretable range
        } else {
            otu_mat <- sweep(otu_mat, 1, sample_factors, FUN = "/")
            otu_mat <- otu_mat * 1e6
        }

        otu_table(physeq) <- otu_table(otu_mat, taxa_are_rows = taxa_are_rows(physeq))
        return(physeq)
    }

    if (method == "RawTSS") {
        return(otu_divide_by_sample_factor(physeq, "Raw_reads"))
    }
    if (method == "HostMapped") {
        return(otu_divide_by_sample_factor(physeq, "Host_mapped_reads"))
    }
}

load(args$phyloseq_data)
read_count_df <- read.delim(args$read_counts, sep = "\t", header = TRUE)
read_count_df <- read_count_df |>
    mutate(
        Sample_Name = sub("_S\\d+$", "", SampleID),
        HostMappedReads = chimera.filtered - HostUnmapped_reads,
        Raw_reads = as.numeric(Raw_reads),
        HostMappedReads = as.numeric(HostMappedReads) 
    ) |>
    select(Sample_Name, Raw_reads = Raw_reads, Host_mapped_reads = HostMappedReads)

physeq <- raw_kraken_phyloseq
if (!is.null(args$tax_agg_level)) {
    physeq <- tax_glom(physeq, taxrank = args$tax_agg_level)
    taxa_names(physeq) <- as.character(tax_table(physeq)[, args$tax_agg_level])
}

print(intersect(rownames(sample_data(physeq)), read_count_df$Sample_Name))
# Adding needed normalization denominators to metadata of phyloseq object
sample_data(physeq) <- as(sample_data(physeq), "data.frame") |>
    rownames_to_column(var = "Sample_Name") |>
    left_join(read_count_df, by = "Sample_Name") |>
    column_to_rownames(var = "Sample_Name")
# Applying normalizations
for (method in args$norm_methods) {
    print(paste("Applying normalization method:", method))
    normalized_physeq <- normalize(physeq, method, pseudocount = args$pseudocount)
    write.table(as(otu_table(normalized_physeq), "matrix"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = TRUE,
                col.names = NA,
                file = file.path(args$out, paste0(args$trialID, "_", method, "_SeqTable.tsv")))
}
