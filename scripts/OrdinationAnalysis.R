library(tibble)
library(dplyr)
library(ggplot2)
library(argparse)
library(phyloseq)

parser <- ArgumentParser()

parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number to attach to output files")
parser$add_argument("--physeqs",
                    type = "character",
                    nargs = "+",
                    help = "list of .RData files contatining phyloseq objects named 'physeq' (one for each sequencing batch)")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Determine which normalized sequence table to feed into limma voom (default = noNorm)")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "pseudocount value to replace zero counts (default: 1.0)")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "taxonomic level to agglomerate to for DA analysis (e.g. Genus, Family, etc.)")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()


counts_normalization <- function(physeq, 
                                 norm_method, 
                                 pseudocount) {


    otu_divide_by_sample_factor <- function(physeq, factor_column) {
      sample_factors <- sample_data(physeq)[[factor_column]]

      otu_mat <- as(otu_table(physeq), "matrix")

      if (taxa_are_rows(physeq)) {
          otu_mat <- sweep(otu_mat, 2, sample_factors, FUN = "/")
      } else {
          otu_mat <- sweep(otu_mat, 1, sample_factors, FUN = "/")
      }
      otu_mat <- otu_mat * 1e6  # scaling factor to bring values back to a more interpretable range
      otu_table(physeq) <- otu_table(otu_mat, taxa_are_rows = taxa_are_rows(physeq))
      return(physeq)
    }

    if (norm_method == "noNorm") {
        return(physeq)
    } else if (norm_method == "log2") {
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else if (norm_method == "RelAbund") {
        return(transform_sample_counts(physeq, function(x) x / sum(x)))
    } else if (norm_method == "RawTSS") {
        return(otu_divide_by_sample_factor(physeq, "Raw_reads"))
    } else if (norm_method == "HostMapped") {
        return(otu_divide_by_sample_factor(physeq, "Host_mapped_reads"))
    } else if (norm_method == "log2HostMapped") {
        physeq <- otu_divide_by_sample_factor(physeq, "Host_mapped_reads")
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else {
        message("Unknown normalization method provided for limma voom pre-normlaization, no normalization used.")
    }

}


#------------------------
# Phyloseq Preprocessing
#------------------------
load_physeq <- function(path) {
  e <- new.env()
  load(path, envir = e)

  if (!exists("physeq", envir = e)) {
    stop("No object named 'physeq' found in: ", path)
  }

  get("physeq", envir = e)
}

physeq_list <- lapply(args$physeqs, load_physeq)
CompPhyseq <- Reduce(phyloseq::merge_phyloseq, physeq_list)

# Glomming to desired taxa level
if (!is.null(args$tax_agg_level)) {
    GlomPhyseq <- tax_glom(CompPhyseq,
                           taxrank = args$tax_agg_level)
    taxa_names(GlomPhyseq) <- as.character(tax_table(GlomPhyseq)[, args$tax_agg_level])
} else {
    GlomPhyseq <- CompPhyseq
}

# Removing Negative Controls
FiltPhyseq <- subset_samples(GlomPhyseq, SampleType != "NegativeControl")

# Normalizing OTU Counts
NormPhyseq <- counts_normalization(physeq = FiltPhyseq, 
                                   norm_method = args$norm_method, 
                                   pseudocount = args$pseudocount)



#--------------------------
# Ordination Plotting
#--------------------------
NormOrd <- ordinate(NormPhyseq, "PCoA", "bray")
SampleTypeOrdPlot <- plot_ordination(NormPhyseq, NormOrd, 
                                     type="samples",
                                     color="SampleType",
                                     shape="SequencingBatch") +
                     geom_point(size = 5) +
                     ggtitle("Composite Batch Sample Ordination")
if (!dir.exists(args$out)) {
    dir.create(args$out,
               recursive = TRUE)
}

ggsave(
    filename = paste0(args$out, "/", args$trialID, "_SampleOrdination.png")
)