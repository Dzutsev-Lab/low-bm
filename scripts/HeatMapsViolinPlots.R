library(phyloseq)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(argparse)


parser <- ArgumentParser()

parser$add_argument("--physeqs",
                    type = "character",
                    nargs = "+",
                    help = "list of .RData files contatining phyloseq objects named 'physeq' (one for each sequencing batch)")
parser$add_argument("--patient-sample-batches",
                    type = "character",
                    nargs = "+",
                    help = "list of sequencing batches included in physeqs that contain patient sample records")
parser$add_argument("--DA-results",
                    type = "character",
                    help = ".tsv file with composite DA results from meta-analysis of the two batchs")
parser$add_argument("--select-taxa-names",
                    type = "character",
                    nargs = "+",
                    help = "optional list of files containing lists of taxa to subset differential abundance analysis to (union of all unique names provided)")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "taxonomic level to agglomerate to for DA analysis (e.g. Genus, Family, etc.)")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Normalization method carried out on abundance counts (default = noNorm)")
parser$add_argument("--limma-voom",
                    action = "store_true",
                    default = FALSE,
                    help = "Optional flag to normalize counts using TMM and Voom normalization (overides --norm-method to noNorm)")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "pseudocount value to replace zero counts (default: 1.0)")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "output directory within Exp_Output")

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
# Step 1: Phyloseq Preprocessing
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
if (args$limma_voom) {
    args$norm_method <- "noNorm"
}

NormPhyseq <- counts_normalization(physeq = FiltPhyseq, 
                                   norm_method = args$norm_method, 
                                   pseudocount = args$pseudocount)


# STEP 2: Optionally use limma/voom normalization (mutually exclusive with other normalization methods)
limma_voom_normalization <- function(prep_physeq) {
    metadata <- as.data.frame(as.matrix(sample_data(prep_physeq)))
    counts <- as.matrix(t(as(otu_table(prep_physeq), "matrix")))
    design <- model.matrix(~ SampleType + PatientID,
                           data = metadata)
    
    counts <- DGEList(counts = counts)
    counts <- calcNormFactors(counts,
                              method = "TMM")
    voom <- voom(counts, design, plot=FALSE)

    norm_counts <- as.matrix(t(voom$E))
    
    otu_table(prep_physeq) <- otu_table(norm_counts, taxa_are_rows = FALSE)

    return(prep_physeq)
}

if (args$limma_voom) {
    NormPhyseq <- limma_voom_normalization(FiltPhyseq)
}

#-------------------------------------
# Step 3: Melt Phyloseq to Data Fame
#-------------------------------------
taxa_of_interest <- unique(unlist(lapply(args$select_taxa_names, function(name_file) {
    read.csv(name_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
})))
MeltPhyseq_df <- psmelt(NormPhyseq) |>
    mutate(
        SequencingBatch = factor(SequencingBatch),
        SampleType = factor(SampleType)
    )

# STEP 4: Import differential abundance results (for statistical context included in violin plot captions)
DA_results_df <- read.delim(paste0("Exp_Output/", args$DA_results), header = TRUE)
str(DA_results_df)
# STEP 5: Violin plotting
out_dir <- paste0("Exp_Output/", args$out_dir, "/", args$norm_method)


for (taxon in taxa_of_interest) {

    single_tax_df <- MeltPhyseq_df |> filter(OTU == taxon)

    violin_plot <- ggplot(single_tax_df, 
                          aes(x = SampleType, y = Abundance, fill = SampleType)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_jitter(width =0.2, size = 1.5, alpha = 0.8) +
        facet_wrap(~SequencingBatch, ncol = 2) +
        labs(
            title = paste("Abundance Comparison for", taxon),
            x = "Sample Type",
            y = paste0("Normalized Abundance (", args$norm_method, ")")
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
        )

    B1_LFC <- round(as.numeric(DA_results_df[DA_results_df$taxon == taxon, "logFC_b1"]), 2)
    B2_LFC <- round(DA_results_df[DA_results_df$taxon == taxon, "logFC_b2"], 2)

    B1_adj_p <- round(DA_results_df[DA_results_df$taxon == taxon, "adj_p_b1"], 3)
    B2_adj_p <- round(DA_results_df[DA_results_df$taxon == taxon, "adj_p_b2"], 3)

    het_Q <- round(DA_results_df[DA_results_df$taxon == taxon, "het_Q"], 3)
    het_Q_q <- round(DA_results_df[DA_results_df$taxon == taxon, "het_Q_q"], 3)


    multibatch_violin_plot <- ggplot(single_tax_df, 
                                     aes(x = SampleType, y = Abundance)) +
        geom_violin(aes(fill = SampleType, group = SampleType),
                        trim = FALSE, alpha = 0.6) +
        geom_jitter(aes(fill = SampleType, shape = SequencingBatch),
                    width =0.2, size = 1.5, alpha = 0.8) +
        labs(
            title = paste("Abundance Comparison for", taxon),
            subtitle = paste0("LFC Exp1: ", B1_LFC, " (q=", B1_adj_p, "), ",
                              "LFC Exp2: ", B2_LFC, " (q=", B2_adj_p, "); ",
                              "HetQ: ", het_Q, " (q=", het_Q_q, ")"),
            x = "Sample Type",
            y = paste0("Normalized Abundance (", args$norm_method, ")")
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
        )
    
    if (!dir.exists(paste0(out_dir, "/", taxon))) {
        dir.create(paste0(out_dir, "/", taxon), recursive = TRUE)
    }

    ggsave(
        filename = paste0(out_dir, "/", taxon, "/",
                          taxon, "_", args$norm_method,"_violin_plot.png"),
        plot = violin_plot,
        width = 8,
        height = 4
    )

    ggsave(
        filename = paste0(out_dir, "/", taxon, "/",
                          taxon, "_", args$norm_method, "_multibatch_violin_plot.png"),
        plot = multibatch_violin_plot,
        width = 8,
        height = 4
    )

}

# STEP 6: Heat map plotting
Heat_Mapping <- function(patient_batch_name, 
                         norm_method,
                         prep_physeq, 
                         taxa_of_interest,
                         out_dir) {
    subset_physeq <- subset_samples(prep_physeq, 
                                    SequencingBatch == patient_batch_name | 
                                    SampleType == "CellLineControl")
    pruned_physeq <- prune_taxa(taxa_of_interest, subset_physeq)
    column_order <- rownames(sample_data(pruned_physeq))[order(sample_data(pruned_physeq)$SampleType,
                                                               sample_data(pruned_physeq)$PatientID)]

    heatmap_plot <- plot_heatmap(pruned_physeq,
                                 sample.order = column_order,
                                 method = "PCoA",
                                 distance = "bray",
                                 title = paste0(patient_batch_name, " Select Taxa Heatmap ",
                                                ifelse(norm_method == "noNorm",
                                                       "(Unnormalized)",
                                                       paste0("(", norm_method, ")")
                                                )
                                         ),
                                 sample.label = "SampleID",
                                 taxa.label = "Genus")
    
    # Add vertical line separating samples by sample type
    type_ordered <- sample_data(pruned_physeq)[column_order, "SampleType"]
    type_block_sizes <- table(type_ordered)
    type_block_cuts <- cumsum(type_block_sizes)
    heatmap_plot <- heatmap_plot + geom_vline(xintercept = type_block_cuts+0.5, linewidth = 1, color = "white")

    ggsave(
        filename = paste0(out_dir, "/", patient_batch_name, "_", norm_method, "_SelectTaxaHeatmap.png"),
        plot = heatmap_plot,
        width = 14,
        height = 12,
        units = "in",
        dpi = 300
    )
}

for (patient_batch_name in args$patient_sample_batches) {
    Heat_Mapping(prep_physeq = NormPhyseq,
                 patient_batch_name = patient_batch_name,
                 norm_method = args$norm_method,
                 taxa_of_interest = taxa_of_interest,
                 out_dir = out_dir)
}