library(phyloseq)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(argparse)


parser <- ArgumentParser()

parser$add_argument("--B1-physeq",
                    type = "character",
                    help = ".RData file with physeq object from batch #1 for comparison")
parser$add_argument("--B2-physeq",
                    type = "character",
                    help = ".RData file with physeq object from batch #2 for comparison")
parser$add_argument("--DA-results",
                    type = "character",
                    help = ".tsv file with composite DA results from meta-analysis of the two batchs")
parser$add_argument("--taxa-of-interest",
                    type = "character",
                    nargs = "+",
                    help = "list of taxa to produce violin plots for")
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

if (args$limma_voom) {
    args$norm_method <- "noNorm"
}

load(paste0("Exp_Output/", args$B1_physeq))
B1_physeq <- physeq

load(paste0("Exp_Output/", args$B2_physeq))
B2_physeq <- physeq

taxa_of_interest <- args$taxa_of_interest


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
    } else {
        message("Unknown normalization method provided for limma voom pre-normlaization, no normalization used.")
    }

}


prep_phyloseq <- function(physeq, norm_method, pseudocount) {
    grouped_physeq <- subset_samples(physeq, SampleType %in% c("NormalTissue", "Tumor"))
    grouped_physeq <- prune_taxa(taxa_sums(grouped_physeq) > 0, grouped_physeq)
   
    glommed_physeq <- tax_glom(grouped_physeq, taxrank = "Genus")
    taxa_names(glommed_physeq) <- as.character(tax_table(glommed_physeq)[, "Genus"])

    normalized_physeq <- counts_normalization(physeq = glommed_physeq,
                                              norm_method = norm_method,
                                              pseudocount = pseudocount)

    return(normalized_physeq)
}

prep_B1_physeq <- prep_phyloseq(
    physeq = B1_physeq,
    norm_method = args$norm_method,
    pseudocount = args$pseudocount
)
prep_B2_physeq <- prep_phyloseq(
    physeq = B2_physeq,
    norm_method = args$norm_method,
    pseudocount = args$pseudocount
)


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
    prep_B1_physeq <- limma_voom_normalization(prep_B1_physeq)
    prep_B2_physeq <- limma_voom_normalization(prep_B2_physeq)
}



melt_and_batch_label <- function(prep_physeq, batch_name) {
    melted_df <- psmelt(prep_physeq)
    melted_df$Batch <- batch_name
    return(melted_df)
}


melt_B1_df <- melt_and_batch_label(prep_B1_physeq, "ItalyExp1")
melt_B2_df <- melt_and_batch_label(prep_B2_physeq, "ItalyExp2")

comp_and_filter_df <- function(df1, df2, taxa_of_interest) {
    comp_df <- bind_rows(df1, df2)

    filt_df <- comp_df |>
        filter(
            OTU %in% taxa_of_interest
        ) |>
        mutate(
            SampleType = factor(SampleType),
            Batch = factor(Batch)
        )

    return(filt_df)
} 

comp_df <- comp_and_filter_df(
    df1 = melt_B1_df, 
    df2 = melt_B2_df,
    taxa_of_interest = args$taxa_of_interest
)


DA_results_df <- read.delim(paste0("Exp_Output/", args$DA_results), header = TRUE) |> filter(taxon %in% args$taxa_of_interest)

out_dir <- paste0("Exp_Output/", args$out_dir)
if (!dir.exists(out_dir)) {
    dir.create(raw_out_dir, recursive = TRUE)
}

for (taxon in taxa_of_interest) {

    single_tax_df <- comp_df |> filter(OTU == taxon)

    violin_plot <- ggplot(single_tax_df, 
                          aes(x = SampleType, y = Abundance, fill = SampleType)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_jitter(width =0.2, size = 1.5, alpha = 0.8) +
        facet_wrap(~Batch, ncol = 2) +
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
    
    B1_LFC <- round(DA_results_df[DA_results_df$taxon == taxon, "logFC_b1"], 2)
    B2_LFC <- round(DA_results_df[DA_results_df$taxon == taxon, "logFC_b2"], 2)
    
    B1_adj_p <- round(DA_results_df[DA_results_df$taxon == taxon, "adj_p_b1"], 3)
    B2_adj_p <- round(DA_results_df[DA_results_df$taxon == taxon, "adj_p_b2"], 3)

    het_Q <- round(DA_results_df[DA_results_df$taxon == taxon, "het_Q"], 3)
    het_Q_q <- round(DA_results_df[DA_results_df$taxon == taxon, "het_Q_q"], 3)


    multibatch_violin_plot <- ggplot(single_tax_df, 
                                     aes(x = SampleType, y = Abundance)) +
        geom_violin(aes(fill = SampleType, group = SampleType),
                        trim = FALSE, alpha = 0.6) +
        geom_jitter(aes(fill = SampleType, shape = Batch),
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
    
    ggsave(
        filename = paste0(out_dir, "/", taxon, "_violin_plot.png"),
        plot = violin_plot,
        width = 8,
        height = 4
    )

    ggsave(
        filename = paste0(out_dir, "/", taxon, "_multibatch_violin_plot.png"),
        plot = multibatch_violin_plot,
        width = 8,
        height = 4
    )

}