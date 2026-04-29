library(phyloseq)
library(limma)
library(edgeR)
library(dplyr)
library(argparse)


parser <- ArgumentParser()

parser$add_argument("--B1-physeq",
                    type = "character",
                    help = "physeq object from batch #1 for comparison")
parser$add_argument("--B2-physeq",
                    type = "character",
                    help = "physeq object from batch #2 for comparison")
parser$add_argument("--taxa-of-interest",
                    type = "character",
                    help = "list of taxa to produce violin plots for")
parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number for the trial to label output plots")
parser$add_argument("--out",
                    type = "character",
                    help = "output directory")

args <- parser$parse_args()

load(args$B1_physeq)
B1_physeq <- physeq

load(args$B2_physeq)
B2_physeq <- physeq

taxa_of_interest <- args$taxa_of_interest



prep_phyloseq <- function(physeq) {
    grouped_physeq <- subset_samples(physeq, SampleType %in% c("NormalTissue", "Tumor"))
    grouped_physeq <- prune_taxa(taxa_sums(grouped_physeq) > 0, grouped_physeq)
   
    glommed_physeq <- tax_glom(grouped_physeq, taxrank = "Genus")
    taxa_names(glommed_physeq) <- as.character(tax_table(glommed_physeq)[, "Genus"])

    return(glommed_physeq)
}

prep_B1_physeq <- prep_phyloseq(B1_physeq)
prep_B2_physeq <- prep_phyloseq(B2_physeq)

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

voom_B1_physeq <- limma_voom_normalization(prep_B1_physeq)
voom_B2_physeq <- limma_voom_normalization(prep_B2_physeq)

melt_and_batch_label <- function(voom_physeq, batch_name) {
    melted_df <- psmelt(voom_physeq)
    melted_df$Batch <- batch_name
    return(melted_df)
}

melt_voom_B1_df <- melt_and_batch_label(voom_B1_physeq, "ItalyExp1")
melt_voom_B2_df <- melt_and_batch_label(voom_B2_physeq, "ItalyExp2")

melt_raw_B1_df <- melt_and_batch_label(prep_B1_physeq, "ItalyExp1")
melt_raw_B2_df <- melt_and_batch_label(prep_B2_physeq, "ItalyExp2")

comp_and_filter_df <- function(df1, df2, taxa_of_interest) {
    comp_df <- bind_rows(df1, df2)

    filt_df <- comp_df |>
        filter(
            OTU %in% taxa_of_interest
        ) |>
        mutate(
            SampleType = factor(SampleType)
        )

    return(filt_df)
} 

comp_raw_df <- comp_and_filter_df(
    df1 = melt_raw_B1_df, 
    df2 = melt_raw_B2_df,
    taxa_of_interest = taxa_of_interest
)

comp_voom_df <- comp_and_filter_df(
    df1 = melt_voom_B1_df, 
    df2 = melt_voom_B2_df,
    taxa_of_interest = taxa_of_interest
)

limma_voom_results_df <- read.delim("/data/taylorng/low-bm/Exp_Output/042226.3_ItalyLungExpComp_RawTMMLimmaVoom/shared_taxa_meta.tsv",
                                    header = TRUE)
filter_limma_voom_results_df <- limma_voom_results_df |>
    filter(taxon in taxa_of_interest)

raw_out_dir <- "/data/taylorng/low-bm/Exp_Output/042726_ItalyLungExpComp_ViolinPlots/Raw"
voom_out_dir <- "/data/taylorng/low-bm/Exp_Output/042726_ItalyLungExpComp_ViolinPlots/Voom"
if (!dir.exists(raw_out_dir)) {
    dir.create(raw_out_dir, recursive = TRUE)
}
if (!dir.exists(voom_out_dir)) {
    dir.create(voom_out_dir, recursive = TRUE)
}

for (taxon in taxa_of_interest) {

    single_tax_raw_df <- comp_raw_df |> filter(OTU == taxon)
    single_tax_voom_df <- comp_voom_df |> filter(OTU == taxon)

    raw_violin_plot <- ggplot(single_tax_raw_df, 
                              aes(x = SampleType, y = Abundance, fill = SampleType)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_jitter(width =0.2, size = 1.5, alpha = 0.8) +
        facet_wrap(~Batch, ncol = 2) +
        labs(
            title = paste("Raw Counts for", taxon),
            x = "Sample Type",
            y = "Count"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
        )


    voom_violin_plot <- ggplot(single_tax_voom_df, 
                               aes(x = SampleType, y = Abundance, fill = SampleType)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_jitter(width =0.2, size = 1.5, alpha = 0.8) +
        facet_wrap(~Batch, ncol = 2) +
        labs(
            title = paste("TMM + limma/voom Normalized Counts for", taxon),
            x = "Sample Type",
            y = "log2 CPM (TMM + voom normalized)"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
        )
    
    B1_LFC <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "logFC_b1"], 2)
    B2_LFC <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "logFC_b2"], 2)
    
    B1_adj_p <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "adj_p_b1"], 3)
    B2_adj_p <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "adj_p_b2"], 3)

    het_Q <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "het_Q"], 3)
    het_Q_q <- round(limma_voom_results_df[limma_voom_results_df$taxon == taxon, "het_Q_q"], 3)


    multibatch_voom_violin_plot <- ggplot(single_tax_voom_df, 
                                          aes(x = SampleType, y = Abundance)) +
        geom_violin(aes(fill = SampleType, group = SampleType),
                        trim = FALSE, alpha = 0.6) +
        geom_jitter(aes(fill = SampleType, shape = Batch),
                    width =0.2, size = 1.5, alpha = 0.8) +
        labs(
            title = paste("TMM + limma/voom Normalized Counts for", taxon),
            subtitle = paste0("LFC Exp1: ", B1_LFC, " (q=", B1_adj_p, "), ",
                              "LFC Exp2: ", B2_LFC, " (q=", B2_adj_p, "); ",
                              "HetQ: ", het_Q, " (q=", het_Q_q, ")"),
            x = "Sample Type",
            y = "log2 CPM (TMM + voom normalized)"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
        )
    
    ggsave(
        filename = paste0(raw_out_dir, "/", taxon, "_raw_violin_plot.png"),
        plot = raw_violin_plot,
        width = 8,
        height = 4
    )
    
    ggsave(
        filename = paste0(voom_out_dir, "/", taxon, "_voom_violin_plot.png"),
        plot = voom_violin_plot,
        width = 8,
        height = 4
    )

    ggsave(
        filename = paste0(voom_out_dir, "/", taxon, "_multibatch_voom_violin_plot.png"),
        plot = multibatch_voom_violin_plot,
        width = 8,
        height = 4
    )

}

write.table(
    as.data.frame(as(otu_table(prep_B1_physeq), "matrix")),
    file = paste0(raw_out_dir, "/raw_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
)

write.table(
    as.data.frame(as(otu_table(voom_B1_physeq), "matrix")),
    file = paste0(voom_out_dir, "/voom_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
)