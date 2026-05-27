library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)
library(phyloseq)
library(limma)
library(sva)
library(vegan)
library(FSA)
library(ggpubr)

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
parser$add_argument("--batch-adj",
                    type = "character",
                    help = "Optional flag with the field in file metadata to use for batch adjustment (default = no batch adjustment)")
parser$add_argument("--techrep-avg",
                    action = "store_true",
                    default = FALSE,
                    help = "Optional flag to average across technical replicates if available across phyloseq objects once compiled")
parser$add_argument("--dist-metric",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "list of desired distance metrics to generate ordination plots for (default is all distance metrics available via distanceMethodList$vegdist)")
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

batch_adjustment <- function(physeq, batch_column) {
    otu_mat <- as(otu_table(physeq), "matrix")
    meta_df <- as(sample_data(physeq), "data.frame")

    meta_df <- meta_df |>
        mutate(
            SampleID = ifelse(ControlStatus == "Control",
                              "Control",
                              SampleID)
        )
    sampleID <- meta_df$SampleID

    if (is.null(batch_column)) {
        message("Skipping Batch Adjustment: No batching column provided.")
        return(physeq)        
    } else if (batch_column %in% names(meta_df)) {
        batch <- meta_df[[batch_column]]
        batch[is.na(batch)] <- "Unknown"
        batch <- factor(batch)
        print(batch)
    } else {
        message("Skipping Batch Adjustment: No appropriate batching column found in metadata.")
        return(physeq)
    }


    # If there is uneven distribution of sample types across batches, causes ComBat to fail
    # Example: all cell controls in one processing batch, then SampleType and Batch are not linearly independent
    #           failing the assumptions of ComBats model, forced to remove covariate model
    mod <- model.matrix(~ SampleType, data = meta_df)

    if (!taxa_are_rows(physeq)) {
        otu_mat <- t(otu_mat)
    }

    # adjusted_otu_mat <- ComBat(
    #     dat = otu_mat,
    #     batch = batch,
    #     mod = mod,
    #     par.prior = TRUE
    # )
    adjusted_otu_mat <- removeBatchEffect(
        x = otu_mat,
        batch = batch,
        design = mod
    )

    otu_table(physeq) <- otu_table(adjusted_otu_mat, taxa_are_rows = TRUE)
    return(physeq)
}

average_by_techrep <- function(physeq) {
    meta_df <- as(sample_data(physeq), "data.frame")
    otu_mat <- as(otu_table(physeq), "matrix")

    if (taxa_are_rows(physeq)) {
        otu_mat <- t(otu_mat)
    }

    # Make vector for grouping by SampleID (should be identical between technical replicates)
    meta_df$SampleID <- as.character(meta_df$SampleID)
    group <- meta_df$SampleID

    # Make explicit factor for SampleID grouping
    group_levels <- unique(group)
    group_factor <- factor(group, levels = group_levels)
    group_counts <- table(group_factor)
    

    # avergae rows within each unique Sample ID
    avg_otu_mat <- rowsum(otu_mat, group = group_factor, reorder = FALSE)
    
    avg_otu_mat <- sweep(
        avg_otu_mat,
        1,
        as.numeric(group_counts[rownames(avg_otu_mat)]),
        FUN = "/"
    )
    

    # filter meta data to representative sample rows (1 for each technical replicate pair)
    avg_meta_df <- meta_df[match(rownames(avg_otu_mat), meta_df$SampleID), , drop = FALSE]
    rownames(avg_meta_df) <- rownames(avg_otu_mat)
    
    phyloseq(
        otu_table(as.matrix(avg_otu_mat), taxa_are_rows = FALSE),
        sample_data(avg_meta_df),
        tax_table(physeq)
    )
}

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
        message("Skipping Normalization: no normalization method selected.")
        return(physeq)
    } else if (norm_method == "log2") {
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else if (norm_method == "RelAbund") {
        return(transform_sample_counts(physeq, function(x) (x + pseudocount) / (sum(x + pseudocount))))
    } else if (norm_method == "RawTSS") {
        return(otu_divide_by_sample_factor(physeq, "Raw_reads"))
    } else if (norm_method == "HostMapped") {
        return(otu_divide_by_sample_factor(physeq, "Host_mapped_reads"))
    } else if (norm_method == "log2HostMapped") {
        physeq <- otu_divide_by_sample_factor(physeq, "Host_mapped_reads")
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else if (norm_method == "log2RelAbund") {
        physeq <- transform_sample_counts(physeq, function(x) (x + pseudocount) / (sum(x + pseudocount)))
        return(transform_sample_counts(physeq, function(x) log2(x)))
    } else {
        message("Skipping Normalization: Unknown normalization method provided.")
        return(physeq)
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

if (!dir.exists(paste0(args$out, "/AlphaDiversity"))) {
    dir.create(paste0(args$out, "/AlphaDiversity"),
               recursive = TRUE)
}
#-----------------------------------
# ALPHA DIVERSITY
#-----------------------------------
alpha_div_results <- estimate_richness(CompPhyseq, measures = c("Shannon", "Simpson"))
sample_data(CompPhyseq)$Shannon <- alpha_div_results$Shannon
sample_data(CompPhyseq)$Simpson <- alpha_div_results$Simpson


# Long-format data
alpha_df <- as(sample_data(CompPhyseq), "data.frame") |>
  dplyr::select(SampleType, Shannon, Simpson) |>
  pivot_longer(
    cols = c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    SampleType = factor(SampleType)
  )

# Build Dunn labels separately for each metric
make_dunn_labels <- function(dat) {
  dunn_tbl <- dunnTest(Value ~ SampleType, data = dat, method = "bh")$res |>
    mutate(
      Metric = unique(dat$Metric),
      group1 = trimws(sub(" - .*", "", Comparison)),
      group2 = trimws(sub(".* - ", "", Comparison)),
      p.adj.signif = case_when(
        P.adj <= 1e-4 ~ "****",
        P.adj <= 1e-3 ~ "***",
        P.adj <= 1e-2 ~ "**",
        P.adj <= 0.05 ~ "*",
        TRUE ~ as.character(round(P.adj, digits = 3))
      )
    )

  y_max <- max(dat$Value, na.rm = TRUE)
  y_rng <- diff(range(dat$Value, na.rm = TRUE))
  if (y_rng == 0) y_rng <- 1

  dunn_tbl |>
    arrange(P.adj) |>
    mutate(y.position = y_max + y_rng * (0.08 * row_number()))
}

dunn_labels <- alpha_df |>
  group_by(Metric) |>
  group_split() |>
  lapply(make_dunn_labels) |>
  bind_rows()

# Plot
alpha_div_plot <- ggplot(alpha_df, aes(x = SampleType, y = Value, fill = SampleType)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  stat_pvalue_manual(
    dunn_labels,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    hide.ns = FALSE,
    bracket.size = 0.4,
    size = 3,
    inherit.aes = FALSE
  )

alpha_div_plot

ggsave(
  filename = paste0(args$out, "/AlphaDiversity/", args$trialID, "_AlphaDivBoxplot.png"),
  plot = alpha_div_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Glomming to desired taxa level
if (!is.null(args$tax_agg_level)) {
    GlomPhyseq <- tax_glom(CompPhyseq,
                           taxrank = args$tax_agg_level)
    taxa_names(GlomPhyseq) <- as.character(tax_table(GlomPhyseq)[, args$tax_agg_level])
} else {
    GlomPhyseq <- CompPhyseq
}

# Removing Negative Controls
FiltPhyseq <- GlomPhyseq
#FiltPhyseq <- subset_samples(GlomPhyseq, SampleType != "NegativeControl" & SampleType %in% c("Tumor", "NormalTissue", "Nontumor"))

# Normalizing OTU Counts
NormPhyseq <- counts_normalization(physeq = FiltPhyseq, 
                                   norm_method = args$norm_method, 
                                   pseudocount = args$pseudocount)

# Pruning phyloseq of taxa and samples with zero counts
NormPhyseq <- prune_samples(sample_sums(NormPhyseq) > 0, NormPhyseq)
NormPhyseq <- prune_taxa(taxa_sums(NormPhyseq) > 0, NormPhyseq)

# Batch Adjustment
AdjPhyseq <- batch_adjustment(physeq = NormPhyseq, batch_column = args$batch_adj)

if (!dir.exists(paste0(args$out, "/DistanceOrdination"))) {
    dir.create(paste0(args$out, "/DistanceOrdination"),
               recursive = TRUE)
}
#------------------------------
# PCA Plotting
#------------------------------
#average across technical replicates if available
if (args$techrep_avg) AdjPhyseq <- average_by_techrep(NormPhyseq)

AdjOTU_mat <- as(otu_table(AdjPhyseq), "matrix")
if (taxa_are_rows(AdjPhyseq)) AdjOTU_mat <- t(AdjOTU_mat)

AdjOTU_mat <- as.matrix(AdjOTU_mat)
storage.mode(AdjOTU_mat) <- "double"

# Removing taxa that have zero variance (cause issue with PCA scaling and provide no information about inter-sample variance)
sds <- apply(AdjOTU_mat, 2, sd, na.rm = TRUE)
keep <- is.finite(sds) & sds > 0
AdjOTU_mat <- AdjOTU_mat[, keep, drop = FALSE]

TechRepAvgPCA <- prcomp(AdjOTU_mat, center = TRUE, scale. = TRUE)

PCAScores <- as.data.frame(TechRepAvgPCA$x)
if ("SampleType" %in% names(sample_data(AdjPhyseq))) {
    PCAScores$SampleType <- sample_data(AdjPhyseq)$SampleType
} else {
    PCAScores$SampleType <- sample_data(AdjPhyseq)$SampleType
}


var_explained <- TechRepAvgPCA$sdev^2
var_explained <- var_explained / sum(var_explained)

pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)


TechRepAvgPCA_plot <- ggplot(PCAScores,
    aes(PC1, PC2, color = SampleType)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
        x = paste0("PC1 (", pc1_var, "%)"),
        y = paste0("PC2 (", pc2_var, "%)")
    )

ggsave(
    filename = paste0(args$out, "/", args$trialID, "_PCA.png")
)

#------------------------------
# Distance Ordination Plotting
#------------------------------
if (is.null(args$dist_metrics)) {
    args$dist_metrics <- distanceMethodList$vegdist
}

for (metric in args$dist_metrics) {
    NormOrd <- ordinate(NormPhyseq, "PCoA", metric)
    SampleTypeOrdPlot <- plot_ordination(NormPhyseq, NormOrd, 
                                        type="samples",
                                        color= "SampleType") +
                        geom_point(size = 5) +
                        ggtitle(paste0("Composite Batch Sample Ordination (", metric, ")"))
    Centroids <- SampleTypeOrdPlot$data |>
        group_by(SampleType) |>
        summarize(
            Axis.1 = mean(Axis.1),
            Axis.2 = mean(Axis.2)
        )


    plot_data <- SampleTypeOrdPlot$data

    plot_data_with_centroids <- plot_data %>%
    left_join(
        Centroids %>% rename(Centroid1 = Axis.1, Centroid2 = Axis.2),
        by = "SampleType"
    )

    SampleTypeOrdPlot <- SampleTypeOrdPlot +
        geom_segment(
            data = plot_data_with_centroids,
            aes(x = Axis.1, y = Axis.2,
                xend = Centroid1, yend = Centroid2,
                color = SampleType),
            alpha = 0.5,
            linewidth = 0.5,
            inherit.aes = FALSE
        ) +
        # stat_ellipse(
        #     aes(x = Axis.1, y = Axis.2, color = SampleType, group = SampleType),
        #     type = "t",
        #     linetype = 2,
        #     linewidth = 1,
        #     inherit.aes = FALSE
        # ) +
        geom_point(
            data = Centroids,
            aes(x = Axis.1, y = Axis.2, color = SampleType),
            size = 6, shape = 9, stroke = 2,
            inherit.aes = FALSE
        ) +
        geom_text(
            data = Centroids,
            aes(x = Axis.1, y = Axis.2, label = SampleType),
            inherit.aes = FALSE,
            vjust = -1)
    ggsave(
        filename = paste0(args$out, "/DistanceOrdination/", args$trialID, "_", metric, "_SampleOrdination.png")
    )
}
