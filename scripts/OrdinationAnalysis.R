library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)
library(phyloseq)
library(vegan)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

parser$add_argument("--trialID",
                    type = "character",
                    default = NULL,
                    help = "ID to attach to output files")
parser$add_argument("--physeqs",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "List of phyloseq RData files named 'physeq'")
parser$add_argument("--compiled-physeq",
                    type = "character",
                    default = NULL,
                    help = "Compiled phyloseq RData file")
parser$add_argument("--batch-table",
                    type = "character",
                    default = NULL,
                    help = "Canonical batch table with include_analysis column")
parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project and ordination settings")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Normalization method")
parser$add_argument("--batch-adj",
                    type = "character",
                    default = NULL,
                    help = "Metadata column to use for batch adjustment")
parser$add_argument("--techrep-avg",
                    action = "store_true",
                    default = FALSE,
                    help = "Average technical replicates before PCA")
parser$add_argument("--dist-metric",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Distance metrics for PCoA plots")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "Pseudocount for normalization")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level to agglomerate to")
parser$add_argument("--out",
                    type = "character",
                    default = NULL,
                    help = "Output directory")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")

args <- parser$parse_args()

project_config <- list()
ordination_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project
  ordination_config <- cfg$ordination
}

trial_id <- if (!is.null(args$trialID)) {
  args$trialID
} else if (!is.null(ordination_config$trialID)) {
  ordination_config$trialID
} else {
  "analysis"
}

base_dir <- if (!is.null(project_config$base_dir)) project_config$base_dir else args$base_dir
out_dir <- if (!is.null(args$out)) {
  args$out
} else if (!is.null(project_config$output_dir)) {
  project_config$output_dir
} else {
  stop("Provide --out or project.output_dir in --analysis-config.", call. = FALSE)
}

norm_method <- if (!identical(args$norm_method, "noNorm")) args$norm_method else ordination_config$norm_method %||% args$norm_method
pseudocount <- if (!is.null(ordination_config$pseudocount)) ordination_config$pseudocount else args$pseudocount
tax_agg_level <- if (!identical(args$tax_agg_level, "Genus")) args$tax_agg_level else ordination_config$tax_agg_level %||% args$tax_agg_level
batch_adj <- if (!is.null(args$batch_adj)) args$batch_adj else ordination_config$batch_adj
dist_metrics <- if (!is.null(args$dist_metric)) args$dist_metric else ordination_config$dist_metrics
techrep_avg <- isTRUE(args$techrep_avg) || truthy_flag(ordination_config$techrep_avg, default = FALSE)

load_input_physeq <- function() {
  if (!is.null(args$compiled_physeq)) {
    return(load_physeq(args$compiled_physeq))
  }
  if (!is.null(project_config$compiled_physeq) && file.exists(project_config$compiled_physeq)) {
    return(load_physeq(project_config$compiled_physeq))
  }
  if (!is.null(args$physeqs)) {
    return(merge_physeqs(load_physeqs(args$physeqs)))
  }
  if (!is.null(args$batch_table) || !is.null(project_config$batch_table)) {
    batch_table <- if (!is.null(args$batch_table)) args$batch_table else project_config$batch_table
    physeq_paths <- resolve_batch_physeqs(batch_table, base_dir = base_dir)
    return(merge_physeqs(load_physeqs(physeq_paths)))
  }
  stop("Provide --compiled-physeq, --physeqs, --batch-table, or project.batch_table.", call. = FALSE)
}

CompPhyseq <- load_input_physeq()

dir.create(file.path(out_dir, "AlphaDiversity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "DistanceOrdination"), recursive = TRUE, showWarnings = FALSE)

alpha_div_results <- estimate_richness(CompPhyseq, measures = c("Shannon", "Simpson"))
sample_data(CompPhyseq)$Shannon <- alpha_div_results$Shannon
sample_data(CompPhyseq)$Simpson <- alpha_div_results$Simpson

alpha_df <- as(sample_data(CompPhyseq), "data.frame") |>
  dplyr::select(SampleType, Shannon, Simpson) |>
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Metric", values_to = "Value") |>
  mutate(SampleType = factor(SampleType))

make_dunn_labels <- function(dat) {
  dunn_tbl <- FSA::dunnTest(Value ~ SampleType, data = dat, method = "bh")$res |>
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
  dunn_tbl |> arrange(P.adj) |> mutate(y.position = y_max + y_rng * (0.08 * row_number()))
}

alpha_div_plot <- ggplot(alpha_df, aes(x = SampleType, y = Value, fill = SampleType)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(legend.position = "none")

if (requireNamespace("FSA", quietly = TRUE) && requireNamespace("ggpubr", quietly = TRUE)) {
  dunn_labels <- alpha_df |>
    group_by(Metric) |>
    group_split() |>
    lapply(make_dunn_labels) |>
    bind_rows()

  alpha_div_plot <- alpha_div_plot +
    ggpubr::stat_pvalue_manual(
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
} else {
  message("Skipping alpha-diversity Dunn labels because FSA and/or ggpubr is unavailable.")
}

ggsave(
  filename = file.path(out_dir, "AlphaDiversity", paste0(trial_id, "_AlphaDivBoxplot.png")),
  plot = alpha_div_plot,
  width = 10,
  height = 8,
  dpi = 300
)

NormPhyseq <- CompPhyseq |>
  tax_glom_rename(tax_agg_level) |>
  counts_normalization(norm_method = norm_method, pseudocount = pseudocount) |>
  prune_empty_physeq()

AdjPhyseq <- batch_adjustment(
  physeq = NormPhyseq,
  batch_column = batch_adj,
  method = "removeBatchEffect"
)

if (techrep_avg) {
  AdjPhyseq <- average_by_techrep(AdjPhyseq)
}

AdjOTU_mat <- otu_samples_by_taxa(AdjPhyseq)
storage.mode(AdjOTU_mat) <- "double"
sds <- apply(AdjOTU_mat, 2, sd, na.rm = TRUE)
AdjOTU_mat <- AdjOTU_mat[, is.finite(sds) & sds > 0, drop = FALSE]

if (ncol(AdjOTU_mat) >= 2 && nrow(AdjOTU_mat) >= 2) {
  pca <- prcomp(AdjOTU_mat, center = TRUE, scale. = TRUE)
  pca_scores <- as.data.frame(pca$x)
  pca_scores$SampleType <- sample_data(AdjPhyseq)$SampleType

  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  pca_plot <- ggplot(pca_scores, aes(PC1, PC2, color = SampleType)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      x = paste0("PC1 (", round(var_explained[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(var_explained[2] * 100, 2), "%)")
    )

  ggsave(
    filename = file.path(out_dir, paste0(trial_id, "_PCA.png")),
    plot = pca_plot
  )
}

if (is.null(dist_metrics) || length(dist_metrics) == 0) {
  dist_metrics <- distanceMethodList$vegdist
}

for (metric in dist_metrics) {
  ord <- ordinate(NormPhyseq, "PCoA", metric)
  ord_plot <- plot_ordination(NormPhyseq, ord, type = "samples", color = "SampleType") +
    geom_point(size = 5) +
    ggtitle(paste0("Composite Batch Sample Ordination (", metric, ")"))

  centroids <- ord_plot$data |>
    group_by(SampleType) |>
    summarize(Axis.1 = mean(Axis.1), Axis.2 = mean(Axis.2), .groups = "drop")

  plot_data <- ord_plot$data |>
    left_join(centroids |> rename(Centroid1 = Axis.1, Centroid2 = Axis.2), by = "SampleType")

  ord_plot <- ord_plot +
    geom_segment(
      data = plot_data,
      aes(x = Axis.1, y = Axis.2, xend = Centroid1, yend = Centroid2, color = SampleType),
      alpha = 0.5,
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(x = Axis.1, y = Axis.2, color = SampleType),
      size = 6,
      shape = 9,
      stroke = 2,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = centroids,
      aes(x = Axis.1, y = Axis.2, label = SampleType),
      inherit.aes = FALSE,
      vjust = -1
    )

  ggsave(
    filename = file.path(out_dir, "DistanceOrdination", paste0(trial_id, "_", metric, "_SampleOrdination.png")),
    plot = ord_plot
  )
}
