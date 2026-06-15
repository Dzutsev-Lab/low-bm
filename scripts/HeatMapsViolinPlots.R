library(phyloseq)
library(dplyr)
library(ggplot2)
library(argparse)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

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
                    help = "Analysis YAML with project and heatmap_violin settings")
parser$add_argument("--patient-sample-batches",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Sequencing batches with patient sample records")
parser$add_argument("--DA-results",
                    type = "character",
                    help = "DA results TSV, either absolute or relative to Exp_Output")
parser$add_argument("--select-taxa-names",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Optional files containing taxa to plot")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level to agglomerate to")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Normalization method")
parser$add_argument("--batch-adj",
                    action = "store_true",
                    default = FALSE,
                    help = "Adjust counts by ProcessingBatch or Batch if available")
parser$add_argument("--limma-voom",
                    action = "store_true",
                    default = FALSE,
                    help = "Normalize counts using TMM and voom")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "Pseudocount for normalization")
parser$add_argument("--out-dir",
                    type = "character",
                    default = NULL,
                    help = "Output directory within Exp_Output")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")

args <- parser$parse_args()

project_config <- list()
plot_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project %||% list()
  plot_config <- cfg$heatmap_violin %||% list()
}

base_dir <- project_config$base_dir %||% args$base_dir
out_root <- if (!is.null(args$out_dir)) args$out_dir else analysis_output_dir(project_config, plot_config, default = "analysis")
out_dir <- resolve_output_path(out_root, base_dir = base_dir)

tax_agg_level <- if (!identical(args$tax_agg_level, "Genus")) args$tax_agg_level else analysis_config_value(project_config, plot_config, "tax_agg_level", args$tax_agg_level)
norm_method <- if (!identical(args$norm_method, "noNorm")) args$norm_method else analysis_config_value(project_config, plot_config, "norm_method", args$norm_method)
pseudocount <- analysis_config_value(project_config, plot_config, "pseudocount", args$pseudocount)
limma_voom <- isTRUE(args$limma_voom) || truthy_flag(plot_config$limma_voom, default = FALSE)
batch_adj_requested <- isTRUE(args$batch_adj) || truthy_flag(plot_config$batch_adj, default = FALSE)
patient_sample_batches <- args$patient_sample_batches %||% plot_config$patient_sample_batches

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
    batch_table <- args$batch_table %||% project_config$batch_table
    physeq_paths <- resolve_batch_physeqs(batch_table, base_dir = base_dir)
    return(merge_physeqs(load_physeqs(physeq_paths)))
  }
  stop("Provide --compiled-physeq, --physeqs, --batch-table, or project.batch_table.", call. = FALSE)
}

resolve_da_results <- function(path) {
  if (is.null(path)) {
    stop("Provide --DA-results for heatmap/violin plotting.", call. = FALSE)
  }
  if (file.exists(path)) {
    return(path)
  }
  candidate <- file.path(base_dir, path)
  if (file.exists(candidate)) {
    return(candidate)
  }
  stop("DA results file not found: ", path, call. = FALSE)
}

choose_batch_column <- function(physeq) {
  metadata_df <- as.data.frame(sample_data(physeq), stringsAsFactors = FALSE)
  if ("ProcessingBatch" %in% names(metadata_df)) return("ProcessingBatch")
  if ("Batch" %in% names(metadata_df)) return("Batch")
  NULL
}

CompPhyseq <- load_input_physeq()

GlomPhyseq <- tax_glom_rename(CompPhyseq, tax_agg_level)
FiltPhyseq <- subset_samples(GlomPhyseq, SampleType != "NegativeControl")

if (limma_voom) {
  norm_method <- "noNorm"
  NormPhyseq <- limma_voom_normalization(FiltPhyseq)
} else {
  NormPhyseq <- counts_normalization(
    physeq = FiltPhyseq,
    norm_method = norm_method,
    pseudocount = pseudocount
  )
}

if (batch_adj_requested) {
  NormPhyseq <- batch_adjustment(
    physeq = NormPhyseq,
    batch_column = choose_batch_column(NormPhyseq),
    method = "ComBat"
  )
}

DA_results_df <- read.delim(resolve_da_results(args$DA_results), header = TRUE)
if (!is.null(args$select_taxa_names)) {
  taxa_of_interest <- unique(unlist(lapply(args$select_taxa_names, function(name_file) {
    read.csv(name_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  })))
} else {
  taxa_of_interest <- DA_results_df$taxon[DA_results_df$significance == "Sig"]
}

if (length(taxa_of_interest) == 0) {
  stop("No taxa selected for heatmap/violin plotting.", call. = FALSE)
}

MeltPhyseq_df <- psmelt(NormPhyseq)
if (!"SequencingBatch" %in% names(MeltPhyseq_df)) {
  MeltPhyseq_df$SequencingBatch <- "Unknown"
}
MeltPhyseq_df <- MeltPhyseq_df |>
  mutate(
    SequencingBatch = factor(SequencingBatch),
    SampleType = factor(SampleType)
  )

multibatch_comparison <- length(patient_sample_batches) > 1
plot_out_dir <- file.path(out_dir, norm_method)

for (taxon in taxa_of_interest) {
  taxon_dir <- file.path(plot_out_dir, taxon)
  dir.create(taxon_dir, recursive = TRUE, showWarnings = FALSE)

  single_tax_df <- MeltPhyseq_df |> filter(OTU == taxon)
  single_tax_DA_results_df <- DA_results_df[DA_results_df$taxon == taxon, , drop = FALSE]

  if (nrow(single_tax_df) == 0) {
    warning("Skipping taxon absent from melted phyloseq object: ", taxon)
    next
  }

  if (multibatch_comparison && all(c("logFC_b1", "logFC_b2", "adj_p_b1", "adj_p_b2") %in% names(single_tax_DA_results_df))) {
    techrep_avg_single_tax_df <- single_tax_df |>
      group_by(SampleID) |>
      summarise(
        Abundance = mean(Abundance, na.rm = FALSE),
        PatientOutcome = if ("PatientOutcome" %in% names(single_tax_df)) first(PatientOutcome) else NA,
        SampleType = first(SampleType),
        PatientID = first(PatientID),
        .groups = "drop"
      )

    b1_lfc <- round(single_tax_DA_results_df[, "logFC_b1"], 2)
    b2_lfc <- round(single_tax_DA_results_df[, "logFC_b2"], 2)
    b1_adj_p <- round(single_tax_DA_results_df[, "adj_p_b1"], 3)
    b2_adj_p <- round(single_tax_DA_results_df[, "adj_p_b2"], 3)
    het_q <- if ("het_Q" %in% names(single_tax_DA_results_df)) round(single_tax_DA_results_df[, "het_Q"], 3) else NA
    het_q_adj <- if ("het_Q_q" %in% names(single_tax_DA_results_df)) round(single_tax_DA_results_df[, "het_Q_q"], 3) else NA
    fisher_q <- if ("fisher_q" %in% names(single_tax_DA_results_df)) round(single_tax_DA_results_df[, "fisher_q"], 4) else NA

    violin_plot <- ggplot(single_tax_df, aes(x = SampleType, y = Abundance)) +
      geom_violin(aes(fill = SampleType, group = SampleType), trim = FALSE, alpha = 0.6) +
      geom_jitter(aes(fill = SampleType, shape = SequencingBatch), width = 0.2, size = 1.5, alpha = 0.8) +
      labs(
        title = paste("Abundance Comparison for", taxon),
        subtitle = paste0(
          "LFC Exp1: ", b1_lfc, " (q=", b1_adj_p, "), ",
          "LFC Exp2: ", b2_lfc, " (q=", b2_adj_p, "); ",
          "HetQ: ", het_q, " (q=", het_q_adj, ")"
        ),
        x = "Sample Type",
        y = paste0("Normalized Abundance (", norm_method, ")")
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    techrep_avg_violin_plot <- ggplot(
      filter(techrep_avg_single_tax_df, SampleType != "CellLineControl"),
      aes(x = SampleType, y = Abundance)
    ) +
      geom_violin(aes(fill = SampleType, group = SampleType), trim = FALSE, alpha = 0.6) +
      geom_jitter(aes(fill = SampleType), width = 0.2, size = 1.5, alpha = 0.8) +
      labs(
        title = paste("Averaged Abundance Comparison for", taxon),
        subtitle = paste0("Fisher q-value: ", fisher_q, "; HetQ: ", het_q, " (q=", het_q_adj, ")"),
        x = "Sample Type",
        y = paste0("Normalized Avg Abund (", norm_method, ")")
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(
      filename = file.path(taxon_dir, paste0(taxon, "_", norm_method, "_techrep_avg_violin_plot.png")),
      plot = techrep_avg_violin_plot,
      width = 8,
      height = 4
    )
  } else {
    lfc <- if ("log2FoldChange" %in% names(single_tax_DA_results_df)) round(single_tax_DA_results_df[, "log2FoldChange"], 2) else NA
    adj_p <- if ("padj" %in% names(single_tax_DA_results_df)) round(single_tax_DA_results_df[, "padj"], 4) else NA

    violin_plot <- ggplot(single_tax_df, aes(x = SampleType, y = Abundance, fill = SampleType)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
      labs(
        title = paste("Abundance Comparison for", taxon),
        subtitle = paste0("LFC = ", lfc, " (adj p-value = ", adj_p, ")"),
        x = "Sample Type",
        y = paste0("Normalized Abundance (", norm_method, ")")
      ) +
      theme_bw()
  }

  ggsave(
    filename = file.path(taxon_dir, paste0(taxon, "_", norm_method, "_violin_plot.png")),
    plot = violin_plot,
    width = 8,
    height = 4
  )
}

heat_mapping <- function(patient_batch_name, subset_physeq) {
  pruned_physeq <- prune_taxa(taxa_of_interest, subset_physeq)
  pruned_physeq <- prune_samples(sample_sums(pruned_physeq) > 0, pruned_physeq)
  column_order <- rownames(sample_data(pruned_physeq))[order(
    sample_data(pruned_physeq)$SampleType,
    sample_data(pruned_physeq)$PatientID
  )]

  heatmap_plot <- plot_heatmap(
    pruned_physeq,
    sample.order = column_order,
    method = "PCoA",
    distance = "bray",
    title = paste0(
      patient_batch_name,
      " Select Taxa Heatmap ",
      ifelse(norm_method == "noNorm", "(Unnormalized)", paste0("(", norm_method, ")"))
    ),
    sample.label = "SampleID",
    taxa.label = tax_agg_level
  )

  type_ordered <- sample_data(pruned_physeq)[column_order, "SampleType"]
  heatmap_plot <- heatmap_plot +
    geom_vline(xintercept = cumsum(table(type_ordered)) + 0.5, linewidth = 1, color = "white")

  ggsave(
    filename = file.path(plot_out_dir, paste0(patient_batch_name, "_", norm_method, "_SelectTaxaHeatmap.png")),
    plot = heatmap_plot,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
  )
}

if (multibatch_comparison) {
  for (patient_batch_name in patient_sample_batches) {
    subset_physeq <- subset_samples(
      NormPhyseq,
      SequencingBatch == patient_batch_name | SampleType == "CellLineControl"
    )
    heat_mapping(patient_batch_name, subset_physeq)
  }
} else {
  subset_physeq <- subset_samples(
    NormPhyseq,
    SampleType %in% c("Tumor", "Nontumor", "NormalTissue", "CellLineControl")
  )
  heat_mapping("TechRep-Averaged", subset_physeq)
}
