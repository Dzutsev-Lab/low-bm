library(argparse)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(dplyr)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with differential_abundance comparison specs")
parser$add_argument("--compiled-physeq",
                    type = "character",
                    default = NULL,
                    help = "Compiled phyloseq RData file")
parser$add_argument("--physeqs",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "List of phyloseq RData files named 'physeq'")
parser$add_argument("--batch-table",
                    type = "character",
                    default = NULL,
                    help = "Canonical batch table")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")

# Legacy-compatible arguments.
parser$add_argument("--trialID",
                    type = "character",
                    default = NULL,
                    help = "ID number to attach to output files")
parser$add_argument("--B1-physeq",
                    type = "character",
                    default = NULL,
                    help = "Legacy first-batch phyloseq RData file")
parser$add_argument("--B2-physeq",
                    type = "character",
                    default = NULL,
                    help = "Legacy optional second-batch phyloseq RData file")
parser$add_argument("--DA-methods",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Legacy method list; config-driven DA uses ANCOMBC2")
parser$add_argument("--DA-comparisons",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Legacy comparison preset names")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Normalization method for DA heatmaps")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "Pseudocount for heatmap normalization")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level to agglomerate to for DA")
parser$add_argument("--tax-label-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level to use for labels")
parser$add_argument("--select-taxa-names",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Optional files containing taxa names to analyze")
parser$add_argument("--alpha",
                    type = "double",
                    default = 0.05,
                    help = "Alpha cutoff for significance")
parser$add_argument("--lfc-cutoff",
                    type = "double",
                    default = 1.0,
                    help = "Log-fold-change cutoff for significance")
parser$add_argument("--out",
                    type = "character",
                    default = NULL,
                    help = "Output directory")

args <- parser$parse_args()

load_project_and_da_config <- function() {
  if (!is.null(args$analysis_config)) {
    full_config <- load_yaml_config(args$analysis_config)
    project_config <- full_config$project %||% list()
    da_config <- apply_project_config_defaults(
      full_config$differential_abundance,
      project_config,
      c("trialID", "output_dir", "norm_method", "pseudocount", "tax_agg_level")
    )
    return(list(
      project = project_config,
      da = normalize_da_config(da_config)
    ))
  }

  if (is.null(args$DA_comparisons)) {
    args$DA_comparisons <- c("NegativeControl", "PatientSample")
  }

  list(
    project = list(),
    da = normalize_da_config(build_legacy_da_config(
      trial_id = args$trialID %||% "analysis",
      comparisons = args$DA_comparisons,
      methods = args$DA_methods %||% c("ANCOMBC2"),
      out_dir = args$out %||% "Exp_Output/analysis",
      norm_method = args$norm_method,
      pseudocount = args$pseudocount,
      tax_agg_level = args$tax_agg_level,
      tax_label_level = args$tax_label_level,
      alpha = args$alpha,
      lfc_cutoff = args$lfc_cutoff
    ))
  )
}

config_bundle <- load_project_and_da_config()
project_config <- config_bundle$project
da_config <- config_bundle$da

trial_id <- args$trialID %||% analysis_config_value(project_config, da_config, "trialID", "analysis")
base_dir <- project_config$base_dir %||% args$base_dir
out_dir <- args$out %||% analysis_output_dir(project_config, da_config, default = file.path(base_dir, "analysis"))

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
  if (!is.null(args$B1_physeq)) {
    paths <- c(args$B1_physeq, args$B2_physeq)
    paths <- paths[!is.na(paths) & !is.null(paths) & nzchar(paths)]
    return(merge_physeqs(load_physeqs(paths)))
  }
  if (!is.null(args$batch_table) || !is.null(project_config$batch_table)) {
    batch_table <- args$batch_table %||% project_config$batch_table
    physeq_paths <- resolve_batch_physeqs(batch_table, base_dir = base_dir)
    return(merge_physeqs(load_physeqs(physeq_paths)))
  }
  stop(
    "Provide --compiled-physeq, --physeqs, --B1-physeq, --batch-table, or project.batch_table.",
    call. = FALSE
  )
}

read_select_taxa <- function(paths) {
  if (is.null(paths)) {
    return(NULL)
  }
  unique(unlist(lapply(paths, function(name_file) {
    read.csv(name_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  })))
}

write_patient_duplicate_policy_audit <- function(audit_df, out_dir, trial_id, method, comparison_name) {
  if (is.null(audit_df)) {
    audit_df <- empty_patient_duplicate_policy_audit()
  }
  comparison_dir <- file.path(out_dir, method, comparison_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(
    comparison_dir,
    paste0(trial_id, "_", comparison_name, "_PatientDuplicatePolicy.tsv")
  )
  write.table(
    audit_df,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  out_file
}

add_result_labels <- function(results_df, grouped_physeq, tax_label_level, tax_agg_level) {
  results_df$label <- results_df$taxon
  tax_df <- as.data.frame(as(tax_table(grouped_physeq), "matrix"), stringsAsFactors = FALSE)
  if ((!tax_label_level %in% names(tax_df)) | (tax_label_level == tax_agg_level)) {
    return(results_df)
  }

  tax_df$taxon <- rownames(tax_df)
  results_df <- left_join(results_df, tax_df[, c("taxon", tax_label_level), drop = FALSE], by = "taxon")
  label_values <- results_df[[tax_label_level]]
  results_df$label <- ifelse(
    is.na(label_values) | label_values == "",
    results_df$taxon,
    paste0(label_values, " (", results_df$taxon, ")")
  )
  results_df
}

plot_da_volcano <- function(results_df, spec, alpha, lfc_cutoff, out_dir, trial_id) {
  if ("test" %in% names(results_df)) {
    plot_df <- results_df |>
      filter(test %in% c("primary", "pairwise", "dunnet", "trend"), !is.na(log2FoldChange), is.na(struc0))
  } else {
    plot_df <- results_df |> filter(is.na(struc0))
  }
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  plot_title <- spec$plot_title %||% spec$name
  direction_labels <- ancombc_direction_labels(spec)
  plot_groups <- split(plot_df, interaction(plot_df$test %||% "primary", plot_df$contrast %||% "contrast", drop = TRUE))

  for (plot_name in names(plot_groups)) {
    single_plot_df <- plot_groups[[plot_name]]
    single_plot_df$volcano_p <- da_volcano_p_value(single_plot_df)
    single_plot_df <- single_plot_df |> filter(!is.na(volcano_p), is.finite(volcano_p))
    if (nrow(single_plot_df) == 0) {
      next
    }
    sig_df <- single_plot_df[da_volcano_highlight_mask(single_plot_df, lfc_cutoff), , drop = FALSE]
    contrast_label <- if ("contrast" %in% names(single_plot_df)) unique(single_plot_df$contrast)[[1]] else NULL
    test_label <- if ("test" %in% names(single_plot_df)) unique(single_plot_df$test)[[1]] else "primary"
    y_label <- if (any(is.na(single_plot_df$padj) & !is.na(single_plot_df$p))) {
      "-log10(adjusted p-value; p-value fallback)"
    } else {
      "-log10(adjusted p-value)"
    }
    volcano <- ggplot(single_plot_df, aes(x = log2FoldChange, y = -log10(volcano_p))) +
      geom_point(alpha = 0.6, size = 4, color = "grey40") +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "darkred") +
      geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkred") +
      labs(
        title = paste("Volcano Plot:", plot_title),
        subtitle = paste(na.omit(c(test_label, contrast_label)), collapse = " | "),
        x = direction_labels$x,
        y = y_label,
        caption = direction_labels$caption
      ) +
      theme_bw() +
      theme(plot.caption = element_text(hjust = 0.5, size = 12))

    if (nrow(sig_df) > 0) {
      volcano <- volcano +
        geom_point(data = sig_df, aes(color = direction), size = 5) +
        scale_color_manual(
          name = direction_labels$legend_title,
          values = c(pos = "firebrick1", neg = "dodgerblue1", none = "grey40"),
          breaks = c("neg", "pos", "none"),
          labels = direction_labels$legend_labels,
          na.translate = FALSE
        ) +
        geom_text_repel(
          data = sig_df,
          aes(label = label, color = direction),
          size = 5,
          max.overlaps = Inf,
          box.padding = 0.5,
          point.padding = 0.4,
          segment.size = 0.5,
          segment.color = "grey40",
          show.legend = FALSE
        )
    }

    suffix <- gsub("[^A-Za-z0-9._-]+", "_", paste(na.omit(c(test_label, contrast_label)), collapse = "_"))
    suffix <- gsub("_+", "_", suffix)
    ggsave(
      filename = file.path(out_dir, "ANCOMBC2", spec$name, paste0(trial_id, "_", spec$name, "_", suffix, "_ANCOMBC2Volcano.png")),
      plot = volcano,
      width = 14,
      height = 12,
      units = "in",
      dpi = 300
    )
  }

  invisible(NULL)
}

plot_da_heatmap <- function(grouped_physeq,
                            results_df,
                            spec,
                            da_config,
                            out_dir,
                            trial_id) {
  sig_taxa <- results_df$taxon[results_df$significance == "Sig"]
  sig_taxa <- intersect(sig_taxa, taxa_names(grouped_physeq))
  if (length(sig_taxa) == 0) {
    return(invisible(NULL))
  }

  norm_physeq <- counts_normalization(
    grouped_physeq,
    norm_method = da_config$norm_method,
    pseudocount = da_config$pseudocount
  )
  pruned_physeq <- prune_taxa(sig_taxa, norm_physeq)
  pruned_physeq <- prune_samples(sample_sums(pruned_physeq) > 0, pruned_physeq)

  metadata_df <- as(sample_data(pruned_physeq), "data.frame")
  order_cols <- intersect(c("SampleType", "PatientID"), names(metadata_df))
  if (length(order_cols) > 0) {
    column_order <- rownames(metadata_df)[do.call(order, metadata_df[order_cols])]
  } else {
    column_order <- rownames(metadata_df)
  }

  tax_label_level <- spec$tax_label_level %||% da_config$tax_label_level
  heatmap <- plot_heatmap(
    pruned_physeq,
    method = "PCoA",
    distance = "bray",
    sample.order = column_order,
    sample.label = "SampleID",
    taxa.label = tax_label_level,
    low = "#000033",
    high = "#FF3300",
    na.value = "black"
  )

  if ("SampleType" %in% names(metadata_df)) {
    type_ordered <- metadata_df[column_order, "SampleType"]
    heatmap <- heatmap +
      geom_vline(xintercept = cumsum(table(type_ordered)) + 0.5, linewidth = 1, color = "white")
  }

  ggsave(
    filename = file.path(out_dir, "ANCOMBC2", spec$name, paste0(trial_id, "_", spec$name, "_ANCOMBC2Heatmap.png")),
    plot = heatmap,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
  )
}

plot_ancombc2_result_heatmap <- function(results_df, spec, out_dir, trial_id) {
  if (!all(c("test", "contrast", "log2FoldChange", "significance") %in% names(results_df))) {
    return(invisible(NULL))
  }
  plot_df <- results_df |>
    filter(test %in% c("primary", "pairwise", "dunnet", "trend"), !is.na(log2FoldChange), significance == "Sig")
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }
  plot_df$display_contrast <- paste(plot_df$test, plot_df$contrast, sep = ": ")
  plot <- ggplot(plot_df, aes(x = display_contrast, y = taxon, fill = log2FoldChange)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", midpoint = 0) +
    labs(
      title = paste("ANCOM-BC2 Significant Effects:", spec$plot_title %||% spec$name),
      x = NULL,
      y = NULL,
      fill = "log2 FC"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    filename = file.path(out_dir, "ANCOMBC2", spec$name, paste0(trial_id, "_", spec$name, "_ANCOMBC2ResultHeatmap.png")),
    plot = plot,
    width = 12,
    height = max(4, min(14, 0.3 * length(unique(plot_df$taxon)) + 3)),
    units = "in",
    dpi = 300
  )
}

select_taxa <- read_select_taxa(args$select_taxa_names)
CompPhyseq <- load_input_physeq()

for (spec in da_config$comparisons) {
  message("Running ANCOMBC2 comparison: ", spec$name)
  comparison_physeq <- prepare_da_physeq(
    physeq = CompPhyseq,
    spec = spec,
    global_config = da_config,
    select_taxa = select_taxa
  )
  resolved_spec <- attr(comparison_physeq, "da_comparison_spec")
  if (is.null(resolved_spec)) {
    resolved_spec <- spec
  }
  duplicate_policy_file <- write_patient_duplicate_policy_audit(
    attr(comparison_physeq, "patient_duplicate_policy_audit"),
    out_dir = out_dir,
    trial_id = trial_id,
    method = "ANCOMBC2",
    comparison_name = resolved_spec$name
  )
  message("Wrote patient duplicate policy audit: ", duplicate_policy_file)

  alpha <- resolved_spec$alpha %||% da_config$alpha
  lfc_cutoff <- resolved_spec$lfc_cutoff %||% da_config$lfc_cutoff
  results_df <- run_ancombc_comparison(comparison_physeq, resolved_spec, da_config)
  results_for_plots <- add_result_labels(
    results_df,
    comparison_physeq,
    resolved_spec$tax_label_level %||% da_config$tax_label_level,
    resolved_spec$tax_agg_level %||% da_config$tax_agg_level
  )

  result_files <- write_ancombc2_results(
    results_df,
    out_dir = out_dir,
    trial_id = trial_id,
    comparison_name = resolved_spec$name
  )
  message("Wrote DA results: ", paste(result_files, collapse = ", "))

  plot_da_volcano(results_for_plots, resolved_spec, alpha, lfc_cutoff, out_dir, trial_id)
  plot_ancombc2_result_heatmap(results_for_plots, resolved_spec, out_dir, trial_id)
  plot_da_heatmap(comparison_physeq, results_for_plots, resolved_spec, da_config, out_dir, trial_id)
}
