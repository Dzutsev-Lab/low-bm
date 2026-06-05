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
parser$add_argument("--group-var",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Metadata grouping variable(s) for alpha diversity, PCA, and PCoA; overrides YAML groupings")
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
  project_config <- cfg$project %||% list()
  ordination_config <- cfg$ordination %||% list()
}

trial_id <- if (!is.null(args$trialID)) {
  args$trialID
} else if (!is.null(ordination_config$trialID)) {
  ordination_config$trialID
} else {
  "analysis"
}

base_dir <- project_config$base_dir %||% args$base_dir
out_dir <- if (!is.null(args$out)) {
  args$out
} else if (!is.null(project_config$output_dir)) {
  project_config$output_dir
} else {
  stop("Provide --out or project.output_dir in --analysis-config.", call. = FALSE)
}

norm_method <- if (!identical(args$norm_method, "noNorm")) args$norm_method else ordination_config$norm_method %||% args$norm_method
pseudocount <- ordination_config$pseudocount %||% args$pseudocount
tax_agg_level <- if (!identical(args$tax_agg_level, "Genus")) args$tax_agg_level else ordination_config$tax_agg_level %||% args$tax_agg_level
batch_adj <- if (!is.null(args$batch_adj)) args$batch_adj else ordination_config$batch_adj
dist_metrics <- if (!is.null(args$dist_metric)) args$dist_metric else ordination_config$dist_metrics
techrep_avg <- isTRUE(args$techrep_avg) || truthy_flag(ordination_config$techrep_avg, default = FALSE)

is_missing_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

sanitize_path_component <- function(x, fallback = "group") {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    x <- fallback
  }
  x <- trimws(x)
  x <- gsub("[[:space:]/\\\\:;,*?\"<>|]+", "_", x)
  x <- gsub("[^A-Za-z0-9._+=@-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (!nzchar(x)) fallback else x
}

normalize_groupings <- function(config_groupings, cli_group_vars = NULL) {
  if (!is.null(cli_group_vars) && length(cli_group_vars) > 0) {
    return(lapply(as.character(cli_group_vars), function(group_var) {
      list(name = group_var, variable = group_var, sample_filter = NULL)
    }))
  }

  if (is.null(config_groupings) || length(config_groupings) == 0) {
    return(list(list(name = "SampleType", variable = "SampleType", sample_filter = NULL)))
  }

  lapply(config_groupings, function(spec) {
    if (is.character(spec)) {
      variable <- as.character(spec[[1]])
      return(list(name = variable, variable = variable, sample_filter = NULL))
    }

    variable <- spec$variable %||% spec$name
    if (is_missing_value(variable)) {
      stop("Each ordination grouping must define a variable.", call. = FALSE)
    }

    list(
      name = spec$name %||% variable,
      variable = variable,
      sample_filter = spec$sample_filter
    )
  })
}

prepare_grouped_physeq <- function(physeq, grouping, data_label) {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  variable <- as.character(grouping$variable)
  group_name <- as.character(grouping$name)
  fail_missing_columns(names(metadata_df), variable, "phyloseq sample_data")

  if (!is.null(grouping$sample_filter) && length(grouping$sample_filter) > 0) {
    physeq <- apply_sample_filter(physeq, grouping$sample_filter)
    metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  }

  group_values <- as.character(metadata_df[[variable]])
  valid_group <- !is.na(group_values) & nzchar(trimws(group_values))
  if (any(!valid_group)) {
    warning(
      "Dropping ",
      sum(!valid_group),
      " sample(s) from ",
      data_label,
      " grouping '",
      group_name,
      "' because ",
      variable,
      " is missing or blank.",
      call. = FALSE
    )
    physeq <- phyloseq::prune_samples(valid_group, physeq)
    metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
    group_values <- as.character(metadata_df[[variable]])
  }

  if (phyloseq::nsamples(physeq) == 0) {
    warning("Skipping ", data_label, " grouping '", group_name, "': no samples remain.", call. = FALSE)
    return(NULL)
  }

  group_values <- factor(group_values)
  phyloseq::sample_data(physeq)$OrdinationGroup <- group_values
  physeq
}

has_enough_groups <- function(physeq, grouping, data_label) {
  if (is.null(physeq)) {
    return(FALSE)
  }

  groups <- phyloseq::sample_data(physeq)$OrdinationGroup
  n_groups <- length(unique(as.character(groups[!is.na(groups)])))
  if (phyloseq::nsamples(physeq) < 2 || n_groups < 2) {
    warning(
      "Skipping ",
      data_label,
      " grouping '",
      grouping$name,
      "': fewer than two samples or two groups remain.",
      call. = FALSE
    )
    return(FALSE)
  }

  TRUE
}

make_group_dirs <- function(out_dir, group_path) {
  dirs <- list(
    alpha = file.path(out_dir, "AlphaDiversity", group_path),
    pca = file.path(out_dir, "PCA", group_path),
    ordination = file.path(out_dir, "DistanceOrdination", group_path)
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

make_dunn_labels <- function(dat) {
  dunn_tbl <- FSA::dunnTest(Value ~ OrdinationGroup, data = dat, method = "bh")$res |>
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

plot_alpha_diversity <- function(alpha_physeq, grouping, group_path, group_dirs) {
  grouped_physeq <- prepare_grouped_physeq(alpha_physeq, grouping, "alpha diversity")
  if (!has_enough_groups(grouped_physeq, grouping, "alpha diversity")) {
    return(invisible(NULL))
  }

  variable <- as.character(grouping$variable)
  group_name <- as.character(grouping$name)
  alpha_df <- as(phyloseq::sample_data(grouped_physeq), "data.frame") |>
    dplyr::select(OrdinationGroup, Shannon, Simpson) |>
    pivot_longer(cols = c(Shannon, Simpson), names_to = "Metric", values_to = "Value") |>
    mutate(OrdinationGroup = factor(OrdinationGroup))

  alpha_div_plot <- ggplot(alpha_df, aes(x = OrdinationGroup, y = Value, fill = OrdinationGroup)) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    facet_wrap(~Metric, scales = "free_y", nrow = 2) +
    labs(x = variable, fill = variable, title = paste0("Alpha diversity by ", group_name)) +
    theme_bw() +
    theme(legend.position = "none")

  if (requireNamespace("FSA", quietly = TRUE) && requireNamespace("ggpubr", quietly = TRUE)) {
    dunn_labels <- tryCatch(
      alpha_df |>
        group_by(Metric) |>
        group_split() |>
        lapply(make_dunn_labels) |>
        bind_rows(),
      error = function(e) {
        message(
          "Skipping alpha-diversity Dunn labels for grouping '",
          group_name,
          "': ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(dunn_labels) && nrow(dunn_labels) > 0) {
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
    }
  } else {
    message("Skipping alpha-diversity Dunn labels because FSA and/or ggpubr is unavailable.")
  }

  ggsave(
    filename = file.path(group_dirs$alpha, paste0(trial_id, "_", group_path, "_AlphaDivBoxplot.png")),
    plot = alpha_div_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

plot_pca <- function(adj_physeq, grouping, group_path, group_dirs) {
  grouped_physeq <- prepare_grouped_physeq(adj_physeq, grouping, "PCA")
  if (!has_enough_groups(grouped_physeq, grouping, "PCA")) {
    return(invisible(NULL))
  }

  otu_mat <- otu_samples_by_taxa(grouped_physeq)
  storage.mode(otu_mat) <- "double"
  sds <- apply(otu_mat, 2, sd, na.rm = TRUE)
  otu_mat <- otu_mat[, is.finite(sds) & sds > 0, drop = FALSE]

  if (ncol(otu_mat) < 2 || nrow(otu_mat) < 2) {
    warning("Skipping PCA grouping '", grouping$name, "': not enough variable taxa or samples remain.", call. = FALSE)
    return(invisible(NULL))
  }

  pca <- prcomp(otu_mat, center = TRUE, scale. = TRUE)
  pca_scores <- as.data.frame(pca$x)
  metadata_df <- as.data.frame(phyloseq::sample_data(grouped_physeq), stringsAsFactors = FALSE)
  group_lookup <- setNames(as.character(metadata_df$OrdinationGroup), rownames(metadata_df))
  pca_scores$OrdinationGroup <- factor(group_lookup[rownames(pca_scores)])

  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  pca_plot <- ggplot(pca_scores, aes(PC1, PC2, color = OrdinationGroup)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste0("PCA colored by ", grouping$name),
      color = grouping$variable,
      x = paste0("PC1 (", round(var_explained[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(var_explained[2] * 100, 2), "%)")
    )

  ggsave(
    filename = file.path(group_dirs$pca, paste0(trial_id, "_", group_path, "_PCA.png")),
    plot = pca_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

plot_distance_ordination <- function(norm_physeq, grouping, group_path, group_dirs, dist_metrics) {
  grouped_physeq <- prepare_grouped_physeq(norm_physeq, grouping, "distance ordination")
  if (!has_enough_groups(grouped_physeq, grouping, "distance ordination")) {
    return(invisible(NULL))
  }

  for (metric in dist_metrics) {
    ord <- tryCatch(
      ordinate(grouped_physeq, "PCoA", metric),
      error = function(e) {
        warning(
          "Skipping distance ordination grouping '",
          grouping$name,
          "' metric '",
          metric,
          "': ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (is.null(ord)) {
      next
    }

    ord_plot <- plot_ordination(grouped_physeq, ord, type = "samples", color = "OrdinationGroup") +
      geom_point(size = 5) +
      ggtitle(paste0("Sample ordination by ", grouping$name, " (", metric, ")")) +
      labs(color = grouping$variable)

    centroids <- ord_plot$data |>
      group_by(OrdinationGroup) |>
      summarize(Axis.1 = mean(Axis.1), Axis.2 = mean(Axis.2), .groups = "drop")

    plot_data <- ord_plot$data |>
      left_join(centroids |> rename(Centroid1 = Axis.1, Centroid2 = Axis.2), by = "OrdinationGroup")

    ord_plot <- ord_plot +
      geom_segment(
        data = plot_data,
        aes(x = Axis.1, y = Axis.2, xend = Centroid1, yend = Centroid2, color = OrdinationGroup),
        alpha = 0.5,
        linewidth = 0.5,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = centroids,
        aes(x = Axis.1, y = Axis.2, color = OrdinationGroup),
        size = 6,
        shape = 9,
        stroke = 2,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = centroids,
        aes(x = Axis.1, y = Axis.2, label = OrdinationGroup),
        inherit.aes = FALSE,
        vjust = -1
      )

    ggsave(
      filename = file.path(group_dirs$ordination, paste0(trial_id, "_", metric, "_", group_path, "_SampleOrdination.png")),
      plot = ord_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

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

groupings <- normalize_groupings(ordination_config$groupings, args$group_var)

if (is.null(dist_metrics) || length(dist_metrics) == 0) {
  dist_metrics <- distanceMethodList$vegdist
}

CompPhyseq <- load_input_physeq()

alpha_div_results <- estimate_richness(CompPhyseq, measures = c("Shannon", "Simpson"))
phyloseq::sample_data(CompPhyseq)$Shannon <- alpha_div_results$Shannon
phyloseq::sample_data(CompPhyseq)$Simpson <- alpha_div_results$Simpson

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

group_paths <- make.unique(
  vapply(groupings, function(grouping) sanitize_path_component(grouping$name), character(1)),
  sep = "_"
)

for (i in seq_along(groupings)) {
  grouping <- groupings[[i]]
  group_path <- group_paths[[i]]
  group_dirs <- make_group_dirs(out_dir, group_path)

  message("Running ordination outputs for grouping: ", grouping$name, " (", grouping$variable, ")")
  plot_alpha_diversity(CompPhyseq, grouping, group_path, group_dirs)
  plot_pca(AdjPhyseq, grouping, group_path, group_dirs)
  plot_distance_ordination(NormPhyseq, grouping, group_path, group_dirs, dist_metrics)
}
