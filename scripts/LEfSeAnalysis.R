library(argparse)
library(phyloseq)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))
source(file.path("scripts", "Rhelpers", "LEfSeAnalysis.R"))

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    required = TRUE,
                    help = "Analysis YAML with project, differential_abundance, and lefse_analysis settings")
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
                    help = "Canonical batch table with include_analysis column")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")
parser$add_argument("--out",
                    type = "character",
                    default = NULL,
                    help = "Output directory")
parser$add_argument("--trialID",
                    type = "character",
                    default = NULL,
                    help = "ID to attach to output files")

args <- parser$parse_args()

if (!requireNamespace("lefser", quietly = TRUE)) {
  stop(
    "The R package 'lefser' is required for LEfSe analysis. ",
    "Install it with BiocManager::install('lefser') in the R environment used by this pipeline.",
    call. = FALSE
  )
}

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
da_config <- normalize_da_config(apply_project_config_defaults(
  cfg$differential_abundance,
  project_config,
  c("trialID", "output_dir", "norm_method", "pseudocount", "tax_agg_level")
))
lefse_config <- normalize_lefse_config(cfg$lefse_analysis %||% list(), project_config)

trial_id <- args$trialID %||% analysis_config_value(project_config, lefse_config, "trialID", "analysis")
base_dir <- project_config$base_dir %||% args$base_dir
out_dir <- args$out %||% analysis_output_dir(project_config, lefse_config, default = file.path(base_dir, "analysis"))

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
  stop(
    "Provide --compiled-physeq, --physeqs, --batch-table, or project.batch_table.",
    call. = FALSE
  )
}

version_or_missing <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    "not installed"
  }
}

make_run_note <- function(spec,
                          status,
                          prepared = NULL,
                          result_file = NULL,
                          plot_file = NULL,
                          skip_reason = NULL) {
  norm_method <- if (is.null(prepared)) resolve_lefse_norm_method(lefse_config, project_config) else prepared$norm_method
  pseudocount <- if (is.null(prepared)) resolve_lefse_pseudocount(lefse_config, project_config) else prepared$pseudocount
  class_levels <- if (is.null(prepared)) character(0) else prepared$class_levels
  class_col <- if (is.null(prepared)) comparison_class_col(spec) else prepared$class_col
  physeq <- prepared$physeq %||% NULL

  lines <- c(
    paste0("Comparison: ", spec$name),
    paste0("Status: ", status),
    paste0("Plot title: ", spec$plot_title %||% spec$name),
    paste0("Class column: ", class_col),
    paste0("Class levels: ", if (length(class_levels) == 0) "not available" else paste(class_levels, collapse = ", ")),
    paste0("Subclass column: ", spec$subclass_col %||% lefse_config$subclass_col %||% "none"),
    paste0("Taxonomic aggregation level: ", spec$tax_agg_level %||% lefse_config$tax_agg_level),
    paste0("Abundance scale: ", lefse_config$abundance_scale),
    paste0("Normalization method: ", norm_method),
    paste0("Pseudocount: ", pseudocount),
    paste0("Filter: ", lefse_config$filter),
    paste0("Kruskal-Wallis threshold: ", lefse_config$kruskal_threshold),
    paste0("Wilcoxon threshold: ", lefse_config$wilcox_threshold),
    paste0("LDA threshold: ", lefse_config$lda_threshold),
    paste0("P-adjust method: ", lefse_config$p_adjust_method),
    paste0("Seed: ", lefse_config$seed),
    paste0("Sample count: ", if (is.null(physeq)) "not available" else nsamples(physeq)),
    paste0("Retained taxon count: ", if (is.null(physeq)) "not available" else ntaxa(physeq)),
    paste0("Relative abundance closure: ", if (is.null(physeq)) "not available" else validate_relative_abundance_closure(physeq)),
    paste0("lefser version: ", version_or_missing("lefser")),
    paste0("SummarizedExperiment version: ", version_or_missing("SummarizedExperiment"))
  )

  if (!is.null(result_file)) {
    lines <- c(lines, paste0("Result file: ", result_file))
  }
  if (!is.null(plot_file)) {
    lines <- c(lines, paste0("Plot file: ", plot_file))
  }
  if (!is.null(skip_reason)) {
    lines <- c(lines, paste0("Skipped reason: ", skip_reason))
  }
  lines
}

run_lefse_comparison <- function(comp_physeq, spec) {
  message("Running LEfSe comparison: ", spec$name)

  prepared <- prepare_lefse_physeq(
    physeq = comp_physeq,
    spec = spec,
    lefse_config = lefse_config,
    project_config = project_config
  )
  if (!validate_relative_abundance_closure(prepared$physeq)) {
    stop("Prepared LEfSe relative abundance matrix does not sum to 1 for every sample.", call. = FALSE)
  }

  se <- physeq_to_lefse_se(prepared$physeq)
  subclass_col <- spec$subclass_col %||% lefse_config$subclass_col
  if (is_missing_lefse_value(subclass_col)) {
    subclass_col <- NULL
  }

  set.seed(as.integer(lefse_config$seed))
  lefse_res <- lefser::lefser(
    relab = se,
    kruskal.threshold = as.numeric(lefse_config$kruskal_threshold),
    wilcox.threshold = as.numeric(lefse_config$wilcox_threshold),
    lda.threshold = as.numeric(lefse_config$lda_threshold),
    classCol = prepared$class_col,
    subclassCol = subclass_col,
    assay = "relative_abundance",
    method = as.character(lefse_config$p_adjust_method)
  )

  results_df <- format_lefse_results(as.data.frame(lefse_res), prepared$class_levels)
  result_file <- write_lefse_results(
    results_df,
    out_dir = out_dir,
    trial_id = trial_id,
    comparison_name = spec$name
  )
  message("Wrote LEfSe results: ", result_file)

  plot_file <- file.path(
    out_dir,
    "LEfSe",
    spec$name,
    paste0(trial_id, "_", spec$name, "_LEfSeLDAPlot.png")
  )
  lda_plot <- if (nrow(results_df) > 0) {
    lefser::lefserPlot(
      lefse_res,
      title = spec$plot_title %||% spec$name,
      trim.names = FALSE
    )
  } else {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = spec$plot_title %||% spec$name,
        subtitle = "No LEfSe biomarkers met the configured thresholds"
      )
  }
  ggplot2::ggsave(
    filename = plot_file,
    plot = lda_plot,
    width = 10,
    height = max(4, min(16, 0.35 * max(1, nrow(results_df)) + 2)),
    units = "in",
    dpi = 300
  )
  message("Wrote LEfSe plot: ", plot_file)

  note_file <- write_lefse_run_note(
    out_dir = out_dir,
    trial_id = trial_id,
    comparison_name = spec$name,
    note_lines = make_run_note(
      spec,
      status = "completed",
      prepared = prepared,
      result_file = result_file,
      plot_file = plot_file
    )
  )
  message("Wrote LEfSe run note: ", note_file)

  invisible(list(result_file = result_file, plot_file = plot_file, note_file = note_file))
}

CompPhyseq <- load_input_physeq()
comparison_specs <- lefse_comparison_specs(lefse_config, da_config)

for (spec in comparison_specs) {
  tryCatch(
    run_lefse_comparison(CompPhyseq, spec),
    error = function(e) {
      warning(
        "Skipping LEfSe comparison '",
        spec$name,
        "': ",
        conditionMessage(e),
        call. = FALSE
      )
      note_file <- write_lefse_run_note(
        out_dir = out_dir,
        trial_id = trial_id,
        comparison_name = spec$name,
        note_lines = make_run_note(
          spec,
          status = "skipped",
          skip_reason = conditionMessage(e)
        )
      )
      message("Wrote skipped LEfSe run note: ", note_file)
    }
  )
}
