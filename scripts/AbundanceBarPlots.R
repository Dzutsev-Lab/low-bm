library(argparse)
library(phyloseq)
library(ggplot2)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "AbundanceBarPlots.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()
parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project and abundance_barplots settings")
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
parser$add_argument("--trialID",
                    type = "character",
                    default = NULL,
                    help = "ID to attach to output files")
parser$add_argument("--out",
                    type = "character",
                    default = NULL,
                    help = "Output directory")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")

args <- parser$parse_args()

if (is.null(args$analysis_config)) {
  stop("Provide --analysis-config with an abundance_barplots section.", call. = FALSE)
}

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
barplot_config <- normalize_abundance_barplot_config(cfg$abundance_barplots %||% list(), project_config)

trial_id <- args$trialID %||% analysis_config_value(project_config, barplot_config, "trialID", "analysis")
base_dir <- project_config$base_dir %||% args$base_dir
out_dir <- args$out %||% analysis_output_dir(project_config, barplot_config, default = file.path(base_dir, "analysis"))

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
    "Provide --compiled-physeq, --physeqs, --batch-table, or project.compiled_physeq/project.batch_table.",
    call. = FALSE
  )
}

write_barplot <- function(plot, spec, plot_dir, safe_name) {
  width <- as.numeric(spec$plot_width %||% 14)
  height <- as.numeric(spec$plot_height %||% 8)
  dpi <- as.numeric(spec$dpi %||% 300)
  png_file <- file.path(plot_dir, paste0(trial_id, "_", safe_name, "_AbundanceBarPlot.png"))
  ggplot2::ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  message("Wrote abundance bar plot: ", png_file)

  if (truthy_flag(spec$write_pdf, default = FALSE)) {
    pdf_file <- file.path(plot_dir, paste0(trial_id, "_", safe_name, "_AbundanceBarPlot.pdf"))
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = "in"
    )
    message("Wrote abundance bar plot PDF: ", pdf_file)
  }
}

run_barplot_spec <- function(comp_physeq, spec) {
  plot_name <- as.character(spec$name)
  safe_name <- sanitize_barplot_path_component(plot_name)
  plot_dir <- file.path(out_dir, "AbundanceBarPlots", safe_name)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running abundance bar plot: ", plot_name)
  plot_physeq <- prepare_abundance_barplot_physeq(comp_physeq, spec)
  plot <- build_abundance_barplot(plot_physeq, spec)
  write_barplot(plot, spec, plot_dir, safe_name)
  invisible(plot)
}

CompPhyseq <- load_input_physeq()
for (spec in barplot_config$plots) {
  run_barplot_spec(CompPhyseq, spec)
}
