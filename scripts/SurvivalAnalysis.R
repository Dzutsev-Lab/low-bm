library(argparse)
library(phyloseq)
library(ggplot2)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is required for survival analysis.", call. = FALSE)
}

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "SurvivalAnalysis.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()
parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project and survival_analysis settings")
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
  stop("Provide --analysis-config with a survival_analysis section.", call. = FALSE)
}

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
survival_config <- normalize_survival_config(cfg$survival_analysis %||% list(), project_config)

trial_id <- args$trialID %||% analysis_config_value(project_config, survival_config, "trialID", "analysis")
base_dir <- project_config$base_dir %||% args$base_dir
out_dir <- args$out %||% analysis_output_dir(project_config, survival_config, default = file.path(base_dir, "analysis"))

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

write_tsv <- function(df, path, columns = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    if (is.null(columns)) {
      columns <- character(0)
    }
    df <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
  } else if (!is.null(columns)) {
    missing <- setdiff(columns, names(df))
    for (column in missing) {
      df[[column]] <- NA
    }
    df <- df[, columns, drop = FALSE]
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

run_survival_analysis <- function(comp_physeq, spec) {
  analysis_name <- as.character(spec$name)
  safe_name <- sanitize_survival_path_component(analysis_name)
  analysis_dir <- file.path(out_dir, "Survival", safe_name)
  km_dir <- file.path(analysis_dir, "KM")
  dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running survival analysis: ", analysis_name)
  analysis_physeq <- apply_sample_filter(comp_physeq, spec$sample_filter)
  metadata_df <- as.data.frame(phyloseq::sample_data(analysis_physeq), stringsAsFactors = FALSE)
  covariates <- as.character(unlist(spec$covariates %||% character(0), use.names = FALSE))
  prepare_survival_metadata(
    metadata_df,
    time_col = survival_config$time_col,
    status_col = survival_config$status_col,
    patient_id_col = survival_config$patient_id_col,
    event_code = survival_config$event_code,
    censor_code = survival_config$censor_code,
    covariates = covariates
  )

  patient_physeq <- collapse_physeq_by_patient(
    analysis_physeq,
    patient_id_col = survival_config$patient_id_col
  )

  feature_bundle <- build_patient_feature_matrix(
    patient_physeq = patient_physeq,
    spec = spec,
    survival_config = survival_config,
    project_config = project_config
  )
  patient_features <- feature_bundle$patient_features
  feature_map <- feature_bundle$feature_map

  cox <- run_survival_cox_models(
    patient_features = patient_features,
    feature_map = feature_map,
    covariates = covariates,
    min_n = spec$min_n,
    min_events = spec$min_events
  )
  missingness <- model_missingness_table(cox$results, cox$skipped)

  result_cols <- c(
    "feature", "feature_family", "label", "hazard_ratio", "conf_low", "conf_high",
    "coef", "se", "p", "fdr", "n_total", "n_used", "events_used",
    "dropped_missing", "ph_p", "global_ph_p"
  )
  skipped_cols <- c(
    "feature", "feature_family", "label", "reason", "n_total", "n_used",
    "events_used", "dropped_missing"
  )
  missingness_cols <- c(
    "feature", "feature_family", "label", "n_total", "n_used",
    "events_used", "dropped_missing"
  )

  result_file <- write_tsv(
    cox$results,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CoxResults.tsv")),
    columns = result_cols
  )
  skipped_file <- write_tsv(
    cox$skipped,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_SkippedFeatures.tsv")),
    columns = skipped_cols
  )
  feature_file <- write_tsv(
    patient_features,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_PatientFeatureMatrix.tsv"))
  )
  missingness_file <- write_tsv(
    missingness,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_ModelMissingness.tsv")),
    columns = missingness_cols
  )
  taxa_filter_file <- write_tsv(
    feature_bundle$taxa_filter_stats,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_TaxaFilterStats.tsv"))
  )

  message("Wrote Cox results: ", result_file)
  message("Wrote skipped feature report: ", skipped_file)
  message("Wrote patient feature matrix: ", feature_file)
  message("Wrote model missingness report: ", missingness_file)
  message("Wrote taxa filter stats: ", taxa_filter_file)

  if (nrow(cox$results) == 0) {
    message("No successful Cox models for Kaplan-Meier plotting.")
    return(invisible(list(results = cox$results, skipped = cox$skipped)))
  }

  top_results <- cox$results[order(cox$results$fdr, cox$results$p), , drop = FALSE]
  top_results <- head(top_results, as.integer(spec$km_top_n))
  for (i in seq_len(nrow(top_results))) {
    row <- top_results[i, , drop = FALSE]
    plot <- tryCatch(
      plot_km_feature(
        patient_features = patient_features,
        feature = row$feature[[1]],
        feature_label = row$label[[1]],
        covariates = covariates,
        cutpoint = spec$km_cutpoint,
        cox_result = row
      ),
      error = function(e) {
        warning("Skipping KM plot for ", row$label[[1]], ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(plot)) {
      next
    }
    km_file <- file.path(
      km_dir,
      paste0(
        trial_id,
        "_",
        safe_name,
        "_KM_",
        sprintf("%02d", i),
        "_",
        sanitize_survival_path_component(row$label[[1]]),
        ".png"
      )
    )
    ggplot2::ggsave(
      filename = km_file,
      plot = plot,
      width = 8,
      height = 6,
      units = "in",
      dpi = 300
    )
    message("Wrote KM plot: ", km_file)
  }

  invisible(list(results = cox$results, skipped = cox$skipped))
}

CompPhyseq <- load_input_physeq()
for (spec in survival_config$analyses) {
  run_survival_analysis(CompPhyseq, spec)
}
