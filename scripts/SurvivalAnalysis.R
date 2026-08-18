library(argparse)
library(phyloseq)

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
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running survival analysis: ", analysis_name)
  analysis_physeq <- apply_sample_filter(comp_physeq, spec$sample_filter)
  duplicate_policy <- apply_survival_patient_duplicate_policy(
    analysis_physeq,
    spec = spec,
    survival_config = survival_config
  )
  analysis_physeq <- duplicate_policy$physeq
  duplicate_policy_file <- write_tsv(
    duplicate_policy$audit,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_PatientDuplicatePolicy.tsv")),
    columns = names(duplicate_policy$audit)
  )
  message("Wrote patient duplicate policy audit: ", duplicate_policy_file)

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

  feature_bundle <- build_patient_feature_matrix(
    sample_physeq = analysis_physeq,
    spec = spec,
    survival_config = survival_config,
    project_config = project_config
  )
  patient_features <- feature_bundle$patient_features
  feature_map <- feature_bundle$feature_map

  feature_file <- write_tsv(
    patient_features,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_PatientFeatureMatrix.tsv"))
  )
  taxa_filter_file <- write_tsv(
    feature_bundle$taxa_filter_stats,
    file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_TaxaFilterStats.tsv"))
  )

  message("Wrote patient feature matrix: ", feature_file)
  message("Wrote taxa filter stats: ", taxa_filter_file)

  cox <- list(results = data.frame(), skipped = data.frame())
  if (survival_method_enabled(spec, "cox")) {
    cox <- run_survival_cox_models(
      patient_features = patient_features,
      feature_map = feature_map,
      covariates = covariates,
      min_n = spec$min_n,
      min_events = spec$min_events
    )
    missingness <- model_missingness_table(cox$results, cox$skipped)

    result_cols <- c(
      "feature", "feature_family", "label", "sample_strata_col", "sample_stratum",
      "hazard_ratio", "conf_low", "conf_high", "coef", "se", "p", "fdr",
      "n_total", "n_used", "events_used", "dropped_missing", "ph_p", "global_ph_p"
    )
    skipped_cols <- c(
      "feature", "feature_family", "label", "sample_strata_col", "sample_stratum",
      "reason", "n_total", "n_used", "events_used", "dropped_missing"
    )
    missingness_cols <- c(
      "feature", "feature_family", "label", "sample_strata_col", "sample_stratum",
      "n_total", "n_used", "events_used", "dropped_missing"
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
    missingness_file <- write_tsv(
      missingness,
      file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_ModelMissingness.tsv")),
      columns = missingness_cols
    )

    message("Wrote Cox results: ", result_file)
    message("Wrote skipped feature report: ", skipped_file)
    message("Wrote model missingness report: ", missingness_file)
  } else {
    message("Skipping standard Cox per-taxon modeling for analysis: ", analysis_name)
  }

  coda <- NULL
  if (survival_method_enabled(spec, "coda4microbiome")) {
    coda <- run_survival_coda_models(
      sample_physeq = analysis_physeq,
      spec = spec,
      survival_config = survival_config,
      covariates = covariates,
      min_n = spec$min_n,
      min_events = spec$min_events
    )
    coda_signature_cols <- c(
      "analysis", "sample_strata_col", "sample_stratum", "taxon", "label",
      "log_contrast_coefficient", "risk_direction", "abs_coefficient"
    )
    coda_risk_score_cols <- unique(c(
      "PatientID", "sample_strata_col", "sample_stratum",
      ".survival_time", ".survival_status", covariates,
      "risk_score", "risk_group_median"
    ))
    coda_metrics_cols <- c(
      "sample_strata_col", "sample_stratum", "status", "reason",
      "n_total", "n_with_stratum", "n_used", "events_used", "dropped_missing",
      "n_taxa_retained", "n_taxa_selected", "lambda", "alpha", "nfolds",
      "apparent_cindex", "mean_cv_cindex", "sd_cv_cindex", "zero_handling"
    )
    coda_plot_manifest_cols <- c(
      "analysis", "sample_strata_col", "sample_stratum",
      "plot_type", "path", "status", "reason"
    )
    coda_assumption_cols <- c(
      "analysis", "sample_strata_col", "sample_stratum", "diagnostic_model",
      "check", "term", "statistic", "df", "p", "status", "reason",
      "n_used", "events_used"
    )

    coda_signature_file <- write_tsv(
      coda$signature,
      file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaSignature.tsv")),
      columns = coda_signature_cols
    )
    coda_risk_file <- write_tsv(
      coda$risk_scores,
      file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaRiskScores.tsv")),
      columns = coda_risk_score_cols
    )
    coda_metrics_file <- write_tsv(
      coda$metrics,
      file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaModelMetrics.tsv")),
      columns = coda_metrics_cols
    )
    if (isTRUE(spec$coda4microbiome$show_plots)) {
      coda_plot_manifest <- write_coda_plot_outputs(
        coda$plots,
        analysis_dir = analysis_dir,
        trial_id = trial_id,
        safe_name = safe_name,
        options = spec$coda4microbiome
      )
      coda_plot_manifest_file <- write_tsv(
        coda_plot_manifest,
        file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaPlotManifest.tsv")),
        columns = coda_plot_manifest_cols
      )
      message("Wrote coda4microbiome plot manifest: ", coda_plot_manifest_file)
    }
    if (isTRUE(spec$coda4microbiome$assumption_checks)) {
      coda_assumption_file <- write_tsv(
        coda$assumption_checks,
        file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaAssumptionChecks.tsv")),
        columns = coda_assumption_cols
      )
      coda_assumption_plot_manifest <- write_coda_plot_outputs(
        coda$assumption_plots,
        analysis_dir = analysis_dir,
        trial_id = trial_id,
        safe_name = safe_name,
        options = spec$coda4microbiome
      )
      coda_assumption_plot_manifest_file <- write_tsv(
        coda_assumption_plot_manifest,
        file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaAssumptionPlotManifest.tsv")),
        columns = coda_plot_manifest_cols
      )
      message("Wrote coda4microbiome assumption checks: ", coda_assumption_file)
      message("Wrote coda4microbiome assumption plot manifest: ", coda_assumption_plot_manifest_file)
    }
    legacy_coda_filter_file <- file.path(analysis_dir, paste0(trial_id, "_", safe_name, "_CodaTaxaFilterStats.tsv"))
    if (file.exists(legacy_coda_filter_file)) {
      unlink(legacy_coda_filter_file)
      message("Removed deprecated duplicate coda4microbiome taxa filter stats: ", legacy_coda_filter_file)
    }

    message("Wrote coda4microbiome signature: ", coda_signature_file)
    message("Wrote coda4microbiome risk scores: ", coda_risk_file)
    message("Wrote coda4microbiome model metrics: ", coda_metrics_file)
  }

  invisible(list(results = cox$results, skipped = cox$skipped, coda = coda))
}

CompPhyseq <- load_input_physeq()
for (spec in survival_config$analyses) {
  run_survival_analysis(CompPhyseq, spec)
}
