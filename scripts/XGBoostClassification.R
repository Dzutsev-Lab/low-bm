library(phyloseq)
library(xgboost)
library(Matrix)
library(pROC)
library(ggplot2)
library(Ckmeans.1d.dp)
library(argparse)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "XGBoostClassification.R"))

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project and xgboost_classification settings")

args <- parser$parse_args()

if (is.null(args$analysis_config)) {
    stop("Please provide valid analysis config file", call. = FALSE)
}

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
xgboost_config <- normalize_xgboost_config(cfg$xgboost_classification %||% list(), project_config)

out_dir <- analysis_output_dir(
    project_config,
    xgboost_config,
    section_keys = c("out_dir", "output_dir"),
    default = file.path(project_config$base_dir %||% "Exp_Output", "analysis")
)
trialID <- analysis_config_value(project_config, xgboost_config, "trialID", "analysis")

plot_tax_level <- function(spec) {
    if (!is.null(spec$tax_agg_level) && nzchar(trimws(as.character(spec$tax_agg_level)))) {
        spec$tax_agg_level
    } else {
        "ASV"
    }
}

metric_fraction <- function(numerator, denominator) {
    if (denominator == 0) {
        return(NA_real_)
    }
    numerator / denominator
}

run_xgboost_model <- function(spec, comp_physeq) {
    message("Running XGBoost classification model: ", spec$name)
    set.seed(spec$seed)

    warn_xgboost_target_leakage(spec)
    batch_adj_formula <- if (is_missing_xgboost_value(spec$batch_adj_formula)) {
        message(
            "XGBoost model '",
            spec$name,
            "' batch adjustment will use native package-default model design for ",
            spec$batch_adj_method,
            "."
        )
        NULL
    } else {
        as.character(spec$batch_adj_formula[[1]])
    }

    model_physeq <- apply_sample_filter(comp_physeq, spec$sample_filter)
    model_meta <- prepare_xgboost_model_metadata(model_physeq, spec)
    model_physeq <- prune_samples(model_meta$.sample_name, model_physeq)
    model_meta <- model_meta[match(sample_names(model_physeq), model_meta$.sample_name), , drop = FALSE]

    glom_physeq <- tax_glom_rename(model_physeq, tax_agg_level = spec$tax_agg_level)
    norm_physeq <- counts_normalization(
        physeq = glom_physeq,
        norm_method = spec$norm_method,
        pseudocount = spec$pseudocount
    )
    adj_physeq <- batch_adjustment(
        physeq = norm_physeq,
        batch_column = spec$batch_adj_covar,
        design_formula = batch_adj_formula,
        method = spec$batch_adj_method
    )

    model_meta <- model_meta[model_meta$.sample_name %in% sample_names(adj_physeq), , drop = FALSE]
    sub_physeq <- prune_samples(model_meta$.sample_name, adj_physeq)
    model_meta <- model_meta[match(sample_names(sub_physeq), model_meta$.sample_name), , drop = FALSE]
    if (nrow(model_meta) == 0 || length(unique(model_meta$.xgb_class)) < 2) {
        stop(
            "XGBoost model '",
            spec$name,
            "' does not leave at least two classes after normalization and batch adjustment.",
            call. = FALSE
        )
    }

    otu <- otu_samples_by_taxa(sub_physeq)
    split <- make_xgboost_split(model_meta, spec)
    X_train <- otu[split$train_idx, , drop = FALSE]
    X_test <- otu[split$test_idx, , drop = FALSE]
    y_train <- model_meta$.xgb_class[split$train_idx]
    y_test <- model_meta$.xgb_class[split$test_idx]

    dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
    folds <- make_xgboost_cv_folds(model_meta, split$train_idx, spec)

    cv_args <- list(
        params = spec$xgb_params,
        data = dtrain,
        nrounds = spec$nrounds,
        folds = folds,
        verbose = 1,
        maximize = TRUE
    )
    if (spec$early_stopping_rounds > 0) {
        cv_args$early_stopping_rounds <- spec$early_stopping_rounds
    }
    cv_fit <- do.call(xgb.cv, cv_args)

    best_nrounds <- cv_fit$best_iteration
    if (is.null(best_nrounds)) {
        elog <- cv_fit$evaluation_log
        auc_col <- grep("^test.*auc", names(elog), value = TRUE)[1]
        if (is.na(auc_col)) {
            stop("Could not identify XGBoost CV AUC column for model '", spec$name, "'.", call. = FALSE)
        }
        best_nrounds <- elog$iter[which.max(elog[[auc_col]])]
    }

    final_fit <- xgb.train(
        params = spec$xgb_params,
        data = dtrain,
        nrounds = best_nrounds,
        verbose = 1
    )

    pred_prob <- predict(final_fit, dtest)
    pred_class <- ifelse(pred_prob >= spec$prediction_threshold, 1L, 0L)

    roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))

    accuracy <- mean(pred_class == y_test)
    sensitivity <- metric_fraction(sum(pred_class == 1L & y_test == 1L), sum(y_test == 1L))
    specificity <- metric_fraction(sum(pred_class == 0L & y_test == 0L), sum(y_test == 0L))

    cat("Model:", spec$name, "\n")
    cat("AUC:", round(auc_val, 3), "\n")
    cat("Accuracy:", round(accuracy, 3), "\n")
    cat("Sensitivity:", round(sensitivity, 3), "\n")
    cat("Specificity:", round(specificity, 3), "\n")

    output_dir <- file.path(out_dir, "XGBoostClassification")
    tax_level <- plot_tax_level(spec)
    roc_file <- file.path(output_dir, paste0(trialID, "_", spec$safe_name, "_ROCPlot.png"))
    png(roc_file, width = 800, height = 800, res = 120)
    plot(
        roc_obj,
        main = paste0(spec$plot_title, " ROC by ", tax_level, " (AUC = ", round(auc_val, 3), ")"),
        col = "blue"
    )
    legend(
        "bottomright",
        legend = c(
            paste0("AUC = ", round(auc_val, 3)),
            paste0("Accuracy = ", round(accuracy, 3)),
            paste0("Sensitivity = ", round(sensitivity, 3)),
            paste0("Specificity = ", round(specificity, 3))
        ),
        bty = "n"
    )
    dev.off()

    imp <- xgb.importance(model = final_fit)
    print(head(imp, 20))
    if (nrow(imp) > 0) {
        imp_plot <- xgb.ggplot.importance(imp[1:min(20, nrow(imp)), ]) +
            labs(title = paste0(spec$plot_title, " Importance by ", tax_level))
        ggsave(
            filename = file.path(output_dir, paste0(trialID, "_", spec$safe_name, "_ImportancePlot.png")),
            plot = imp_plot
        )
    } else {
        warning("No XGBoost feature importance values returned for model '", spec$name, "'.", call. = FALSE)
    }

    shap_top_n <- min(20L, ncol(X_test))
    shap_plot_file <- file.path(output_dir, paste0(trialID, "_", spec$safe_name, "_SHAPSummaryPlot.png"))

    png(shap_plot_file, width = 1100, height = 900, res = 140)
    shap_plot <- xgb.plot.shap.summary(
        data = as.matrix(X_test),
        model = final_fit,
        top_n = shap_top_n
    )
    if (inherits(shap_plot, "ggplot")) {
        print(shap_plot + labs(title = paste0(spec$plot_title, " SHAP summary by ", tax_level)))
    }
    dev.off()

    cat("Saved ROC plot:", roc_file, "\n")
    cat("Saved SHAP summary plot:", shap_plot_file, "\n")
    invisible(TRUE)
}

CompPhyseq <- load_physeq(project_config$compiled_physeq)
dir.create(file.path(out_dir, "XGBoostClassification"), recursive = TRUE, showWarnings = FALSE)

success_count <- 0L
for (spec in xgboost_config$models) {
    completed <- tryCatch(
        run_xgboost_model(spec, CompPhyseq),
        error = function(e) {
            warning(
                "Skipping XGBoost model '",
                spec$name,
                "': ",
                conditionMessage(e),
                call. = FALSE
            )
            FALSE
        }
    )
    if (isTRUE(completed)) {
        success_count <- success_count + 1L
    }
}

if (success_count == 0L) {
    stop("No XGBoost classification models completed successfully.", call. = FALSE)
}
