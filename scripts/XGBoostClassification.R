library(phyloseq)
library(dplyr)
library(xgboost)
library(Matrix)
library(pROC)
library(ggplot2)
library(Ckmeans.1d.dp)
library(argparse)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project and ordination settings")

args <- parser$parse_args()

project_config <- list()
xgboost_config <- list()
if (!is.null(args$analysis_config)) {
    cfg <- load_yaml_config(args$analysis_config)
    project_config <- cfg$project %||% list()
    xgboost_config <- cfg$xgboost_classification %||% list()
} else {
    stop("Please provide valid analysis config file")
}

config_value <- function(config, name) {
    config[[name, exact = TRUE]]
}

out_dir <- config_value(xgboost_config, "out_dir")
trialID <- config_value(xgboost_config, "trialID")
tax_agg_level <- config_value(xgboost_config, "tax_agg_level")
norm_method <- config_value(xgboost_config, "norm_method")
pseudocount <- config_value(xgboost_config, "pseudocount")
class_factors_config <- config_value(xgboost_config, "class_factors") %||% config_value(xgboost_config, "class_factor")
batch_adj_covar <- config_value(xgboost_config, "batch_adj_covar")
batch_adj_formula <- config_value(xgboost_config, "batch_adj_formula")
batch_adj_method <- config_value(xgboost_config, "batch_adj_method")

is_missing_config_value <- function(x) {
    is.null(x) ||
        length(x) == 0 ||
        all(is.na(x)) ||
        all(!nzchar(trimws(as.character(x))))
}

target_metadata_vars <- function(class_factor) {
    target_vars <- c("class", "class_factor", class_factor)
    if (identical(class_factor, "PatientSample")) {
        target_vars <- c(target_vars, "TumorType", "SampleType")
    }
    unique(target_vars[!is.na(target_vars) & nzchar(target_vars)])
}

if (is_missing_config_value(class_factors_config)) {
    stop("Provide xgboost_classification.class_factors or xgboost_classification.class_factor.", call. = FALSE)
}
class_factors <- unique(as.character(unlist(class_factors_config, use.names = FALSE)))
class_factors <- class_factors[nzchar(trimws(class_factors))]
if (length(class_factors) == 0) {
    stop("No usable XGBoost class factors were provided.", call. = FALSE)
}

if (is_missing_config_value(batch_adj_method)) {
    batch_adj_method <- "removeBatchEffect"
}

if (is_missing_config_value(batch_adj_formula)) {
    batch_adj_formula <- NULL
    message(
        "XGBoost batch adjustment will use native package-default model design for ",
        batch_adj_method,
        "."
    )
} else {
    batch_adj_formula <- as.character(batch_adj_formula[[1]])
    formula_vars <- tryCatch(
        all.vars(stats::as.formula(batch_adj_formula)),
        error = function(e) {
            warning(
                "Could not parse xgboost_classification.batch_adj_formula for target-leakage screening: ",
                conditionMessage(e),
                call. = FALSE
            )
            character()
        }
    )
    target_vars <- unique(unlist(lapply(class_factors, target_metadata_vars), use.names = FALSE))
    risky_vars <- intersect(formula_vars, target_vars)
    if (length(risky_vars) > 0) {
        warning(
            "xgboost_classification.batch_adj_formula references target-related variable(s): ",
            paste(risky_vars, collapse = ", "),
            ". This can inflate classifier performance because labels influence preprocessing.",
            call. = FALSE
        )
    }
    message("XGBoost batch adjustment formula: ", batch_adj_formula)
}

plot_title_for_class <- function(class_factor) {
    if (class_factor == "PatientSample") {
        "Tumor vs Nontumor"
    } else if (class_factor == "TumorType") {
        "Tumor Type Comparison"
    } else {
        class_factor
    }
}

make_class_metadata <- function(physeq, class_factor) {
    meta <- as(sample_data(physeq), "data.frame")
    meta$SampleID <- rownames(meta)

    if (class_factor == "PatientSample") {
        meta <- meta %>%
            mutate(
                class = case_when(
                    TumorType %in% c("HCC", "iCC") ~ 1L,
                    TumorType %in% c("HCC-NT", "iCC-NT") ~ 0L,
                    TRUE ~ NA_integer_
                )
            ) %>%
            filter(!is.na(class), !is.na(PatientID))
    } else if (class_factor == "TumorType") {
        meta <- meta %>%
            mutate(
                class = case_when(
                    TumorType == "HCC" ~ 1L,
                    TumorType == "iCC" ~ 0L,
                    TRUE ~ NA_integer_
                )
            ) %>%
            filter(!is.na(class), !is.na(PatientID))
    } else if (class_factor %in% names(meta)) {
        values <- as.character(meta[[class_factor]])
        keep <- !is.na(values) & nzchar(trimws(values)) & !is.na(meta$PatientID)
        meta <- meta[keep, , drop = FALSE]
        values <- as.character(meta[[class_factor]])
        levels <- sort(unique(values))
        if (length(levels) != 2) {
            stop(
                "Generic XGBoost class factor '",
                class_factor,
                "' must have exactly two non-missing levels after filtering; found: ",
                paste(levels, collapse = ", "),
                call. = FALSE
            )
        }
        meta$class <- ifelse(values == levels[[2]], 1L, 0L)
        message("Generic class mapping for ", class_factor, ": ", levels[[2]], " = 1, ", levels[[1]], " = 0")
    } else {
        stop("Unsupported XGBoost class factor: ", class_factor, call. = FALSE)
    }

    if (nrow(meta) == 0 || length(unique(meta$class)) < 2) {
        stop("Class factor '", class_factor, "' does not leave at least two classes.", call. = FALSE)
    }

    meta
}

make_patient_split <- function(meta, test_fraction = 0.2, max_attempts = 100L) {
    patients <- unique(meta$PatientID)
    if (length(patients) < 2) {
        stop("Need at least two patients for a patient-level train/test split.", call. = FALSE)
    }

    test_size <- min(max(1L, ceiling(test_fraction * length(patients))), length(patients) - 1L)
    last_split <- NULL
    for (attempt in seq_len(max_attempts)) {
        test_patients <- sample(patients, size = test_size)
        train_idx <- !(meta$PatientID %in% test_patients)
        test_idx <- meta$PatientID %in% test_patients
        last_split <- list(train_idx = train_idx, test_idx = test_idx, test_patients = test_patients)
        if (length(unique(meta$class[train_idx])) == 2 && length(unique(meta$class[test_idx])) == 2) {
            return(last_split)
        }
    }

    stop(
        "Could not create a patient-level train/test split with both classes in train and test after ",
        max_attempts,
        " attempts.",
        call. = FALSE
    )
}

run_xgboost_model <- function(class_factor, adj_physeq) {
    message("Running XGBoost classification for class factor: ", class_factor)
    meta <- make_class_metadata(adj_physeq, class_factor)
    keep_samples <- intersect(meta$SampleID, sample_names(adj_physeq))
    sub_physeq <- prune_samples(keep_samples, adj_physeq)
    meta <- meta[match(sample_names(sub_physeq), meta$SampleID), , drop = FALSE]

    otu <- as(otu_table(sub_physeq), "matrix")
    if (taxa_are_rows(sub_physeq)) {
        otu <- t(otu)
    }

    split <- make_patient_split(meta)
    X_train <- otu[split$train_idx, , drop = FALSE]
    X_test <- otu[split$test_idx, , drop = FALSE]
    y_train <- meta$class[split$train_idx]
    y_test <- meta$class[split$test_idx]
    patient_train <- meta$PatientID[split$train_idx]

    dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

    k <- min(5L, length(unique(patient_train)))
    if (k < 2) {
        stop("Need at least two training patients for grouped cross-validation.", call. = FALSE)
    }
    train_patients <- unique(patient_train)
    patient_fold <- sample(rep(seq_len(k), length.out = length(train_patients)))
    names(patient_fold) <- train_patients
    folds <- lapply(seq_len(k), function(f) {
        which(patient_train %in% names(patient_fold)[patient_fold == f])
    })

    params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = 4,
        eta = 0.05,
        subsample = 0.8,
        colsample_bytree = 0.8,
        min_child_weight = 1,
        gamma = 0
    )

    cv_fit <- xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 500,
        folds = folds,
        early_stopping_rounds = 20,
        verbose = 1,
        maximize = TRUE
    )

    best_nrounds <- cv_fit$best_iteration
    if (is.null(best_nrounds)) {
        elog <- cv_fit$evaluation_log
        auc_col <- grep("^test.*auc", names(elog), value = TRUE)[1]
        best_nrounds <- elog$iter[which.max(elog[[auc_col]])]
    }

    final_fit <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = best_nrounds,
        verbose = 1
    )

    pred_prob <- predict(final_fit, dtest)
    pred_class <- ifelse(pred_prob >= 0.5, 1L, 0L)

    roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))

    accuracy <- mean(pred_class == y_test)
    sensitivity <- sum(pred_class == 1L & y_test == 1L) / sum(y_test == 1L)
    specificity <- sum(pred_class == 0L & y_test == 0L) / sum(y_test == 0L)

    cat("Class factor:", class_factor, "\n")
    cat("AUC:", round(auc_val, 3), "\n")
    cat("Accuracy:", round(accuracy, 3), "\n")
    cat("Sensitivity:", round(sensitivity, 3), "\n")
    cat("Specificity:", round(specificity, 3), "\n")

    plot_title <- plot_title_for_class(class_factor)
    png(file.path(out_dir, "XGBoostClassification", paste0(trialID, "_", class_factor, "_ROCPlot.png")), width = 800, height = 800, res = 120)
    plot(
        roc_obj,
        main = paste0(plot_title, " ROC by ", if (!is.null(tax_agg_level)) tax_agg_level else "ASV", " (AUC = ", round(auc_val, 3), ")"),
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
    imp_plot <- xgb.ggplot.importance(imp[1:min(20, nrow(imp)), ]) +
        labs(title = paste0(plot_title, " Importance by ", if (!is.null(tax_agg_level)) tax_agg_level else "ASV"))
    ggsave(
        filename = file.path(out_dir, "XGBoostClassification", paste0(trialID, "_", class_factor, "_ImportancePlot.png")),
        plot = imp_plot
    )

    shap_top_n <- min(20L, ncol(X_test))
    shap_plot_file <- file.path(
        out_dir,
        "XGBoostClassification",
        paste0(trialID, "_", class_factor, "_SHAPSummaryPlot.png")
    )

    png(shap_plot_file, width = 1100, height = 900, res = 140)
    shap_plot <- xgb.plot.shap.summary(
        data = as.matrix(X_test),
        model = final_fit,
        top_n = shap_top_n
    )
    if (inherits(shap_plot, "ggplot")) {
        print(
            shap_plot +
                labs(title = paste0(plot_title, " SHAP summary by ", if (!is.null(tax_agg_level)) tax_agg_level else "ASV"))
        )
    }
    dev.off()

    cat("Saved SHAP summary plot:", shap_plot_file, "\n")
}

CompPhyseq <- load_physeq(project_config$compiled_physeq)
dir.create(file.path(out_dir, "XGBoostClassification"), recursive = TRUE, showWarnings = FALSE)
set.seed(42)

GlomPhyseq <- tax_glom_rename(CompPhyseq, tax_agg_level = tax_agg_level)
NormPhyseq <- counts_normalization(
    physeq = GlomPhyseq,
    norm_method = norm_method,
    pseudocount = pseudocount
)
AdjPhyseq <- batch_adjustment(
    physeq = NormPhyseq,
    batch_column = batch_adj_covar,
    design_formula = batch_adj_formula,
    method = batch_adj_method
)

for (class_factor in class_factors) {
    tryCatch(
        run_xgboost_model(class_factor, AdjPhyseq),
        error = function(e) {
            warning(
                "Skipping XGBoost class factor '",
                class_factor,
                "': ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )
}
