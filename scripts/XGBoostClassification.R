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

io_dir <- xgboost_config$io_dir
trialID <- xgboost_config$trialID
tax_agg_level <- xgboost_config$tax_agg_level
norm_method <- xgboost_config$norm_method
pseudocount <- xgboost_config$pseudocount
class_factor <- xgboost_config$class_factor
batch_adj_covar <- xgboost_config$batch_adj_covar
batch_adj_formula <- xgboost_config$batch_adj_formula
batch_adj_method <- xgboost_config$batch_adj_method

CompPhyseq <- load_physeq(project_config$compiled_physeq)

dir.create(file.path(io_dir,"XGBoostClassification"), recursive = TRUE, showWarnings = FALSE)

set.seed(32)

# ----------------------------
# 1) Extract metadata and make binary label
# ----------------------------
meta <- as(sample_data(CompPhyseq), "data.frame")
meta$SampleID <- rownames(meta)

if (class_factor == "PatientSample") {
    meta <- meta %>%
        mutate(
            class = case_when(
            TumorType %in% c("HCC", "iCC") ~ 1L, # tumor
            TumorType %in% c("HCC-NT", "iCC-NT") ~ 0L, # nontumor
            TRUE ~ NA_integer_
            )
        ) %>%
        filter(!is.na(class), !is.na(PatientID))
} else if (class_factor == "TumorType") {
    meta <- meta %>%
        mutate(
            class = case_when(
            TumorType == "HCC" ~ 1L, # HCC tumors
            TumorType == "iCC" ~ 0L, # iCC tumors
            TRUE ~ NA_integer_
            )
        ) %>%
        filter(!is.na(class), !is.na(PatientID))
}

# Keep only samples in the filtered metadata
SubPhyseq <- prune_samples(meta$SampleID, CompPhyseq)

# Reorder metadata to match sample order in phyloseq object
meta <- meta[match(sample_names(SubPhyseq), meta$SampleID), ]

# -----------------------------------------
# 2) Normalize and adjust OTU
# -----------------------------------------
GlomPhyseq <- tax_glom_rename(SubPhyseq, tax_agg_level = tax_agg_level)
NormPhyseq <- counts_normalization(
    physeq = GlomPhyseq, 
    norm_method = norm_method, 
    pseudocount = pseudocount
)
# TODO: fix to account for samples dropped due to zero divisor problems in normalization 
AdjPhyseq <- batch_adjustment(
    physeq = NormPhyseq, 
    batch_column = batch_adj_covar, 
    design_formula = batch_adj_formula,
    method = batch_adj_method
)

# ----------------------------------------------
# 3) Extract OTU table as samples x taxa matrix
# ----------------------------------------------
otu <- as(otu_table(AdjPhyseq), "matrix")
if (taxa_are_rows(AdjPhyseq)) {
  otu <- t(otu)
}


# ------------------------------------------
# 4) Patient-level train/test split
#   -  Staring 80:20 training:testing split
# ------------------------------------------
patients <- unique(meta$PatientID)
test_patients <- sample(patients, size = ceiling(0.2 * length(patients)))

train_idx <- !(meta$PatientID %in% test_patients)
test_idx  <-  (meta$PatientID %in% test_patients)

X_train <- otu[train_idx, , drop = FALSE]
X_test  <- otu[test_idx, , drop = FALSE]

y_train <- meta$class[train_idx]
y_test  <- meta$class[test_idx]

patient_train <- meta$PatientID[train_idx]

# Convert to sparse matrices for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
dtest  <- xgb.DMatrix(data = as.matrix(X_test),  label = y_test)


# ---------------------------------------------
# 5) Grouped CV by patient on the training set
# ---------------------------------------------
k <- 5
train_patients <- unique(patient_train)

patient_fold <- sample(rep(1:k, length.out = length(train_patients)))
names(patient_fold) <- train_patients

folds <- lapply(1:k, function(f) {
  which(patient_train %in% names(patient_fold)[patient_fold == f])
})


# ----------------------------
# 6) Train XGBoost
# ----------------------------
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

# Fall back to the evaluation log if needed
if (is.null(best_nrounds)) {
  elog <- cv_fit$evaluation_log
  
  # Going off generic naming that should apply to many different versions of xgboost
  auc_col <- grep("^test.*auc", names(elog), value = TRUE)[1]
  
  best_nrounds <- elog$iter[which.max(elog[[auc_col]])]
}

best_nrounds


final_fit <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds,
  verbose = 1
)

# -----------------------------------
# 7) Evaluate on Testing Patient Set
# -----------------------------------
pred_prob <- predict(final_fit, dtest)
pred_class <- ifelse(pred_prob >= 0.5, 1L, 0L)

roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

accuracy <- mean(pred_class == y_test)
sensitivity <- sum(pred_class == 1L & y_test == 1L) / sum(y_test == 1L)
specificity <- sum(pred_class == 0L & y_test == 0L) / sum(y_test == 0L)

cat("AUC:", round(auc_val, 3), "\n")
cat("Accuracy:", round(accuracy, 3), "\n")
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")

png(file.path(io_dir, "XGBoostClassification", paste0(trialID, "_", class_factor, "_ROCPlot.png")), width = 800, height = 800, res = 120)
if (class_factor == "PatientSample") {
    plot_title <- "Tumor vs Nontumor"
} else if (class_factor == "TumorType") {
    plot_title <- "Tumor Type Comparison"
} else {
    plot_title <- class_factor
}
plot(roc_obj, 
     main = paste0(plot_title," ROC by ", if (!is.null(tax_agg_level)) tax_agg_level else "ASV", " (AUC = ", round(auc_val, 3), ")"),
     col = "blue")
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

# ----------------------------
# 8) Feature importance
# ----------------------------
imp <- xgb.importance(model = final_fit)
print(head(imp, 20))
imp_plot <- xgb.ggplot.importance(imp[1:min(20, nrow(imp)), ]) +
            labs(title = paste0(plot_title, " Importance by ", if (!is.null(tax_agg_level)) tax_agg_level else "ASV"))
ggsave(
    filename = file.path(io_dir, "XGBoostClassification", paste0(trialID, "_", class_factor, "_ImportancePlot.png")),
    plot = imp_plot
)

# ----------------------------
# 9) SHAP summary plot
# ----------------------------
shap_top_n <- min(20L, ncol(X_test))
shap_plot_file <- file.path(
    io_dir,
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
