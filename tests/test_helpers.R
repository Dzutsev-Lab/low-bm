library(phyloseq)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))
source(file.path("scripts", "Rhelpers", "LEfSeAnalysis.R"))
source(file.path("scripts", "Rhelpers", "SurvivalAnalysis.R"))

expect_error <- function(expr, pattern = NULL) {
  err <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )
  if (is.null(err)) {
    stop("Expected error but expression succeeded.", call. = FALSE)
  }
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(err))) {
    stop(
      "Expected error matching '",
      pattern,
      "' but got: ",
      conditionMessage(err),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_test_physeq <- function() {
  otu <- matrix(
    c(
      10, 0, 5,
      20, 0, 6,
      0, 30, 10,
      0, 40, 12,
      0, 40, 10
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      c("SampleA_rep1", "SampleA_rep2", "SampleB", "SampleC", "NegCtl"),
      c("ASV1", "ASV2", "ASV3")
    )
  )

  metadata <- data.frame(
    SampleName = rownames(otu),
    SampleID = c("SampleA", "SampleA", "SampleB", "SampleC", "NegCtl"),
    SampleType = c("Tumor", "Tumor", "Nontumor", "Nontumor", "NegativeControl"),
    PatientID = c("P1", "P1", "P2", "P3", "Control"),
    ProcessingBatch = c("B1", "B1", "B2", "B2", "B2"),
    Host_mapped_reads = c(100, 200, 0, 300, 50),
    SurvivalStatus = c(1, 1, 0, 1, 0),
    SurvivalDays = c(100, 100, 200, 300, 400),
    Age = c("60", "NA", "65", "70", "80"),
    Gender = c("Female", "Female", "Male", "Female", "NA"),
    stringsAsFactors = FALSE,
    row.names = rownames(otu)
  )

  tax <- matrix(
    c(
      "d__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Order1", "f__Family1", "g__GenusA", "s__one",
      "d__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Order1", "f__Family1", "g__GenusA", "s__two",
      "d__Bacteria", "p__Proteobacteria", "c__Gamma", "o__Order2", "f__Family2", "g__GenusB", "s__three"
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("ASV1", "ASV2", "ASV3"),
      c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    )
  )

  phyloseq(
    otu_table(otu, taxa_are_rows = FALSE),
    sample_data(validate_metadata_df(metadata)),
    tax_table(tax)
  )
}

physeq <- make_test_physeq()

validated <- validate_metadata_df(as(sample_data(physeq), "data.frame"))
stopifnot("ControlStatus" %in% names(validated))
stopifnot(is.na(validated$Age[[2]]))
stopifnot(is.na(validated$Gender[[5]]))
expect_error(validate_metadata_df(data.frame(SampleName = "A")), "missing required")

sample_level_diff <- as(sample_data(physeq), "data.frame")
sample_level_diff$SampleType[[2]] <- "Nontumor"
sample_level_diff$TumorType <- c("HCC", "HCC-NT", "iCC-NT", "iCC-NT", "NegativeControl")
invisible(validate_metadata_df(sample_level_diff))

bad_patient_meta <- as(sample_data(physeq), "data.frame")
bad_patient_meta$Age[[2]] <- "61"
expect_error(validate_metadata_df(bad_patient_meta), "inconsistent patient-level")

tmp_file <- file.path(tempdir(), "test_physeq.RData")
saved_file <- save_physeq(physeq, tmp_file)
stopifnot(saved_file == tmp_file)
loaded <- load_physeq(tmp_file)
stopifnot(inherits(loaded, "phyloseq"))
stopifnot(nsamples(loaded) == nsamples(physeq))

legacy_physeq <- physeq
legacy_meta <- as(sample_data(legacy_physeq), "data.frame")
legacy_meta$Age[[1]] <- "N/A"
sample_data(legacy_physeq) <- sample_data(legacy_meta)
legacy_loaded <- validate_physeq_metadata(legacy_physeq)
legacy_loaded_meta <- as(sample_data(legacy_loaded), "data.frame")
stopifnot(is.na(legacy_loaded_meta$Age[[1]]))

merged <- merge_physeqs(list(physeq, physeq))
stopifnot(inherits(merged, "phyloseq"))

avg <- average_by_techrep(physeq)
stopifnot("SampleA" %in% sample_names(avg))
avg_mat <- otu_samples_by_taxa(avg)
stopifnot(avg_mat["SampleA", "ASV1"] == 15)

glommed <- tax_glom_rename(physeq, "Genus")
stopifnot("g__GenusA" %in% taxa_names(glommed))
stopifnot(ntaxa(glommed) == 2)

rel <- counts_normalization(physeq, "RelAbund")
rel_mat <- otu_samples_by_taxa(rel)
stopifnot(all(abs(rowSums(rel_mat) - 1) < 1e-8))

expect_error(divide_by_sample_factor(physeq, "MissingColumn"), "missing required")
host_norm <- suppressWarnings(divide_by_sample_factor(physeq, "Host_mapped_reads"))
stopifnot(!"SampleB" %in% sample_names(host_norm))
stopifnot(nsamples(host_norm) == nsamples(physeq) - 1)
host_norm_mat <- otu_samples_by_taxa(host_norm)
stopifnot(host_norm_mat["SampleA_rep1", "ASV1"] == 10 / 100 * 1e6)
expect_error(
  divide_by_sample_factor(physeq, "Host_mapped_reads", drop_invalid = FALSE),
  "all values must be positive"
)

lefse_config <- normalize_lefse_config(list(), project_config = list(norm_method = "log2HostMapped", tax_agg_level = "Genus"))
stopifnot(identical(lefse_config$source_comparisons, "differential_abundance"))
stopifnot(identical(lefse_config$abundance_scale, "project_norm_to_relative_abundance"))
stopifnot(identical(lefse_config$p_adjust_method, "BH"))
stopifnot(identical(lefse_config$filter, "nonzero"))
stopifnot(isTRUE(lefse_config$lda_safe_filter))
stopifnot(identical(as.numeric(lefse_config$relative_abundance_scale), 1e6))
stopifnot(identical(resolve_lefse_norm_method(lefse_config, list(norm_method = "log2HostMapped")), "log2HostMapped"))

spec <- legacy_comparison_spec("PatientSample")
prepared <- prepare_da_physeq(physeq, spec, build_legacy_da_config(
  trial_id = "test",
  comparisons = "PatientSample",
  out_dir = tempdir()
))
stopifnot(inherits(prepared, "phyloseq"))

bad_spec <- spec
bad_spec$factor_levels <- list(SampleType = c("Tumor"))
expect_error(prepare_da_physeq(physeq, bad_spec, list(tax_agg_level = "Genus")), "absent")

surv_meta <- prepare_survival_metadata(
  as(sample_data(physeq), "data.frame"),
  time_col = "SurvivalDays",
  status_col = "SurvivalStatus",
  patient_id_col = "PatientID",
  covariates = c("Age", "Gender")
)
stopifnot(all(na.omit(surv_meta$.survival_status) %in% c(0, 1)))
stopifnot(all(na.omit(surv_meta$.survival_time) > 0))
stopifnot(is.na(surv_meta$Age[[2]]))
stopifnot(is.na(surv_meta$Gender[[5]]))

bad_surv_meta <- as(sample_data(physeq), "data.frame")
bad_surv_meta$SurvivalDays[[1]] <- -1
expect_error(
  prepare_survival_metadata(bad_surv_meta, "SurvivalDays", "SurvivalStatus", "PatientID"),
  "positive survival times"
)
bad_surv_meta <- as(sample_data(physeq), "data.frame")
bad_surv_meta$SurvivalStatus[[1]] <- 2
expect_error(
  prepare_survival_metadata(bad_surv_meta, "SurvivalDays", "SurvivalStatus", "PatientID"),
  "event/censor"
)
bad_surv_meta <- as(sample_data(physeq), "data.frame")
bad_surv_meta$SurvivalDays[[2]] <- 101
expect_error(
  validate_metadata_df(bad_surv_meta),
  "inconsistent patient-level"
)

patient_physeq <- collapse_physeq_by_patient(physeq, patient_id_col = "PatientID")
stopifnot(nsamples(patient_physeq) == length(unique(as(sample_data(physeq), "data.frame")$PatientID)))
stopifnot("P1" %in% sample_names(patient_physeq))
patient_mat <- otu_samples_by_taxa(patient_physeq)
stopifnot(patient_mat["P1", "ASV1"] == 15)

survival_config <- normalize_survival_config(
  list(
    analyses = list(list(
      name = "TestSurvival",
      sample_filter = list(SampleType = "*"),
      covariates = c("Age"),
      tax_agg_level = "Genus",
      pcoa_distance = "bray",
      pcoa_axes = 2,
      taxa_min_prevalence = 0,
      taxa_min_mean_relative_abundance = 0.5,
      min_n = 2,
      min_events = 1
    ))
  ),
  project_config = list()
)
survival_spec <- survival_config$analyses[[1]]
feature_bundle <- build_patient_feature_matrix(patient_physeq, survival_spec, survival_config)
stopifnot("alpha_Shannon" %in% names(feature_bundle$patient_features))
stopifnot(any(grepl("^pcoa_bray_axis", names(feature_bundle$patient_features))))
stopifnot(nrow(feature_bundle$taxa_filter_stats) == 2)
retained_taxa <- feature_bundle$taxa_filter_stats$taxon[feature_bundle$taxa_filter_stats$retained]
stopifnot(identical(retained_taxa, "g__GenusA"))

if (requireNamespace("vegan", quietly = TRUE)) {
  bray_pcoa <- build_pcoa_features(patient_physeq, tax_agg_level = "Genus", pcoa_distance = "bray", pcoa_axes = 1)
  euclidean_pcoa <- build_pcoa_features(patient_physeq, tax_agg_level = "Genus", pcoa_distance = "euclidean", pcoa_axes = 1)
  stopifnot(!isTRUE(all.equal(bray_pcoa$matrix[[1]], euclidean_pcoa$matrix[[1]])))
}

cox_test_df <- data.frame(
  PatientID = paste0("P", 1:6),
  .survival_time = c(100, 130, 160, 210, 250, 300),
  .survival_status = c(1, 0, 1, 0, 1, 0),
  FeatureA = c(0.2, 1.4, 0.8, 1.8, 1.1, 0.5),
  Age = c("50", "NA", "61", "55", "72", "68"),
  stringsAsFactors = FALSE
)
cox_without_age <- fit_cox_feature(
  cox_test_df,
  feature = "FeatureA",
  feature_family = "test",
  covariates = character(0),
  min_n = 6,
  min_events = 1
)
stopifnot(!is.null(cox_without_age$result))
stopifnot(cox_without_age$result$n_used == 6)
cox_with_age <- fit_cox_feature(
  cox_test_df,
  feature = "FeatureA",
  feature_family = "test",
  covariates = "Age",
  min_n = 6,
  min_events = 1
)
stopifnot(is.null(cox_with_age$result))
stopifnot(cox_with_age$skip$n_used == 5)
stopifnot(cox_with_age$skip$dropped_missing == 1)

km_df <- km_group_data(cox_test_df, "FeatureA", cutpoint = "median")
stopifnot(identical(levels(km_df$.km_group), c("Low", "High")))
expect_error(km_group_data(transform(cox_test_df, FeatureA = 1), "FeatureA"), "cannot be split")

lefse_prepared <- prepare_lefse_physeq(
  physeq,
  spec,
  modifyList(lefse_config, list(norm_method = "noNorm")),
  project_config = list()
)
stopifnot(inherits(lefse_prepared$physeq, "phyloseq"))
stopifnot(identical(lefse_prepared$class_levels, c("Nontumor", "Tumor")))
stopifnot(validate_relative_abundance_closure(lefse_prepared$physeq))
stopifnot(all(abs(rowSums(otu_samples_by_taxa(lefse_prepared$physeq)) - 1e6) < 1e-2))
stopifnot(isTRUE(lefse_prepared$lda_filter$enabled))
stopifnot(lefse_prepared$lda_filter$retained_taxa == ntaxa(lefse_prepared$physeq))
stopifnot(ntaxa(filter_lefse_taxa(physeq, "nonzero")) == 3)

three_level_spec <- spec
three_level_spec$sample_filter <- list(SampleType = "*")
three_level_spec$factor_levels <- NULL
expect_error(
  prepare_lefse_physeq(physeq, three_level_spec, modifyList(lefse_config, list(norm_method = "noNorm")), project_config = list()),
  "exactly two observed"
)

mapped_spec <- three_level_spec
mapped_spec$class_map <- list(Control = "NegativeControl", Case = c("Tumor", "Nontumor"))
mapped_prepared <- prepare_lefse_physeq(
  physeq,
  mapped_spec,
  modifyList(lefse_config, list(norm_method = "noNorm", lda_safe_filter = FALSE)),
  project_config = list()
)
stopifnot(identical(mapped_prepared$class_levels, c("Control", "Case")))
stopifnot(validate_relative_abundance_closure(mapped_prepared$physeq))

rel_se <- physeq_to_lefse_se(lefse_prepared$physeq)
stopifnot(inherits(rel_se, "SummarizedExperiment"))
stopifnot("relative_abundance" %in% names(SummarizedExperiment::assays(rel_se)))

fake_lefse <- data.frame(features = c("g__B", "g__A"), scores = c(-2.5, 3.1))
formatted_lefse <- format_lefse_results(fake_lefse, c("Nontumor", "Tumor"))
stopifnot(identical(formatted_lefse$taxon[[1]], "g__A"))
stopifnot(identical(formatted_lefse$enriched_class[[1]], "Tumor"))

fake_ancombc <- list(res = list(lfc = data.frame(taxon = "a", other = 1)))
expect_error(extract_ancombc_table(fake_ancombc, "lfc", "missing_coef"), "missing required")

batch_table_file <- file.path(tempdir(), "batch_table.tsv")
writeLines(
  c(
    paste(
      c(
        "trialID",
        "trial_descript",
        "exp_dir",
        "metadata",
        "batch_label",
        "include_processing",
        "include_analysis"
      ),
      collapse = "\t"
    ),
    paste(
      c(
        "051926.1",
        "TIGER062822_PrelimAnalysis",
        "exp",
        "metadata.xlsx",
        "TIGER062822",
        "false",
        "true"
      ),
      collapse = "\t"
    )
  ),
  batch_table_file
)
batch_df <- read_batch_table(batch_table_file)
stopifnot(identical(batch_df$trialID[[1]], "051926.1"))
stopifnot(identical(batch_row_to_trial_name(batch_df[1, , drop = FALSE]), "051926.1_TIGER062822_PrelimAnalysis"))

if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  message("ANCOMBC is available; full model smoke tests can be added for project fixtures.")
} else {
  message("Skipping ANCOMBC model smoke test because ANCOMBC is not installed.")
}

if (requireNamespace("lefser", quietly = TRUE)) {
  set.seed(as.integer(lefse_config$seed))
  smoke_res <- lefser::lefser(
    relab = rel_se,
    classCol = lefse_prepared$class_col,
    kruskal.threshold = lefse_config$kruskal_threshold,
    wilcox.threshold = lefse_config$wilcox_threshold,
    lda.threshold = lefse_config$lda_threshold,
    assay = "relative_abundance",
    checkAbundances = TRUE,
    method = lefse_config$p_adjust_method
  )
  stopifnot(is.data.frame(smoke_res))
} else {
  message("Skipping lefser smoke test because lefser is not installed.")
}

message("helper tests passed")
