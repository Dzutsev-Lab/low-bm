library(phyloseq)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "ASVFasta.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "TaxaSelection.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))
source(file.path("scripts", "Rhelpers", "LEfSeAnalysis.R"))
source(file.path("scripts", "Rhelpers", "SurvivalAnalysis.R"))
source(file.path("scripts", "Rhelpers", "AbundanceBarPlots.R"))

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

no_patient_meta <- as(sample_data(physeq), "data.frame")
no_patient_meta$PatientID <- NULL
no_patient_validated <- validate_metadata_df(no_patient_meta)
stopifnot(!"PatientID" %in% names(no_patient_validated))

no_patient_physeq <- physeq
sample_data(no_patient_physeq) <- sample_data(no_patient_meta)
no_patient_physeq <- validate_physeq_metadata(no_patient_physeq)
stopifnot(inherits(no_patient_physeq, "phyloseq"))

sample_level_diff <- as(sample_data(physeq), "data.frame")
sample_level_diff$SampleType[[2]] <- "Nontumor"
sample_level_diff$TumorType <- c("HCC", "HCC-NT", "iCC-NT", "iCC-NT", "NegativeControl")
invisible(validate_metadata_df(sample_level_diff))

control_diff <- as(sample_data(physeq), "data.frame")
control_extra <- control_diff["NegCtl", , drop = FALSE]
rownames(control_extra) <- "NegCtl_rep2"
control_extra$SampleName <- "NegCtl_rep2"
control_extra$SampleID <- "NegCtl_rep2"
control_extra$Age <- "90"
control_extra$Gender <- "Male"
control_diff <- rbind(control_diff, control_extra)
control_validated <- validate_metadata_df(control_diff)
stopifnot(nrow(control_validated) == nrow(control_diff))
stopifnot(sum(as.character(control_validated$PatientID) == "Control") == 2)
assert_patient_metadata_consistency(control_diff)

bad_patient_meta <- as(sample_data(physeq), "data.frame")
bad_patient_meta$Age[[2]] <- "61"
expect_error(assert_patient_metadata_consistency(bad_patient_meta), "inconsistent patient-level")
dropped_patient_meta <- suppressWarnings(validate_metadata_df(bad_patient_meta))
stopifnot(!"P1" %in% as.character(dropped_patient_meta$PatientID))
stopifnot(nrow(dropped_patient_meta) == nrow(bad_patient_meta) - 2)

tmp_file <- file.path(tempdir(), "test_physeq.RData")
saved_file <- save_physeq(physeq, tmp_file)
stopifnot(saved_file == tmp_file)
loaded <- load_physeq(tmp_file)
stopifnot(inherits(loaded, "phyloseq"))
stopifnot(nsamples(loaded) == nsamples(physeq))

asv_test_dir <- tempfile("asv_fasta_paths_")
dir.create(asv_test_dir, recursive = TRUE)
write_test_fasta <- function(path, seqs) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fasta <- Biostrings::DNAStringSet(unname(seqs))
  names(fasta) <- names(seqs)
  Biostrings::writeXStringSet(fasta, filepath = path)
}

explicit_compiled_fasta <- file.path(asv_test_dir, "ExplicitMerged.fasta")
explicit_asv_fasta <- file.path(asv_test_dir, "ExplicitSingle.fasta")
trial_dir <- file.path(asv_test_dir, "010126.1_batch1")
canonical_asv_fasta <- file.path(trial_dir, "010126.1_ASV.fasta")
override_asv_fasta <- file.path(asv_test_dir, "override", "batch2.ASV.fasta")
write_test_fasta(explicit_compiled_fasta, c(ASV_1 = "ACGT"))
write_test_fasta(explicit_asv_fasta, c(ASV_2 = "TGCA"))
write_test_fasta(canonical_asv_fasta, c(ASV_3 = "AAAA"))
write_test_fasta(override_asv_fasta, c(ASV_4 = "CCCC"))

stopifnot(identical(
  as.character(resolve_project_asv_fastas(list(compiled_asv_fasta = explicit_compiled_fasta), base_dir = asv_test_dir)),
  explicit_compiled_fasta
))
stopifnot(identical(
  as.character(resolve_project_asv_fastas(list(asv_fasta = explicit_asv_fasta), base_dir = asv_test_dir)),
  explicit_asv_fasta
))
stopifnot(identical(
  as.character(resolve_project_asv_fastas(
    list(trialID = "010126.1", compiled_physeq = file.path(trial_dir, "010126.1_physeq.RData")),
    base_dir = asv_test_dir
  )),
  canonical_asv_fasta
))

asv_batch_table <- file.path(asv_test_dir, "batch_table.tsv")
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
        "include_analysis",
        "asv_fasta_path"
      ),
      collapse = "\t"
    ),
    paste(c("010126.1", "batch1", "exp1", "metadata1.xlsx", "batch1", "true", "true", ""), collapse = "\t"),
    paste(c("010126.2", "batch2", "exp2", "metadata2.xlsx", "batch2", "true", "true", override_asv_fasta), collapse = "\t")
  ),
  asv_batch_table
)
resolved_batch_fastas <- resolve_project_asv_fastas(
  list(batch_table = asv_batch_table),
  base_dir = asv_test_dir
)
stopifnot(identical(as.character(resolved_batch_fastas), c(canonical_asv_fasta, override_asv_fasta)))

dup_fasta <- file.path(asv_test_dir, "dup.ASV.fasta")
bad_dup_fasta <- file.path(asv_test_dir, "bad_dup.ASV.fasta")
write_test_fasta(dup_fasta, c(ASV_1 = "ACGT", ASV_5 = "GGGG"))
write_test_fasta(bad_dup_fasta, c(ASV_1 = "TTTT"))
merged_fasta <- read_merged_asv_fasta(c(explicit_compiled_fasta, dup_fasta))
stopifnot(identical(names(merged_fasta), c("ASV_1", "ASV_5")))
expect_error(
  read_merged_asv_fasta(c(explicit_compiled_fasta, bad_dup_fasta)),
  "conflicting sequences"
)

legacy_physeq <- physeq
legacy_meta <- as(sample_data(legacy_physeq), "data.frame")
legacy_meta$Age[[1]] <- "N/A"
sample_data(legacy_physeq) <- sample_data(legacy_meta)
legacy_loaded <- validate_physeq_metadata(legacy_physeq)
legacy_loaded_meta <- as(sample_data(legacy_loaded), "data.frame")
stopifnot(is.na(legacy_loaded_meta$Age[[1]]))

inconsistent_physeq <- physeq
inconsistent_meta <- as(sample_data(inconsistent_physeq), "data.frame")
inconsistent_meta$Age[[2]] <- "61"
sample_data(inconsistent_physeq) <- sample_data(inconsistent_meta)
pruned_inconsistent <- suppressWarnings(validate_physeq_metadata(inconsistent_physeq))
pruned_meta <- as(sample_data(pruned_inconsistent), "data.frame")
stopifnot(!"P1" %in% as.character(pruned_meta$PatientID))
stopifnot(nsamples(pruned_inconsistent) == nsamples(inconsistent_physeq) - 2)
stopifnot(!any(c("SampleA_rep1", "SampleA_rep2") %in% sample_names(pruned_inconsistent)))

merged <- merge_physeqs(list(physeq, physeq))
stopifnot(inherits(merged, "phyloseq"))

avg <- average_by_techrep(physeq)
stopifnot("SampleA" %in% sample_names(avg))
avg_mat <- otu_samples_by_taxa(avg)
stopifnot(avg_mat["SampleA", "ASV1"] == 15)

glommed <- tax_glom_rename(physeq, "Genus")
stopifnot("g__GenusA" %in% taxa_names(glommed))
stopifnot(ntaxa(glommed) == 2)
ranked_glommed <- rank_taxa_by_abundance(glommed)
stopifnot(identical(ranked_glommed$taxon, c("g__GenusA", "g__GenusB")))
stopifnot(identical(select_top_taxa_by_abundance(glommed, top_n = 1), "g__GenusA"))
expect_error(select_top_taxa_by_abundance(glommed, top_n = 0), "positive integer")

rel <- counts_normalization(physeq, "RelAbund")
rel_mat <- otu_samples_by_taxa(rel)
stopifnot(all(abs(rowSums(rel_mat) - 1) < 1e-8))

barplot_config <- normalize_abundance_barplot_config(
  list(plots = list(list(name = "TestPlot", x = "SampleID"))),
  project_config = list()
)
stopifnot(identical(barplot_config$plots[[1]]$norm_method, "noNorm"))

barplot_spec <- list(
  name = "TestPlot",
  x = "SampleID",
  facet = "SampleType",
  sample_filter = list(SampleType = "*"),
  tax_agg_level = "Genus",
  fill_tax_level = "Genus",
  taxa_display = "top_n",
  top_n = 1,
  norm_method = "noNorm",
  pseudocount = 1
)
validate_barplot_metadata_columns(physeq, barplot_spec)
expect_error(
  validate_barplot_metadata_columns(physeq, modifyList(barplot_spec, list(facet = "MissingFacet"))),
  "missing required"
)

glommed_for_barplot <- tax_glom_rename(physeq, "Genus")
collapsed_barplot <- collapse_top_taxa(glommed_for_barplot, top_n = 1)
stopifnot(ntaxa(collapsed_barplot) == 2)
stopifnot("Other" %in% taxa_names(collapsed_barplot))
stopifnot(all(rowSums(otu_samples_by_taxa(collapsed_barplot)) == rowSums(otu_samples_by_taxa(glommed_for_barplot))))

prepared_top_n <- prepare_abundance_barplot_physeq(physeq, barplot_spec)
stopifnot("Other" %in% taxa_names(prepared_top_n))
prepared_all <- prepare_abundance_barplot_physeq(physeq, modifyList(barplot_spec, list(taxa_display = "all")))
stopifnot(!"Other" %in% taxa_names(prepared_all))
stopifnot(ntaxa(prepared_all) == ntaxa(glommed_for_barplot))

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
dropped_surv_meta <- suppressWarnings(validate_metadata_df(bad_surv_meta))
stopifnot(!"P1" %in% as.character(dropped_surv_meta$PatientID))

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
