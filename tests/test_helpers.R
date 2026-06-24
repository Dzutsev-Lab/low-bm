library(phyloseq)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))
source(file.path("scripts", "Rhelpers", "LEfSeAnalysis.R"))

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
expect_error(validate_metadata_df(data.frame(SampleName = "A")), "missing required")

tmp_file <- file.path(tempdir(), "test_physeq.RData")
saved_file <- save_physeq(physeq, tmp_file)
stopifnot(saved_file == tmp_file)
loaded <- load_physeq(tmp_file)
stopifnot(inherits(loaded, "phyloseq"))
stopifnot(nsamples(loaded) == nsamples(physeq))

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
