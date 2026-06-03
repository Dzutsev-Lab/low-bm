library(phyloseq)

source(file.path("scripts", "Rhelpers", "MetadataSchema.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "DifferentialAbundance.R"))

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
      20, 0, 5,
      0, 30, 10,
      0, 40, 10
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("SampleA_rep1", "SampleA_rep2", "SampleB", "NegCtl"),
      c("ASV1", "ASV2", "ASV3")
    )
  )

  metadata <- data.frame(
    SampleName = rownames(otu),
    SampleID = c("SampleA", "SampleA", "SampleB", "NegCtl"),
    SampleType = c("Tumor", "Tumor", "Nontumor", "NegativeControl"),
    PatientID = c("P1", "P1", "P2", "Control"),
    ProcessingBatch = c("B1", "B1", "B2", "B2"),
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

fake_ancombc <- list(res = list(lfc = data.frame(taxon = "a", other = 1)))
expect_error(extract_ancombc_table(fake_ancombc, "lfc", "missing_coef"), "missing required")

if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  message("ANCOMBC is available; full model smoke tests can be added for project fixtures.")
} else {
  message("Skipping ANCOMBC model smoke test because ANCOMBC is not installed.")
}

message("helper tests passed")
