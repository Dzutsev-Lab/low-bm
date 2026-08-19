library(phyloseq)
library(Biostrings)
library(data.table)
library(tibble)
library(argparse)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "TaxaSelection.R"))
source(file.path("scripts", "Rhelpers", "ASVFasta.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with blast_confirmation settings")
parser$add_argument("--DA-comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Differential abundance comparisons to pull significant taxa from")
parser$add_argument("--candidate-comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Candidate comparison directories to create; defaults to --DA-comparisons")
parser$add_argument("--candidate-source",
                    type = "character",
                    default = NULL,
                    choices = c("differential_abundance", "top_abundance", "taxa_file"),
                    help = "How to select taxa for BLAST confirmation")
parser$add_argument("--DA-method",
                    type = "character",
                    help = "Differential abundance method used to perform comparisons")
parser$add_argument("--DA-result-test",
                    type = "character",
                    default = NULL,
                    help = "ANCOM-BC2 result family to select taxa from, such as primary, global, pairwise, dunnet, or trend")
parser$add_argument("--DA-result-contrast",
                    type = "character",
                    default = NULL,
                    help = "Optional ANCOM-BC2 contrast label to select taxa from")
parser$add_argument("--top-n-taxa",
                    type = "integer",
                    default = NULL,
                    help = "Number of top-abundance taxa to select when --candidate-source top_abundance")
parser$add_argument("--candidate-norm-method",
                    type = "character",
                    default = NULL,
                    help = "Normalization method used only for top-abundance taxon selection")
parser$add_argument("--candidate-pseudocount",
                    type = "double",
                    default = NULL,
                    help = "Pseudocount used with --candidate-norm-method for top-abundance selection")
parser$add_argument("--candidate-taxa-file",
                    type = "character",
                    default = NULL,
                    help = "Taxa-selection TSV used when --candidate-source taxa_file")
parser$add_argument("--compiled-asv-fasta",
                    type = "character",
                    help = "FASTA file with all ASV sequences using ASV IDs as headers")
parser$add_argument("--io-dir",
                    type = "character",
                    help = "Analysis output directory")
parser$add_argument("--trialID",
                    type = "character")
parser$add_argument("--taxa-level",
                    type = "character",
                    default = "Genus",
                    help = "Taxonomic level used by DA results and BLAST candidate selection")
parser$add_argument("--compiled-physeq",
                    type = "character")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing per-batch phyloseq outputs")

args <- parser$parse_args()

cfg <- list()
project_config <- list()
blast_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project %||% list()
  blast_config <- cfg$blast_confirmation %||% list()
}

is_missing_value <- function(x) {
  is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(!nzchar(trimws(as.character(x))))
}

fail_missing_input <- function(required_input, context) {
  missing_input <- names(required_input)[vapply(required_input, is_missing_value, logical(1))]
  if (length(missing_input) > 0) {
    stop(
      "Missing required ",
      context,
      ": ",
      paste(missing_input, collapse = ", "),
      call. = FALSE
    )
  }
}

sanitize_path_component <- function(x, fallback = "taxon") {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    x <- fallback
  }
  x <- trimws(x)
  x <- gsub("[[:space:]/\\\\:;,*?\"<>|]+", "_", x)
  x <- gsub("[^A-Za-z0-9._+=@-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  if (!nzchar(x)) fallback else x
}

normalize_da_method <- function(method) {
  method <- toupper(as.character(method))
  method[method == "ANCOMBC"] <- "ANCOMBC2"
  method
}

load_input_physeq <- function() {
  if (!is.null(args$compiled_physeq)) {
    return(load_physeq(args$compiled_physeq))
  }
  if (!is.null(project_config$compiled_physeq) && file.exists(project_config$compiled_physeq)) {
    return(load_physeq(project_config$compiled_physeq))
  }
  if (!is.null(project_config$batch_table)) {
    base_dir <- project_config$base_dir %||% args$base_dir
    physeq_paths <- resolve_batch_physeqs(project_config$batch_table, base_dir = base_dir)
    return(merge_physeqs(load_physeqs(physeq_paths)))
  }
  stop(
    "Provide --compiled-physeq, project.compiled_physeq, or project.batch_table in --analysis-config.",
    call. = FALSE
  )
}

da_results_file_candidates <- function(io_dir, DA_method, trialID, comparison) {
  method <- normalize_da_method(DA_method)
  candidates <- c(file.path(
    io_dir,
    method,
    comparison,
    paste0(trialID, "_", comparison, "_", method, "Results.tsv")
  ))
  if (identical(method, "ANCOMBC2")) {
    candidates <- c(candidates, file.path(
      io_dir,
      "ANCOMBC",
      comparison,
      paste0(trialID, "_", comparison, "_ANCOMBCResults.tsv")
    ))
  }
  candidates
}

read_significant_taxa <- function(io_dir,
                                  DA_method,
                                  trialID,
                                  comparisons,
                                  result_test = NULL,
                                  result_contrast = NULL) {
  sig_taxa_by_comparison <- list()

  for (comparison in comparisons) {
    results_candidates <- da_results_file_candidates(io_dir, DA_method, trialID, comparison)
    existing <- results_candidates[file.exists(results_candidates)]

    if (length(existing) == 0) {
      stop(
        "Missing differential abundance results file. Tried: ",
        paste(results_candidates, collapse = "; "),
        call. = FALSE
      )
    }
    results_file <- existing[[1]]

    DA_results_df <- read.delim(results_file, stringsAsFactors = FALSE)
    if (!is_missing_value(result_test) && "test" %in% names(DA_results_df)) {
      DA_results_df <- DA_results_df[as.character(DA_results_df$test) == as.character(result_test), , drop = FALSE]
    }
    if (!is_missing_value(result_contrast) && "contrast" %in% names(DA_results_df)) {
      DA_results_df <- DA_results_df[as.character(DA_results_df$contrast) == as.character(result_contrast), , drop = FALSE]
    }
    missing_columns <- setdiff(c("taxon", "significance"), names(DA_results_df))
    if (length(missing_columns) > 0) {
      stop(
        "Missing required column(s) in ",
        results_file,
        ": ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    if (nrow(DA_results_df) == 0) {
      stop(
        "No rows remained in ",
        results_file,
        " after applying DA result filters.",
        call. = FALSE
      )
    }

    sig_mask <- !is.na(DA_results_df$significance) & DA_results_df$significance == "Sig"
    sig_taxa <- DA_results_df[sig_mask, "taxon"]
    sig_taxa <- unique(as.character(sig_taxa[!is.na(sig_taxa) & nzchar(trimws(sig_taxa))]))
    sig_taxa_by_comparison[[comparison]] <- sig_taxa
  }

  sig_taxa_by_comparison
}

resolve_candidate_comparisons <- function(candidate_comparisons, DA_comparisons) {
  comparisons <- candidate_comparisons %||% DA_comparisons
  comparisons <- as.character(unlist(comparisons, use.names = FALSE))
  comparisons <- comparisons[!is.na(comparisons) & nzchar(trimws(comparisons))]
  unique(comparisons)
}

apply_candidate_sample_filter <- function(ps, sample_filter = NULL) {
  if (is.null(sample_filter) || length(sample_filter) == 0) {
    return(prune_empty_physeq(ps))
  }

  prune_empty_physeq(apply_sample_filter(ps, sample_filter))
}

read_top_abundance_taxa <- function(ps,
                                    comparisons,
                                    taxa_level_col = "Genus",
                                    top_n_taxa = 10,
                                    norm_method = "noNorm",
                                    pseudocount = 1) {
  selection_ps <- tax_glom_rename(ps, taxa_level_col)
  selection_ps <- prune_empty_physeq(selection_ps)
  selection_ps <- counts_normalization(
    selection_ps,
    norm_method = norm_method,
    pseudocount = pseudocount
  )
  selection_ps <- prune_empty_physeq(selection_ps)

  top_taxa <- select_top_taxa_by_abundance(
    selection_ps,
    top_n = top_n_taxa,
    positive_only = TRUE
  )

  if (length(top_taxa) == 0) {
    stop("No positive-abundance taxa were available for BLAST candidate selection.", call. = FALSE)
  }

  message(
    "Selected top ",
    length(top_taxa),
    " ",
    taxa_level_col,
    " candidate(s): ",
    paste(top_taxa, collapse = ", ")
  )

  out <- vector("list", length(comparisons))
  names(out) <- comparisons
  for (comparison in comparisons) {
    out[[comparison]] <- top_taxa
  }
  out
}

export_significant_taxon_asvs <- function(ps,
                                          asv_fasta,
                                          sig_taxa_by_comparison,
                                          out_dir,
                                          taxa_level_col = "Genus") {
  all_asv_fasta <- read_merged_asv_fasta(asv_fasta)

  tax_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  tax_df <- rownames_to_column(tax_df, "ASVid")

  if (!taxa_level_col %in% colnames(tax_df)) {
    stop("Taxonomic level '", taxa_level_col, "' not found in tax_table(ps).", call. = FALSE)
  }

  abund <- otu_totals_by_taxa(ps)
  abund <- abund[tax_df$ASVid]
  if (any(is.na(abund))) {
    stop("Some ASV IDs in the tax table were not found in the OTU table.", call. = FALSE)
  }

  tax_df$TotalCount <- as.numeric(abund)
  all_manifests <- list()

  for (comparison_name in names(sig_taxa_by_comparison)) {
    comp_taxa <- unique(sig_taxa_by_comparison[[comparison_name]])
    comp_taxa <- comp_taxa[!is.na(comp_taxa) & nzchar(trimws(comp_taxa))]
    comp_dir <- file.path(out_dir, "BlastAnalysis", comparison_name)
    dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

    if (length(comp_taxa) == 0) {
      warning("No significant taxa found for comparison '", comparison_name, "'.", call. = FALSE)
      next
    }

    taxon_path_names <- make.unique(
      vapply(comp_taxa, sanitize_path_component, character(1), fallback = "taxon"),
      sep = "_"
    )
    names(taxon_path_names) <- comp_taxa
    comp_manifest_list <- list()

    for (taxon in comp_taxa) {
      taxon_path_name <- taxon_path_names[[taxon]]
      taxon_dir <- file.path(comp_dir, taxon_path_name)
      dir.create(taxon_dir, recursive = TRUE, showWarnings = FALSE)

      taxon_values <- as.character(tax_df[[taxa_level_col]])
      taxon_asvs <- tax_df$ASVid[!is.na(taxon_values) & taxon_values == taxon]

      if (length(taxon_asvs) == 0) {
        warning(
          "No ASVs found for ",
          taxa_level_col,
          " '",
          taxon,
          "' in comparison '",
          comparison_name,
          "'.",
          call. = FALSE
        )
        next
      }

      missing_fasta_asvs <- setdiff(taxon_asvs, names(all_asv_fasta))
      if (length(missing_fasta_asvs) > 0) {
        warning(
          "Skipping ",
          length(missing_fasta_asvs),
          " ASV(s) for ",
          taxa_level_col,
          " '",
          taxon,
          "' because they are missing from the ASV FASTA input(s).",
          call. = FALSE
        )
      }
      taxon_asvs <- intersect(taxon_asvs, names(all_asv_fasta))
      if (length(taxon_asvs) == 0) {
        warning(
          "No FASTA sequences available for ",
          taxa_level_col,
          " '",
          taxon,
          "' in comparison '",
          comparison_name,
          "'.",
          call. = FALSE
        )
        next
      }

      taxon_fasta <- all_asv_fasta[taxon_asvs]
      m <- tax_df[match(taxon_asvs, tax_df$ASVid), , drop = FALSE]
      m$Comparison <- comparison_name
      m$TaxaLevel <- taxa_level_col
      m$Taxon <- taxon
      m$TaxonPathName <- taxon_path_name
      m$Sequence <- as.character(taxon_fasta)
      taxon_total_count <- sum(m$TotalCount, na.rm = TRUE)
      if (taxon_total_count <= 0) {
        warning(
          "Skipping ",
          taxa_level_col,
          " '",
          taxon,
          "' because selected ASVs have no positive total counts.",
          call. = FALSE
        )
        next
      }
      m$RelativeAbundanceWithinTaxon <- m$TotalCount / taxon_total_count
      if (taxa_level_col == "Genus") {
        m$RelativeAbundanceWithinGenus <- m$RelativeAbundanceWithinTaxon
      }

      priority_cols <- c(
        "Comparison",
        "TaxaLevel",
        "Taxon",
        "TaxonPathName",
        taxa_level_col,
        "ASVid",
        "TotalCount",
        "RelativeAbundanceWithinTaxon",
        if (taxa_level_col == "Genus") "RelativeAbundanceWithinGenus",
        "Sequence"
      )
      priority_cols <- unique(priority_cols[priority_cols %in% colnames(m)])
      m <- m[, c(priority_cols, setdiff(colnames(m), priority_cols)), drop = FALSE]

      manifest_file <- file.path(taxon_dir, paste0(taxon_path_name, "_ASV_manifest.tsv"))
      fasta_file <- file.path(taxon_dir, paste0(taxon_path_name, "_ASV.fasta"))

      fwrite(as.data.table(m), file = manifest_file, sep = "\t", quote = FALSE)
      writeXStringSet(taxon_fasta, filepath = fasta_file)

      comp_manifest_list[[taxon_path_name]] <- m
    }

    if (length(comp_manifest_list) > 0) {
      comp_manifest <- rbindlist(lapply(comp_manifest_list, as.data.table), fill = TRUE)
      comp_manifest_file <- file.path(comp_dir, paste0(comparison_name, "_AllTaxa_ASV_manifest.tsv"))
      fwrite(comp_manifest, file = comp_manifest_file, sep = "\t", quote = FALSE)

      if (taxa_level_col == "Genus") {
        legacy_manifest_file <- file.path(comp_dir, paste0(comparison_name, "_AllGenera_ASV_manifest.tsv"))
        fwrite(comp_manifest, file = legacy_manifest_file, sep = "\t", quote = FALSE)
      }

      all_manifests[[comparison_name]] <- comp_manifest
    }
  }

  invisible(all_manifests)
}

if (!is.null(args$analysis_config)) {
  trialID <- analysis_config_value(project_config, blast_config, "trialID")
  io_dir <- analysis_output_dir(project_config, blast_config, section_keys = c("io_dir", "output_dir", "out_dir"))
  candidate_source <- blast_config$candidate_source %||% "differential_abundance"
  candidate_comparisons <- resolve_candidate_comparisons(blast_config$candidate_comparisons, blast_config$DA_comparisons)
  DA_comparisons <- blast_config$DA_comparisons
  DA_method <- blast_config$DA_method
  taxa_level <- analysis_config_value(project_config, blast_config, "taxa_level", analysis_config_value(project_config, blast_config, "tax_agg_level", "Genus"))
  compiled_asv_fasta <- project_config$compiled_asv_fasta
  top_n_taxa <- blast_config$top_n_taxa %||% 10
  candidate_taxa_file <- blast_config$candidate_taxa_file %||% blast_config$taxa_file
  candidate_sample_filter <- blast_config$candidate_sample_filter %||% blast_config$sample_filter
  candidate_norm_method <- blast_config$candidate_norm_method %||% "noNorm"
  candidate_pseudocount <- blast_config$candidate_pseudocount %||% project_config$pseudocount %||% 1
  DA_result_test <- blast_config$DA_result_test %||% blast_config$da_result_test
  DA_result_contrast <- blast_config$DA_result_contrast %||% blast_config$da_result_contrast
  asv_fasta_project_config <- project_config
  asv_fasta_base_dir <- project_config$base_dir %||% args$base_dir
} else {
  trialID <- args$trialID
  io_dir <- args$io_dir
  candidate_source <- args$candidate_source %||% "differential_abundance"
  candidate_comparisons <- resolve_candidate_comparisons(args$candidate_comparisons, args$DA_comparisons)
  DA_comparisons <- args$DA_comparisons
  DA_method <- args$DA_method
  taxa_level <- args$taxa_level
  compiled_asv_fasta <- args$compiled_asv_fasta
  top_n_taxa <- args$top_n_taxa %||% 10
  candidate_taxa_file <- args$candidate_taxa_file
  candidate_sample_filter <- NULL
  candidate_norm_method <- args$candidate_norm_method %||% "noNorm"
  candidate_pseudocount <- args$candidate_pseudocount %||% 1
  DA_result_test <- args$DA_result_test
  DA_result_contrast <- args$DA_result_contrast
  asv_fasta_project_config <- list(
    trialID = trialID,
    compiled_physeq = args$compiled_physeq,
    compiled_asv_fasta = compiled_asv_fasta
  )
  asv_fasta_base_dir <- args$base_dir
}

asv_fasta_paths <- resolve_project_asv_fastas(
  project_config = asv_fasta_project_config,
  base_dir = asv_fasta_base_dir
)

required_input <- list(
  trialID = trialID,
  io_dir = io_dir,
  candidate_comparisons = candidate_comparisons,
  candidate_source = candidate_source,
  taxa_level = taxa_level
)
fail_missing_input(required_input, "BLAST candidate preparation input")

if (candidate_source == "differential_abundance") {
  fail_missing_input(
    list(
      DA_comparisons = DA_comparisons,
      DA_method = DA_method
    ),
    "differential-abundance BLAST candidate preparation input"
  )
}
if (candidate_source == "taxa_file") {
  fail_missing_input(
    list(candidate_taxa_file = candidate_taxa_file),
    "taxa-file BLAST candidate preparation input"
  )
}

if (length(asv_fasta_paths) == 0) {
  stop(format_asv_fasta_resolution_error(asv_fasta_paths), call. = FALSE)
}

CompPhyseq <- load_input_physeq()
CompPhyseq <- apply_candidate_sample_filter(CompPhyseq, candidate_sample_filter)

if (candidate_source == "differential_abundance") {
  sig_taxa_by_comparison <- read_significant_taxa(
    io_dir = io_dir,
    DA_method = DA_method,
    trialID = trialID,
    comparisons = DA_comparisons,
    result_test = DA_result_test,
    result_contrast = DA_result_contrast
  )
} else if (candidate_source == "top_abundance") {
  sig_taxa_by_comparison <- read_top_abundance_taxa(
    ps = CompPhyseq,
    comparisons = candidate_comparisons,
    taxa_level_col = taxa_level,
    top_n_taxa = top_n_taxa,
    norm_method = candidate_norm_method,
    pseudocount = candidate_pseudocount
  )
} else if (candidate_source == "taxa_file") {
  sig_taxa_by_comparison <- read_taxa_selection_table(
    path = candidate_taxa_file,
    comparisons = candidate_comparisons,
    taxa_level = taxa_level
  )
} else {
  stop("Unsupported BLAST candidate_source: ", candidate_source, call. = FALSE)
}

selection_name <- blast_config$candidate_selection_name %||%
  blast_config$selection_name %||%
  paste0(candidate_source, "_", taxa_level)
selection_name <- sanitize_path_component(selection_name, fallback = "selected_taxa")
selection_file <- file.path(io_dir, "SelectedTaxa", paste0(selection_name, "_taxa.tsv"))
selection_reason <- if (candidate_source == "differential_abundance") {
  paste0(
    "significance == Sig; method = ",
    normalize_da_method(DA_method),
    if (!is_missing_value(DA_result_test)) paste0("; test = ", DA_result_test) else "",
    if (!is_missing_value(DA_result_contrast)) paste0("; contrast = ", DA_result_contrast) else ""
  )
} else if (candidate_source == "top_abundance") {
  paste0("top_n_taxa = ", top_n_taxa, "; norm_method = ", candidate_norm_method)
} else {
  paste0("taxa_file = ", candidate_taxa_file)
}
write_taxa_selection_table(
  selection = sig_taxa_by_comparison,
  out_file = selection_file,
  source = candidate_source,
  taxa_level = taxa_level,
  selection_reason = selection_reason
)
message("Wrote taxa-selection table: ", selection_file)

export_significant_taxon_asvs(
  ps = CompPhyseq,
  asv_fasta = asv_fasta_paths,
  sig_taxa_by_comparison = sig_taxa_by_comparison,
  out_dir = io_dir,
  taxa_level_col = taxa_level
)
