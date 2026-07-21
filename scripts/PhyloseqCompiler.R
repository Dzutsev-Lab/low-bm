library(argparse)
library(phyloseq)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "ASVFasta.R"))

parser <- ArgumentParser()

parser$add_argument("--trial-list",
                    type = "character",
                    default = NULL,
                    help = "Legacy text file with one trial output directory name per line")
parser$add_argument("--batch-table",
                    type = "character",
                    default = NULL,
                    help = "Canonical batch table")
parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with project.batch_table/base_dir/output_dir settings")
parser$add_argument("--physeqs",
                    type = "character",
                    nargs = "+",
                    default = NULL,
                    help = "Explicit list of phyloseq RData paths to compile")
parser$add_argument("--techrep-avg",
                    action = "store_true",
                    default = FALSE,
                    help = "Average OTU counts between technical replicates after merging")
parser$add_argument("--out",
                    type = "character",
                    default = NULL,
                    help = "Output directory for CompPhyseq.RData")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing trial output folders")

args <- parser$parse_args()

project_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project
}

base_dir <- if (!is.null(project_config$base_dir)) project_config$base_dir else args$base_dir
out_dir <- if (!is.null(args$out)) {
  args$out
} else if (!is.null(project_config$output_dir)) {
  project_config$output_dir
} else {
  stop("Provide --out or project.output_dir in --analysis-config.", call. = FALSE)
}

techrep_avg <- isTRUE(args$techrep_avg) || truthy_flag(project_config$techrep_avg, default = FALSE)
asv_fasta_paths <- character(0)

if (!is.null(args$physeqs)) {
  physeq_paths <- args$physeqs
  inferred_asv_paths <- lapply(
    physeq_paths,
    function(path) {
      candidates <- asv_fasta_candidates_from_physeq(path)
      existing_candidates <- candidates[file.exists(candidates)]
      if (length(existing_candidates) == 0) {
        return(NA_character_)
      }
      existing_candidates[[1]]
    }
  )
  asv_fasta_paths <- unlist(inferred_asv_paths, use.names = FALSE)
  asv_fasta_paths <- asv_fasta_paths[!is.na(asv_fasta_paths)]
} else if (!is.null(args$batch_table) || !is.null(project_config$batch_table)) {
  batch_table <- if (!is.null(args$batch_table)) args$batch_table else project_config$batch_table
  physeq_paths <- resolve_batch_physeqs(
    batch_table = batch_table,
    base_dir = base_dir
  )
  asv_fasta_paths <- resolve_batch_asv_fastas(
    batch_table = batch_table,
    base_dir = base_dir
  )
} else if (!is.null(args$trial_list)) {
  trials <- readLines(args$trial_list)
  trials <- trimws(trials)
  trials <- trials[nzchar(trials)]
  if (length(trials) == 0) {
    stop("No trial names found in --trial-list.", call. = FALSE)
  }
  physeq_paths <- vapply(trials, trial_to_rdata, character(1), base_dir = base_dir)
  asv_fasta_paths <- vapply(trials, trial_to_asv_fasta, character(1), base_dir = base_dir)
} else {
  stop("Provide --physeqs, --batch-table, --analysis-config, or --trial-list.", call. = FALSE)
}

physeq <- merge_physeqs(load_physeqs(physeq_paths))

if (techrep_avg) {
  physeq <- average_by_techrep(physeq)
}

out_file <- file.path(out_dir, "CompPhyseq.RData")
save_physeq(physeq, out_file)

message("Merged phyloseq object saved to: ", out_file)

asv_fasta_paths <- unique(asv_fasta_paths)
if (length(asv_fasta_paths) > 0) {
  missing_asv_fastas <- asv_fasta_paths[!file.exists(asv_fasta_paths)]
  if (length(missing_asv_fastas) > 0) {
    warning(
      "Skipping MergedASV.fasta because ASV FASTA file(s) were not found: ",
      paste(missing_asv_fastas, collapse = ", "),
      call. = FALSE
    )
  } else {
    merged_asv_fasta <- file.path(out_dir, "MergedASV.fasta")
    write_merged_asv_fasta(asv_fasta_paths, merged_asv_fasta)
    message("Merged ASV FASTA saved to: ", merged_asv_fasta)
  }
}
