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
                    help = "Meta-analysis YAML with project and meta_compile settings")
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
compile_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project %||% list()
  compile_config <- cfg$meta_compile %||% cfg$compile_phyloseq %||% list()
}

present_values <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  values <- as.character(unlist(x, use.names = FALSE))
  values <- values[!is.na(values) & nzchar(trimws(values))]
  if (length(values) == 0) NULL else values
}

path_label <- function(path) {
  parent <- basename(dirname(path))
  if (!is.na(parent) && nzchar(parent) && parent != ".") {
    return(parent)
  }
  sub("\\.RData$", "", basename(path))
}

normalize_physeq_path_labels <- function(paths) {
  paths <- as.character(paths)
  labels <- names(paths)
  if (is.null(labels)) {
    labels <- rep("", length(paths))
  }
  missing_labels <- is.na(labels) | !nzchar(trimws(labels))
  labels[missing_labels] <- vapply(paths[missing_labels], path_label, character(1))
  names(paths) <- labels
  paths
}

load_annotated_physeqs <- function(paths) {
  paths <- normalize_physeq_path_labels(paths)
  lapply(seq_along(paths), function(i) {
    annotate_physeq_source(
      load_physeq(paths[[i]]),
      source_label = names(paths)[[i]],
      source_path = paths[[i]],
      source_index = i
    )
  })
}

base_dir <- compile_config$base_dir %||% project_config$base_dir %||% args$base_dir
out_dir <- if (!is.null(args$out)) {
  args$out
} else if (!is.null(compile_config$output_dir)) {
  compile_config$output_dir
} else if (!is.null(compile_config$out_dir)) {
  compile_config$out_dir
} else if (!is.null(project_config$output_dir)) {
  project_config$output_dir
} else {
  stop("Provide --out, meta_compile.output_dir, or project.output_dir in --analysis-config.", call. = FALSE)
}

output_physeq <- compile_config$output_physeq %||% "CompPhyseq.RData"
output_asv_fasta <- compile_config$output_asv_fasta %||% "MergedASV.fasta"
techrep_avg <- isTRUE(args$techrep_avg) ||
  truthy_flag(compile_config$techrep_avg, default = truthy_flag(project_config$techrep_avg, default = FALSE))
asv_fasta_paths <- character(0)

config_physeqs <- present_values(compile_config$physeqs)
config_batch_table <- compile_config$batch_table %||% project_config$batch_table
config_trial_list <- compile_config$trial_list

if (!is.null(args$physeqs) || !is.null(config_physeqs)) {
  physeq_paths <- args$physeqs %||% config_physeqs
  physeq_paths <- normalize_physeq_path_labels(physeq_paths)
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
} else if (!is.null(args$batch_table) || !is.null(config_batch_table)) {
  batch_table <- if (!is.null(args$batch_table)) args$batch_table else config_batch_table
  physeq_paths <- resolve_batch_physeqs(
    batch_table = batch_table,
    base_dir = base_dir
  )
  asv_fasta_paths <- resolve_batch_asv_fastas(
    batch_table = batch_table,
    base_dir = base_dir
  )
} else if (!is.null(args$trial_list) || !is.null(config_trial_list)) {
  trial_list <- args$trial_list %||% config_trial_list
  trials <- readLines(trial_list)
  trials <- trimws(trials)
  trials <- trials[nzchar(trials)]
  if (length(trials) == 0) {
    stop("No trial names found in --trial-list.", call. = FALSE)
  }
  physeq_paths <- vapply(trials, trial_to_rdata, character(1), base_dir = base_dir)
  names(physeq_paths) <- trials
  asv_fasta_paths <- vapply(trials, trial_to_asv_fasta, character(1), base_dir = base_dir)
} else {
  stop("Provide --physeqs, --batch-table, --analysis-config, or --trial-list.", call. = FALSE)
}

physeq <- merge_physeqs(load_annotated_physeqs(physeq_paths))

if (techrep_avg) {
  physeq <- average_by_techrep(physeq)
}

out_file <- file.path(out_dir, output_physeq)
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
    merged_asv_fasta <- file.path(out_dir, output_asv_fasta)
    write_merged_asv_fasta(asv_fasta_paths, merged_asv_fasta)
    message("Merged ASV FASTA saved to: ", merged_asv_fasta)
  }
}
