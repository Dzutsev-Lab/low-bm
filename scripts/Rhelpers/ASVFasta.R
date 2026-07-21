source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))

is_present_value <- function(x) {
  !is.null(x) &&
    length(x) > 0 &&
    !all(is.na(x)) &&
    any(nzchar(trimws(as.character(x))))
}

first_present_value <- function(x) {
  if (!is_present_value(x)) {
    return(NULL)
  }
  values <- as.character(unlist(x, use.names = FALSE))
  values <- values[!is.na(values) & nzchar(trimws(values))]
  values[[1]]
}

trial_to_asv_fasta <- function(trial, base_dir = "Exp_Output") {
  trial_dir <- file.path(base_dir, trial)
  trial_id <- sub("_.*$", "", trial)
  file.path(trial_dir, paste0(trial_id, "_ASV.fasta"))
}

batch_row_to_asv_fasta_path <- function(batch_row, base_dir = "Exp_Output") {
  asv_fasta_path <- first_present_value(batch_row[["asv_fasta_path"]])
  if (!is.null(asv_fasta_path)) {
    return(asv_fasta_path)
  }

  trial_to_asv_fasta(batch_row_to_trial_name(batch_row), base_dir)
}

resolve_batch_asv_fastas <- function(batch_table,
                                     base_dir = "Exp_Output",
                                     include_column = NULL) {
  batch_df <- read_batch_table(batch_table, include_column = include_column)
  paths <- vapply(
    seq_len(nrow(batch_df)),
    function(i) batch_row_to_asv_fasta_path(batch_df[i, , drop = FALSE], base_dir),
    character(1)
  )
  names(paths) <- vapply(
    seq_len(nrow(batch_df)),
    function(i) batch_row_to_label(batch_df[i, , drop = FALSE]),
    character(1)
  )
  paths
}

asv_fasta_candidates_from_physeq <- function(physeq_path, trial_id = NULL) {
  physeq_path <- first_present_value(physeq_path)
  if (is.null(physeq_path)) {
    return(character(0))
  }

  physeq_dir <- dirname(physeq_path)
  candidates <- file.path(physeq_dir, "MergedASV.fasta")

  trial_id <- first_present_value(trial_id)
  if (!is.null(trial_id)) {
    candidates <- c(candidates, file.path(physeq_dir, paste0(trial_id, "_ASV.fasta")))
  }

  physeq_file <- basename(physeq_path)
  if (grepl("_physeq\\.RData$", physeq_file)) {
    candidates <- c(
      candidates,
      file.path(physeq_dir, sub("_physeq\\.RData$", "_ASV.fasta", physeq_file))
    )
  }

  unique(candidates)
}

project_asv_fasta_candidate_groups <- function(project_config = list(),
                                               base_dir = "Exp_Output",
                                               include_column = NULL) {
  project_config <- project_config %||% list()
  groups <- list()

  compiled_asv_fasta <- first_present_value(config_value(project_config, "compiled_asv_fasta"))
  if (!is.null(compiled_asv_fasta)) {
    groups$compiled_asv_fasta <- compiled_asv_fasta
  }

  asv_fasta <- first_present_value(config_value(project_config, "asv_fasta"))
  if (!is.null(asv_fasta)) {
    groups$asv_fasta <- asv_fasta
  }

  compiled_physeq <- first_present_value(config_value(project_config, "compiled_physeq"))
  if (!is.null(compiled_physeq)) {
    compiled_physeq_candidates <- asv_fasta_candidates_from_physeq(
      physeq_path = compiled_physeq,
      trial_id = config_value(project_config, "trialID")
    )
    for (i in seq_along(compiled_physeq_candidates)) {
      groups[[paste0("compiled_physeq_directory_", i)]] <- compiled_physeq_candidates[[i]]
    }
  }

  batch_table <- first_present_value(config_value(project_config, "batch_table"))
  if (!is.null(batch_table)) {
    groups$batch_table <- resolve_batch_asv_fastas(
      batch_table = batch_table,
      base_dir = base_dir,
      include_column = include_column
    )
  }

  groups
}

resolve_project_asv_fastas <- function(project_config = list(),
                                       base_dir = "Exp_Output",
                                       include_column = NULL) {
  groups <- project_asv_fasta_candidate_groups(
    project_config = project_config,
    base_dir = base_dir,
    include_column = include_column
  )
  searched_paths <- unique(unname(unlist(groups, use.names = FALSE)))

  for (group_name in names(groups)) {
    paths <- unique(as.character(groups[[group_name]]))
    paths <- paths[!is.na(paths) & nzchar(trimws(paths))]
    if (length(paths) > 0 && all(file.exists(paths))) {
      attr(paths, "source") <- group_name
      attr(paths, "searched_paths") <- searched_paths
      return(paths)
    }
  }

  out <- character(0)
  attr(out, "searched_paths") <- searched_paths
  out
}

format_asv_fasta_resolution_error <- function(asv_fasta_paths) {
  searched_paths <- attr(asv_fasta_paths, "searched_paths") %||% character(0)
  searched_text <- if (length(searched_paths) == 0) {
    "  (no candidate paths could be inferred)"
  } else {
    paste0("  - ", searched_paths, collapse = "\n")
  }

  paste(
    "Missing required BLAST candidate preparation input: ASV FASTA.",
    "Searched paths:",
    searched_text,
    "Provide project.asv_fasta or project.compiled_asv_fasta, or rerun the processing pipeline to create <trialID>_ASV.fasta beside the phyloseq output.",
    sep = "\n"
  )
}

read_merged_asv_fasta <- function(fasta_files) {
  fasta_files <- unique(as.character(fasta_files))
  fasta_files <- fasta_files[!is.na(fasta_files) & nzchar(trimws(fasta_files))]

  if (length(fasta_files) == 0) {
    stop("No ASV FASTA files supplied.", call. = FALSE)
  }

  missing_files <- fasta_files[!file.exists(fasta_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing ASV FASTA file(s): ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  seqs_list <- lapply(fasta_files, Biostrings::readDNAStringSet)
  merged_seqs <- do.call(c, seqs_list)

  if (length(merged_seqs) == 0) {
    stop("No sequences found in ASV FASTA file(s): ", paste(fasta_files, collapse = ", "), call. = FALSE)
  }

  asv_ids <- names(merged_seqs)
  seq_values <- as.character(merged_seqs)
  if (any(is.na(asv_ids) | !nzchar(trimws(asv_ids)))) {
    stop("ASV FASTA files contain sequence(s) without ASV IDs.", call. = FALSE)
  }

  duplicated_ids <- unique(asv_ids[duplicated(asv_ids)])
  if (length(duplicated_ids) > 0) {
    conflicting_ids <- duplicated_ids[vapply(
      duplicated_ids,
      function(id) length(unique(seq_values[asv_ids == id])) > 1,
      logical(1)
    )]

    if (length(conflicting_ids) > 0) {
      stop(
        "Found duplicated ASV IDs with conflicting sequences: ",
        paste(conflicting_ids, collapse = ", "),
        call. = FALSE
      )
    }
  }

  keep <- !duplicated(asv_ids)
  out <- Biostrings::DNAStringSet(seq_values[keep])
  names(out) <- asv_ids[keep]
  out
}

write_merged_asv_fasta <- function(fasta_files, out_file) {
  merged_asvs <- read_merged_asv_fasta(fasta_files)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  Biostrings::writeXStringSet(merged_asvs, filepath = out_file)
  out_file
}
