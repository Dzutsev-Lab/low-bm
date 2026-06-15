source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "MetadataSchema.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

config_value <- function(config, name) {
  if (is.null(config)) {
    return(NULL)
  }
  config[[name, exact = TRUE]]
}

analysis_config_value <- function(project_config,
                                  section_config,
                                  name,
                                  default = NULL) {
  config_value(section_config, name) %||% config_value(project_config, name) %||% default
}

apply_project_config_defaults <- function(section_config,
                                          project_config,
                                          keys) {
  section_config <- section_config %||% list()
  for (key in keys) {
    if (is.null(config_value(section_config, key)) && !is.null(config_value(project_config, key))) {
      section_config[[key]] <- config_value(project_config, key)
    }
  }
  section_config
}

analysis_output_dir <- function(project_config,
                                section_config = list(),
                                section_keys = c("output_dir", "io_dir", "out_dir"),
                                default = NULL) {
  for (key in section_keys) {
    value <- config_value(section_config, key)
    if (!is.null(value)) {
      return(value)
    }
  }
  config_value(project_config, "output_dir") %||% default
}

resolve_output_path <- function(path, base_dir = "Exp_Output") {
  if (is.null(path)) {
    return(NULL)
  }
  path <- as.character(path[[1]])
  if (grepl("^/", path)) {
    return(path)
  }
  base_dir <- sub("/+$", "", base_dir)
  if (identical(path, base_dir) || startsWith(path, paste0(base_dir, .Platform$file.sep))) {
    return(path)
  }
  file.path(base_dir, path)
}

truthy_flag <- function(value, default = TRUE) {
  if (is.null(value) || is.na(value) || trimws(as.character(value)) == "") {
    return(default)
  }
  !tolower(trimws(as.character(value))) %in% c("0", "false", "f", "no", "n")
}

trial_to_rdata <- function(trial, base_dir = "Exp_Output") {
  trial_dir <- file.path(base_dir, trial)
  trial_id <- sub("_.*$", "", trial)
  file.path(trial_dir, paste0(trial_id, "_physeq.RData"))
}

batch_row_to_trial_name <- function(batch_row) {
  paste0(batch_row[["trialID"]], "_", batch_row[["trial_descript"]])
}

batch_row_to_physeq_path <- function(batch_row, base_dir = "Exp_Output") {
  if ("physeq_path" %in% names(batch_row) &&
      !is.na(batch_row[["physeq_path"]]) &&
      nzchar(trimws(as.character(batch_row[["physeq_path"]])))) {
    return(as.character(batch_row[["physeq_path"]]))
  }

  trial_name <- batch_row_to_trial_name(batch_row)
  trial_to_rdata(trial_name, base_dir)
}

read_batch_table <- function(path,
                             include_column = NULL,
                             require_canonical = TRUE) {
  batch_df <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    check.names = FALSE
  )
  canonical <- c(
    "trialID",
    "trial_descript",
    "exp_dir",
    "metadata",
    "batch_label",
    "include_processing",
    "include_analysis"
  )

  if (require_canonical) {
    fail_missing_columns(names(batch_df), canonical, paste0("batch table ", path))
  }

  if (!is.null(include_column) && include_column %in% names(batch_df)) {
    keep <- vapply(batch_df[[include_column]], truthy_flag, logical(1), default = TRUE)
    batch_df <- batch_df[keep, , drop = FALSE]
  }

  if (nrow(batch_df) == 0) {
    stop("No batch rows selected from: ", path, call. = FALSE)
  }

  batch_df
}

resolve_batch_physeqs <- function(batch_table,
                                  base_dir = "Exp_Output",
                                  include_column = "include_analysis") {
  batch_df <- read_batch_table(batch_table, include_column = include_column)
  paths <- vapply(
    seq_len(nrow(batch_df)),
    function(i) batch_row_to_physeq_path(batch_df[i, , drop = FALSE], base_dir),
    character(1)
  )
  names(paths) <- batch_df$batch_label
  paths
}

load_physeq <- function(path, object_name = "physeq") {
  if (!file.exists(path)) {
    stop("Missing phyloseq RData file: ", path, call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  load(path, envir = env)

  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop(
      "File does not contain an object named '",
      object_name,
      "': ",
      path,
      call. = FALSE
    )
  }

  physeq <- get(object_name, envir = env)
  if (!inherits(physeq, "phyloseq")) {
    stop("Object '", object_name, "' is not a phyloseq object: ", path, call. = FALSE)
  }

  validate_physeq_metadata(physeq)
}

load_physeqs <- function(paths) {
  lapply(paths, load_physeq)
}

merge_physeqs <- function(physeq_list) {
  if (length(physeq_list) == 0) {
    stop("No phyloseq objects supplied for merging.", call. = FALSE)
  }
  if (length(physeq_list) == 1) {
    return(physeq_list[[1]])
  }
  Reduce(function(x, y) phyloseq::merge_phyloseq(x, y), physeq_list)
}

save_physeq <- function(physeq, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  save(physeq, file = out_file)
  out_file
}

load_yaml_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The R package 'yaml' is required to read config file: ", path, call. = FALSE)
  }
  yaml::read_yaml(path)
}
