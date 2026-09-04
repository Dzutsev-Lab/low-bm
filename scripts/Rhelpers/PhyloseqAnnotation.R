annotation_sample_name_levels <- c(".sample_name", "sample_names")
annotation_missing_tokens <- c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "NULL", "null", "None", "none")

annotation_scalar <- function(value, default = NULL) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]]) ||
      !nzchar(trimws(as.character(value[[1]])))) {
    return(default)
  }
  as.character(value[[1]])
}

annotation_vector <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }
  values <- as.character(unlist(value, use.names = FALSE))
  values <- trimws(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) NULL else values
}

normalize_annotation_key_values <- function(values) {
  values <- as.character(values)
  values <- trimws(values)
  values[is.na(values) | values %in% annotation_missing_tokens] <- NA_character_
  values
}

normalize_annotation_overwrite <- function(overwrite = NULL) {
  if (is.null(overwrite) || length(overwrite) == 0 || is.na(overwrite[[1]])) {
    return("error")
  }
  if (is.logical(overwrite)) {
    return(if (isTRUE(overwrite[[1]])) "replace" else "error")
  }
  value <- tolower(trimws(as.character(overwrite[[1]])))
  if (value %in% c("replace", "overwrite", "true", "t", "yes", "y", "1")) {
    return("replace")
  }
  if (value %in% c("error", "false", "f", "no", "n", "0")) {
    return("error")
  }
  stop("meta_annotation.overwrite must be one of: error, replace.", call. = FALSE)
}

resolve_annotation_level_keys <- function(physeq, level) {
  level <- annotation_scalar(level)
  if (is.null(level)) {
    stop("meta_annotation.level must be set.", call. = FALSE)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  if (level %in% annotation_sample_name_levels) {
    return(list(
      keys = normalize_annotation_key_values(phyloseq::sample_names(physeq)),
      metadata = metadata_df,
      level = level,
      sample_name_level = TRUE
    ))
  }

  if (!level %in% names(metadata_df)) {
    available <- paste(names(metadata_df), collapse = ", ")
    reserved <- paste(annotation_sample_name_levels, collapse = ", ")
    stop(
      "Annotation level '",
      level,
      "' was not found in phyloseq sample_data. Available annotation levels: ",
      available,
      ". Reserved sample-name levels: ",
      reserved,
      ".",
      call. = FALSE
    )
  }

  list(
    keys = normalize_annotation_key_values(metadata_df[[level]]),
    metadata = metadata_df,
    level = level,
    sample_name_level = FALSE
  )
}

annotation_default_key <- function(level, annotation_key = NULL) {
  annotation_key <- annotation_scalar(annotation_key)
  if (!is.null(annotation_key)) {
    return(annotation_key)
  }
  level <- annotation_scalar(level)
  if (level %in% annotation_sample_name_levels) {
    return("SampleName")
  }
  level
}

annotation_row_signature <- function(df) {
  if (ncol(df) == 0) {
    return(rep("", nrow(df)))
  }
  apply(df, 1, function(row) {
    row <- as.character(row)
    row[is.na(row)] <- "<NA>"
    paste(row, collapse = "\r")
  })
}

collapse_annotation_rows <- function(annotation_df,
                                     annotation_key,
                                     columns,
                                     context = "annotation table") {
  key_values <- normalize_annotation_key_values(annotation_df[[annotation_key]])
  valid_rows <- !is.na(key_values)
  selected <- annotation_df[valid_rows, columns, drop = FALSE]
  selected$.annotation_key <- key_values[valid_rows]

  if (nrow(selected) == 0) {
    return(selected[, c(".annotation_key", columns), drop = FALSE])
  }

  duplicate_keys <- unique(selected$.annotation_key[duplicated(selected$.annotation_key)])
  for (key in duplicate_keys) {
    rows <- selected$.annotation_key == key
    signatures <- unique(annotation_row_signature(selected[rows, columns, drop = FALSE]))
    if (length(signatures) > 1) {
      stop(
        context,
        " contains duplicate annotation key '",
        key,
        "' with non-identical annotation values.",
        call. = FALSE
      )
    }
  }

  selected <- selected[!duplicated(selected$.annotation_key), , drop = FALSE]
  selected[, c(".annotation_key", columns), drop = FALSE]
}

make_physeq_annotation_audit <- function(physeq_keys,
                                         annotation_keys,
                                         added_columns,
                                         overwritten_columns,
                                         overwrite) {
  physeq_nonmissing <- physeq_keys[!is.na(physeq_keys)]
  annotation_nonmissing <- annotation_keys[!is.na(annotation_keys)]
  physeq_unique <- unique(physeq_nonmissing)
  annotation_unique <- unique(annotation_nonmissing)
  matched <- intersect(physeq_unique, annotation_unique)
  unmatched_physeq <- setdiff(physeq_unique, annotation_unique)
  unmatched_annotation <- setdiff(annotation_unique, physeq_unique)

  rows <- list(data.frame(
    status = "summary",
    key = NA_character_,
    n_samples = length(physeq_keys),
    n_annotation_rows = length(annotation_keys),
    added_columns = paste(added_columns, collapse = ","),
    overwritten_columns = paste(overwritten_columns, collapse = ","),
    overwrite = overwrite,
    stringsAsFactors = FALSE
  ))

  key_rows <- function(status, keys) {
    if (length(keys) == 0) {
      return(NULL)
    }
    do.call(rbind, lapply(keys, function(key) {
      data.frame(
        status = status,
        key = key,
        n_samples = sum(physeq_keys == key, na.rm = TRUE),
        n_annotation_rows = sum(annotation_keys == key, na.rm = TRUE),
        added_columns = paste(added_columns, collapse = ","),
        overwritten_columns = paste(overwritten_columns, collapse = ","),
        overwrite = overwrite,
        stringsAsFactors = FALSE
      )
    }))
  }

  rows[[length(rows) + 1]] <- key_rows("matched_key", matched)
  rows[[length(rows) + 1]] <- key_rows("unmatched_physeq_key", unmatched_physeq)
  rows[[length(rows) + 1]] <- key_rows("unmatched_annotation_key", unmatched_annotation)
  if (any(is.na(physeq_keys))) {
    rows[[length(rows) + 1]] <- data.frame(
      status = "missing_physeq_key",
      key = NA_character_,
      n_samples = sum(is.na(physeq_keys)),
      n_annotation_rows = 0L,
      added_columns = paste(added_columns, collapse = ","),
      overwritten_columns = paste(overwritten_columns, collapse = ","),
      overwrite = overwrite,
      stringsAsFactors = FALSE
    )
  }
  if (any(is.na(annotation_keys))) {
    rows[[length(rows) + 1]] <- data.frame(
      status = "missing_annotation_key",
      key = NA_character_,
      n_samples = 0L,
      n_annotation_rows = sum(is.na(annotation_keys)),
      added_columns = paste(added_columns, collapse = ","),
      overwritten_columns = paste(overwritten_columns, collapse = ","),
      overwrite = overwrite,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

annotate_physeq_sample_data <- function(physeq,
                                        annotation_df,
                                        level,
                                        annotation_key = NULL,
                                        columns = NULL,
                                        prefix = NULL,
                                        overwrite = NULL) {
  annotation_df <- as.data.frame(annotation_df, stringsAsFactors = FALSE, check.names = FALSE)
  resolved <- resolve_annotation_level_keys(physeq, level)
  metadata_df <- resolved$metadata
  physeq_keys <- resolved$keys
  level <- resolved$level
  annotation_key <- annotation_default_key(level, annotation_key)

  if (!annotation_key %in% names(annotation_df)) {
    stop(
      "Annotation table is missing annotation_key column: ",
      annotation_key,
      call. = FALSE
    )
  }

  columns <- annotation_vector(columns)
  if (is.null(columns)) {
    columns <- setdiff(names(annotation_df), annotation_key)
  }
  if (length(columns) == 0) {
    stop("No annotation columns selected to add.", call. = FALSE)
  }
  missing_columns <- setdiff(columns, names(annotation_df))
  if (length(missing_columns) > 0) {
    stop(
      "Annotation table is missing selected annotation column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  prefix <- annotation_scalar(prefix, default = "")
  destination_columns <- paste0(prefix, columns)
  duplicate_destinations <- unique(destination_columns[duplicated(destination_columns)])
  if (length(duplicate_destinations) > 0) {
    stop(
      "Annotation destination column names are not unique after prefixing: ",
      paste(duplicate_destinations, collapse = ", "),
      call. = FALSE
    )
  }

  overwrite <- normalize_annotation_overwrite(overwrite)
  existing_destinations <- intersect(destination_columns, names(metadata_df))
  if (length(existing_destinations) > 0 && identical(overwrite, "error")) {
    stop(
      "Annotation destination column(s) already exist in phyloseq sample_data: ",
      paste(existing_destinations, collapse = ", "),
      ". Set meta_annotation.overwrite to 'replace' to overwrite them.",
      call. = FALSE
    )
  }

  collapsed <- collapse_annotation_rows(annotation_df, annotation_key, columns)
  match_idx <- match(physeq_keys, collapsed$.annotation_key)

  for (i in seq_along(columns)) {
    values <- collapsed[[columns[[i]]]][match_idx]
    metadata_df[[destination_columns[[i]]]] <- values
  }

  phyloseq::sample_data(physeq) <- phyloseq::sample_data(metadata_df)
  audit <- make_physeq_annotation_audit(
    physeq_keys = physeq_keys,
    annotation_keys = normalize_annotation_key_values(annotation_df[[annotation_key]]),
    added_columns = destination_columns,
    overwritten_columns = existing_destinations,
    overwrite = overwrite
  )

  list(
    physeq = physeq,
    audit = audit,
    added_columns = destination_columns,
    overwritten_columns = existing_destinations
  )
}
