required_metadata_columns <- c("SampleName", "SampleID", "SampleType")
required_physeq_metadata_columns <- c("SampleID", "SampleType")
metadata_missing_tokens <- c("", "NA", "N/A", "na", "n/a", "NaN", "nan", "NULL", "null", "None", "none")
patient_consistency_excluded_columns <- c(
  "SampleName",
  "SampleID",
  "SampleType",
  "PatientID",
  "IsControl",
  "ControlStatus",
  "ProcessingBatch",
  "SequencingBatch",
  "Raw_reads",
  "Host_mapped_reads",
  "TumorType",
  "cCluster",
  "Microbiome",
  "Expression (Microarray)",
  "Expression..Microarray.",
  "TumorMetabolomics",
  "Panel WES deep-seq (mutation - 500 genes)",
  "Panel.WES.deep.seq..mutation...500.genes.",
  "CNV SNP",
  "CNV.SNP",
  "Reason to exclude"
)

fail_missing_columns <- function(columns, required, context = "data") {
  missing <- setdiff(required, columns)
  if (length(missing) > 0) {
    stop(
      context,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

is_metadata_missing_like <- function(x) {
  is.na(x) | trimws(as.character(x)) %in% metadata_missing_tokens
}

standardize_metadata_missing <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.character(x)) {
    x <- trimws(x)
    x[x %in% metadata_missing_tokens] <- NA_character_
  }
  x
}

standardize_metadata_missing_df <- function(metadata_df, columns = names(metadata_df)) {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  columns <- intersect(columns, names(metadata_df))
  for (column in columns) {
    metadata_df[[column]] <- standardize_metadata_missing(metadata_df[[column]])
  }
  metadata_df
}

fail_empty_required_values <- function(metadata_df, required, context = "metadata") {
  for (column in required) {
    empty <- is_metadata_missing_like(metadata_df[[column]])
    if (any(empty)) {
      stop(
        context,
        " contains empty value(s) in required column: ",
        column,
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

patient_consistency_columns <- function(metadata_df,
                                        patient_id_col = "PatientID",
                                        exclude = patient_consistency_excluded_columns,
                                        columns = NULL) {
  if (!is.null(columns)) {
    return(intersect(as.character(columns), names(metadata_df)))
  }
  setdiff(names(metadata_df), unique(c(patient_id_col, exclude)))
}

patient_consistency_sample_mask <- function(metadata_df) {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)

  if ("SampleType" %in% names(metadata_df)) {
    sample_type <- standardize_metadata_missing(metadata_df$SampleType)
    return(!is_metadata_missing_like(sample_type) & !grepl("Control", as.character(sample_type)))
  }

  if ("ControlStatus" %in% names(metadata_df)) {
    control_status <- standardize_metadata_missing(metadata_df$ControlStatus)
    return(!is_metadata_missing_like(control_status) & as.character(control_status) == "PatientSample")
  }

  if ("IsControl" %in% names(metadata_df)) {
    is_control <- as.logical(metadata_df$IsControl)
    return(is.na(is_control) | !is_control)
  }

  rep(TRUE, nrow(metadata_df))
}

patient_metadata_inconsistencies <- function(metadata_df,
                                             patient_id_col = "PatientID",
                                             columns = NULL,
                                             exclude = patient_consistency_excluded_columns,
                                             context = "metadata") {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), patient_id_col, context)

  metadata_df[[patient_id_col]] <- standardize_metadata_missing(metadata_df[[patient_id_col]])
  if (any(is_metadata_missing_like(metadata_df[[patient_id_col]]))) {
    stop(context, " contains empty value(s) in patient ID column: ", patient_id_col, call. = FALSE)
  }

  sample_mask <- patient_consistency_sample_mask(metadata_df)
  metadata_df <- metadata_df[sample_mask, , drop = FALSE]
  if (nrow(metadata_df) == 0) {
    return(data.frame(
      PatientID = character(0),
      column = character(0),
      values = character(0),
      stringsAsFactors = FALSE
    ))
  }

  check_columns <- patient_consistency_columns(metadata_df, patient_id_col, exclude, columns)
  if (length(check_columns) == 0) {
    return(data.frame(
      PatientID = character(0),
      column = character(0),
      values = character(0),
      stringsAsFactors = FALSE
    ))
  }

  patient_ids <- as.character(metadata_df[[patient_id_col]])
  duplicated_patients <- unique(patient_ids[duplicated(patient_ids) | duplicated(patient_ids, fromLast = TRUE)])
  if (length(duplicated_patients) == 0) {
    return(data.frame(
      PatientID = character(0),
      column = character(0),
      values = character(0),
      stringsAsFactors = FALSE
    ))
  }

  issues <- list()
  for (column in check_columns) {
    values <- standardize_metadata_missing(metadata_df[[column]])
    for (patient_id in duplicated_patients) {
      patient_values <- values[patient_ids == patient_id]
      patient_values <- patient_values[!is_metadata_missing_like(patient_values)]
      unique_values <- unique(as.character(patient_values))
      if (length(unique_values) > 1) {
        issues[[length(issues) + 1]] <- data.frame(
          PatientID = patient_id,
          column = column,
          values = paste(unique_values, collapse = ", "),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(issues) == 0) {
    return(data.frame(
      PatientID = character(0),
      column = character(0),
      values = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, issues)
}

format_patient_consistency_issues <- function(issue_df, preview_n = 10) {
  if (is.null(issue_df) || nrow(issue_df) == 0) {
    return("")
  }
  preview_n <- min(preview_n, nrow(issue_df))
  preview <- paste(
    paste0(
      issue_df$PatientID[seq_len(preview_n)],
      " / ",
      issue_df$column[seq_len(preview_n)],
      " = [",
      issue_df$values[seq_len(preview_n)],
      "]"
    ),
    collapse = "; "
  )
  more <- if (nrow(issue_df) > preview_n) paste0("; ... and ", nrow(issue_df) - preview_n, " more") else ""
  paste0(preview, more)
}

assert_patient_metadata_consistency <- function(metadata_df,
                                                patient_id_col = "PatientID",
                                                columns = NULL,
                                                exclude = patient_consistency_excluded_columns,
                                                context = "metadata") {
  issue_df <- patient_metadata_inconsistencies(
    metadata_df,
    patient_id_col = patient_id_col,
    columns = columns,
    exclude = exclude,
    context = context
  )
  if (nrow(issue_df) > 0) {
    stop(
      context,
      " contains inconsistent patient-level value(s) within ",
      patient_id_col,
      ": ",
      format_patient_consistency_issues(issue_df),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

drop_inconsistent_patient_metadata <- function(metadata_df,
                                               patient_id_col = "PatientID",
                                               columns = NULL,
                                               exclude = patient_consistency_excluded_columns,
                                               context = "metadata") {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  issue_df <- patient_metadata_inconsistencies(
    metadata_df,
    patient_id_col = patient_id_col,
    columns = columns,
    exclude = exclude,
    context = context
  )
  if (nrow(issue_df) == 0) {
    return(metadata_df)
  }

  patient_ids <- as.character(metadata_df[[patient_id_col]])
  sample_mask <- patient_consistency_sample_mask(metadata_df)
  dropped_patients <- unique(issue_df$PatientID)
  drop_rows <- sample_mask & patient_ids %in% dropped_patients
  dropped_samples <- rownames(metadata_df)[drop_rows]
  warning(
    context,
    " contains inconsistent patient-level value(s) within ",
    patient_id_col,
    ". Dropping ",
    length(dropped_patients),
    " patient(s) and ",
    length(dropped_samples),
    " associated sample(s): ",
    format_patient_consistency_issues(issue_df),
    call. = FALSE
  )

  metadata_df <- metadata_df[!drop_rows, , drop = FALSE]
  if (nrow(metadata_df) == 0) {
    stop(
      context,
      " has no samples remaining after dropping inconsistent patient-level entries.",
      call. = FALSE
    )
  }
  metadata_df
}

derive_control_status <- function(metadata_df) {
  fail_missing_columns(names(metadata_df), "SampleType", "metadata")

  metadata_df$IsControl <- grepl("Control", as.character(metadata_df$SampleType))
  metadata_df$ControlStatus <- ifelse(
    metadata_df$IsControl,
    "Control",
    "PatientSample"
  )
  metadata_df$IsControl <- as.logical(metadata_df$IsControl)
  metadata_df
}

coerce_metadata_schema <- function(metadata_df) {
  if ("SampleName" %in% names(metadata_df)) {
    metadata_df$SampleName <- as.character(metadata_df$SampleName)
  }
  if ("SampleID" %in% names(metadata_df)) {
    metadata_df$SampleID <- as.character(metadata_df$SampleID)
  }
  if ("SampleType" %in% names(metadata_df)) {
    metadata_df$SampleType <- factor(metadata_df$SampleType)
  }
  if ("PatientID" %in% names(metadata_df)) {
    metadata_df$PatientID <- factor(metadata_df$PatientID)
  }

  if ("ProcessingBatch" %in% names(metadata_df)) {
    metadata_df$ProcessingBatch <- factor(metadata_df$ProcessingBatch)
  }
  if ("SequencingBatch" %in% names(metadata_df)) {
    metadata_df$SequencingBatch <- factor(metadata_df$SequencingBatch)
  }
  if ("TumorType" %in% names(metadata_df)) {
    metadata_df$TumorType <- factor(metadata_df$TumorType)
  }
  if ("PatientOutcome" %in% names(metadata_df)) {
    metadata_df$PatientOutcome <- factor(metadata_df$PatientOutcome)
  }
  if ("ControlStatus" %in% names(metadata_df)) {
    metadata_df$ControlStatus <- factor(
      metadata_df$ControlStatus,
      levels = c("Control", "PatientSample")
    )
  }
  metadata_df
}

validate_metadata_df <- function(metadata_df,
                                 required = required_metadata_columns,
                                 context = "metadata") {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), required, context)
  metadata_df <- standardize_metadata_missing_df(metadata_df)
  fail_empty_required_values(metadata_df, required, context)
  if ("PatientID" %in% names(metadata_df)) {
    metadata_df <- drop_inconsistent_patient_metadata(metadata_df, context = context)
  }

  metadata_df <- derive_control_status(metadata_df)
  coerce_metadata_schema(metadata_df)
}

validate_physeq_metadata <- function(physeq,
                                     required = required_physeq_metadata_columns,
                                     context = "phyloseq sample_data") {
  if (is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))) {
    stop("phyloseq object is missing sample_data.", call. = FALSE)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), required, context)
  metadata_df <- standardize_metadata_missing_df(metadata_df)
  fail_empty_required_values(metadata_df, required, context)
  if ("PatientID" %in% names(metadata_df)) {
    metadata_df <- drop_inconsistent_patient_metadata(metadata_df, context = context)
  }

  if (!"ControlStatus" %in% names(metadata_df)) {
    metadata_df <- derive_control_status(metadata_df)
  }
  metadata_df <- coerce_metadata_schema(metadata_df)
  if (!all(phyloseq::sample_names(physeq) %in% rownames(metadata_df))) {
    physeq <- phyloseq::prune_samples(rownames(metadata_df), physeq)
  }
  phyloseq::sample_data(physeq) <- phyloseq::sample_data(metadata_df)

  physeq
}
