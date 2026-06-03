required_metadata_columns <- c("SampleName", "SampleID", "SampleType", "PatientID")

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
  metadata_df$SampleName <- as.character(metadata_df$SampleName)
  metadata_df$SampleID <- as.character(metadata_df$SampleID)
  metadata_df$SampleType <- factor(metadata_df$SampleType)
  metadata_df$PatientID <- factor(metadata_df$PatientID)

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
  metadata_df$ControlStatus <- factor(
    metadata_df$ControlStatus,
    levels = c("Control", "PatientSample")
  )
  metadata_df
}

validate_metadata_df <- function(metadata_df,
                                 required = required_metadata_columns,
                                 context = "metadata") {
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), required, context)

  for (column in required) {
    empty <- is.na(metadata_df[[column]]) | trimws(as.character(metadata_df[[column]])) == ""
    if (any(empty)) {
      stop(
        context,
        " contains empty value(s) in required column: ",
        column,
        call. = FALSE
      )
    }
  }

  metadata_df <- derive_control_status(metadata_df)
  coerce_metadata_schema(metadata_df)
}

validate_physeq_metadata <- function(physeq,
                                     required = c("SampleID", "SampleType", "PatientID"),
                                     context = "phyloseq sample_data") {
  if (is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))) {
    stop("phyloseq object is missing sample_data.", call. = FALSE)
  }

  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  fail_missing_columns(names(metadata_df), required, context)

  if (!"ControlStatus" %in% names(metadata_df)) {
    metadata_df <- derive_control_status(metadata_df)
    phyloseq::sample_data(physeq) <- phyloseq::sample_data(metadata_df)
  }

  physeq
}

