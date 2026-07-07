source_if_needed <- function(path) {
  if (file.exists(path)) {
    source(path)
  }
}

source_if_needed(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

sanitize_barplot_path_component <- function(x, fallback = "plot") {
  x <- as.character(x[[1]])
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

normalize_abundance_barplot_config <- function(config = list(), project_config = list()) {
  config <- config %||% list()
  if (is.null(config$norm_method)) config$norm_method <- project_config$norm_method %||% "noNorm"
  if (is.null(config$pseudocount)) config$pseudocount <- project_config$pseudocount %||% 1
  if (is.null(config$tax_agg_level)) config$tax_agg_level <- project_config$tax_agg_level %||% "Genus"
  if (is.null(config$fill_tax_level)) config$fill_tax_level <- config$tax_agg_level
  if (is.null(config$taxa_display)) config$taxa_display <- "top_n"
  if (is.null(config$top_n)) config$top_n <- 20
  if (is.null(config$write_pdf)) config$write_pdf <- FALSE

  if (is.null(config$plots) || length(config$plots) == 0) {
    stop("abundance_barplots must define at least one plot spec.", call. = FALSE)
  }

  config$plots <- lapply(config$plots, function(spec) {
    spec <- spec %||% list()
    if (is.null(spec$name) || !nzchar(trimws(as.character(spec$name)))) {
      stop("Each abundance bar plot spec must define a non-empty name.", call. = FALSE)
    }
    if (is.null(spec$x) || !nzchar(trimws(as.character(spec$x)))) {
      stop("Abundance bar plot spec '", spec$name, "' must define x.", call. = FALSE)
    }
    if (is.null(spec$norm_method)) spec$norm_method <- config$norm_method %||% project_config$norm_method %||% "noNorm"
    if (is.null(spec$pseudocount)) spec$pseudocount <- config$pseudocount %||% project_config$pseudocount %||% 1
    if (is.null(spec$tax_agg_level)) spec$tax_agg_level <- config$tax_agg_level %||% project_config$tax_agg_level %||% "Genus"
    if (is.null(spec$fill_tax_level)) spec$fill_tax_level <- config$fill_tax_level %||% spec$tax_agg_level
    if (is.null(spec$taxa_display)) spec$taxa_display <- config$taxa_display %||% "top_n"
    if (is.null(spec$top_n)) spec$top_n <- config$top_n %||% 20
    if (is.null(spec$write_pdf)) spec$write_pdf <- config$write_pdf %||% FALSE
    spec
  })

  config
}

validate_barplot_metadata_columns <- function(physeq, spec) {
  metadata_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  required <- c(as.character(spec$x))
  if (!is.null(spec$facet) && nzchar(trimws(as.character(spec$facet)))) {
    required <- c(required, as.character(spec$facet))
  }
  if (!is.null(spec$sample_filter) && length(spec$sample_filter) > 0) {
    required <- c(required, names(spec$sample_filter))
  }
  fail_missing_columns(names(metadata_df), unique(required), "phyloseq sample_data")
  invisible(TRUE)
}

resolve_barplot_fill <- function(physeq, fill_tax_level = "Genus") {
  fill_tax_level <- as.character(fill_tax_level[[1]])
  tax_df <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  if (fill_tax_level %in% names(tax_df)) {
    return(fill_tax_level)
  }
  "OTU"
}

collapse_top_taxa <- function(physeq,
                              top_n = 20,
                              other_label = "Other",
                              other_tax_table = NULL) {
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  if (phyloseq::ntaxa(physeq) <= top_n) {
    return(physeq)
  }

  otu_mat <- otu_samples_by_taxa(physeq)
  taxa_totals <- colSums(otu_mat, na.rm = TRUE)
  taxa_order <- names(sort(taxa_totals, decreasing = TRUE))
  top_taxa <- taxa_order[seq_len(min(top_n, length(taxa_order)))]
  other_taxa <- setdiff(colnames(otu_mat), top_taxa)

  other_counts <- rowSums(otu_mat[, other_taxa, drop = FALSE], na.rm = TRUE)
  collapsed_mat <- cbind(otu_mat[, top_taxa, drop = FALSE], other_counts)
  colnames(collapsed_mat)[ncol(collapsed_mat)] <- other_label

  tax_df <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  tax_df <- tax_df[top_taxa, , drop = FALSE]
  if (is.null(other_tax_table)) {
    other_row <- as.list(rep(other_label, ncol(tax_df)))
    names(other_row) <- colnames(tax_df)
  } else {
    other_row <- as.list(other_tax_table)
  }
  tax_df[other_label, ] <- other_row[colnames(tax_df)]

  phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(collapsed_mat), taxa_are_rows = FALSE),
    phyloseq::sample_data(physeq),
    phyloseq::tax_table(as.matrix(tax_df))
  )
}

prepare_abundance_barplot_physeq <- function(physeq, spec) {
  validate_barplot_metadata_columns(physeq, spec)
  physeq <- apply_sample_filter(physeq, spec$sample_filter)
  physeq <- tax_glom_rename(physeq, spec$tax_agg_level)
  physeq <- prune_empty_physeq(physeq)
  physeq <- counts_normalization(
    physeq,
    norm_method = spec$norm_method,
    pseudocount = spec$pseudocount
  )
  physeq <- prune_empty_physeq(physeq)

  taxa_display <- as.character(spec$taxa_display %||% "top_n")
  if (identical(taxa_display, "top_n")) {
    physeq <- collapse_top_taxa(physeq, top_n = spec$top_n)
  } else if (!identical(taxa_display, "all")) {
    stop("Unsupported taxa_display value: ", taxa_display, call. = FALSE)
  }

  validate_barplot_metadata_columns(physeq, spec)
  physeq
}

build_abundance_barplot <- function(physeq, spec) {
  fill_var <- resolve_barplot_fill(physeq, spec$fill_tax_level)
  plot <- phyloseq::plot_bar(
    physeq,
    x = as.character(spec$x),
    fill = fill_var
  ) +
    ggplot2::labs(
      title = spec$plot_title %||% spec$name,
      x = as.character(spec$x),
      y = "Abundance",
      fill = fill_var
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  if (!is.null(spec$facet) && nzchar(trimws(as.character(spec$facet)))) {
    plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", as.character(spec$facet))), scales = "free_x")
  }

  plot
}
