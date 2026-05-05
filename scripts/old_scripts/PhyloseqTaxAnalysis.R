library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggrepel)
library(argparse)
library(dplyr)
library(taxonomizr)
library(DESeq2)
library(tidyr)
library(ANCOMBC)
library(microbiome)
library(readr)
library(tibble)
library(readxl)
library(stringr)
library(limma)
library(edgeR)
library(purrr)

library(ComplexHeatmap)
library(circlize)


# Keep ggplot from producing Rplots.pdf
if (!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--kraken-file",
                    type = "character",
                    help = "kraken taxonomy classification file path")
parser$add_argument("--dump-dir",
                    type = "character",
                    help = "directory containing nodes.dmp and names.dmp 
                            to be used in database construction")
parser$add_argument("--raw-seq-table",
                    type = "character",
                    help = "un-normalized seq table tsv file path")
parser$add_argument("--norm-method",
                    type = "character",
                    help = "Determine which normalized sequence table to use for limma voom")
parser$add_argument("--bacterial-names",
                    type = "character",
                    help = "names of bacterial ASVs")
parser$add_argument("--library-counts",
                    type = "character",
                    help = "tsv file with read counts at various pipeline stages to use as normalization denominators if needed (e.g. raw read counts, host read counts, etc.)")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "pseudocount value to replace zero counts (default: 1.0)")
parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number for the trial")
parser$add_argument("--add-unclassified-prefix",
                    help = "Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action = "store_true",
                    default = FALSE)
parser$add_argument("--metadata",
                    type = "character",
                    help = "standardized metadata sheet as .xlsx file")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    help = "taxonomic level to agglomerate to for DA analysis (e.g. Genus, Family, etc.)")
parser$add_argument("--tax-label-level",
                    type = "character",
                    help = "taxonomic level to use for taxon labels in plots (e.g. Genus, Family, etc.)")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()


#----------------------------------------
# Positive Bacterial ASVID Extraction
#----------------------------------------
bacterial_IDs <- readLines(args$bacterial_names) # split IDs by line # nolint
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines # nolint


#-----------------------------------------
# Sequence Table Construction
#-----------------------------------------
# Need to remove S## label from sample names in the sequence table to match with metadata sample names
# UN-NORMALIZED
raw_seq_table <- read.delim(args$raw_seq_table,
                            header = TRUE,
                            row.names = 1)
raw_seq_table <- raw_seq_table[, bacterial_IDs] |>
  rownames_to_column(var = "Sample_ID") |>
  mutate(Sample_ID = sub("_S\\d+$", "", Sample_ID)) |>
  column_to_rownames(var = "Sample_ID")


#--------------------------------
# Sample Data Table Construction
#--------------------------------
sample_meta_data_df <- read_excel(args$metadata)
sample_meta_data_df <- sample_meta_data_df |>
  rename(Batch = "ProcessingBatch") |>
  mutate(
    IsControl = str_detect(SampleType, "Control"),
    IsControl = as.logical(IsControl),
    PatientID = ifelse(IsControl,
                       NA,
                       parse_number(SampleID)),
    ControlStatus = ifelse(IsControl,
                           "Control",
                           "PatientSample"),
    SampleType = as.factor(SampleType),
    Batch = as.factor(Batch),
    PatientID = as.factor(PatientID)) 

library_counts_df <- read.delim(args$library_counts, sep = "\t", header = TRUE)
library_counts_df <- library_counts_df |>
  mutate(
      SampleName = sub("_S\\d+$", "", SampleID),
      Raw_reads = as.numeric(Raw_reads),
      HostMappedReads = chimera.filtered - HostUnmapped_reads,
      HostMappedReads = as.numeric(HostMappedReads) 
  ) |>
  select(SampleName, Raw_reads = Raw_reads, Host_mapped_reads = HostMappedReads)  

sample_meta_data_df <- sample_meta_data_df |>
  left_join(library_counts_df, by = "SampleName") |>
  column_to_rownames(var = "SampleName")

#------------------------------------
# Kraken Taxonomy Table Construction
#------------------------------------
# Taxonomy Database File Paths
names_dmp <- file.path(args$dump_dir, "names.dmp")
nodes_dmp <- file.path(args$dump_dir, "nodes.dmp")
sql_db <- file.path(args$dump_dir, "tax_db.sqlite")

# Taxonomy Database Construction
if (!file.exists(sql_db)) {
  read.names.sql(names_dmp, sql_db, overwrite = TRUE)
  read.nodes.sql(nodes_dmp, sql_db, overwrite = TRUE)
}


# Import ASVid -> TaxID info from .kraken2 file
kraken_info <- read.delim(args$kraken_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(kraken_info)[1:3] <- c("status", "ASVid", "taxid")

# trimming kraken information to ASVid and TaxID for all with positive classifications
kraken_info <- kraken_info[kraken_info$status == "C" & kraken_info$taxid != 0, c("ASVid", "taxid")]
kraken_info$taxid <- as.integer(kraken_info$taxid)

desiredTaxa <- c("superkingdom",
                 "phylum",
                 "class",
                 "order",
                 "family",
                 "genus",
                 "species")

kraken_tax_df <- getTaxonomy(ids = kraken_info$taxid,
                             sqlFile = sql_db,
                             desiredTaxa = desiredTaxa)
kraken_tax_df <- as.data.frame(kraken_tax_df, stringsAsFactors = FALSE)
kraken_tax_df$ASVid <- kraken_info$ASVid

kraken_tax_df <- rename(kraken_tax_df,
                        Domain  = superkingdom,
                        Phylum  = phylum,
                        Class   = class,
                        Order   = order,
                        Family  = family,
                        Genus   = genus,
                        Species = species)


kraken_tax_matrix <- as.matrix(kraken_tax_df[, c("Domain",
                                                 "Phylum",
                                                 "Class",
                                                 "Order",
                                                 "Family",
                                                 "Genus",
                                                 "Species")])
rownames(kraken_tax_matrix) <- kraken_tax_df$ASVid

# Taxanomic tag added to all classified taxa
taxa_level_prefix_addition <- function(tax_matrix) {

  tax_prefixes <- c(
    "Domain" = "d",
    "Phylum" = "p",
    "Class" = "c",
    "Order" = "o",
    "Family" = "f",
    "Genus" = "g",
    "Species" = "s",
    "Strain" = "st",
    "Substrain" = "sst"
  )
  for (i in seq_len(ncol(tax_matrix))) {
    for (j in seq_len(nrow(tax_matrix))) {
      if (!is.na(tax_matrix[j, i])) {
        tax_matrix[j, i] <- paste0(tax_prefixes[colnames(tax_matrix)[i]],
                                   "__",
                                   tax_matrix[j, i])
      } else {
        tax_matrix[j, i] <- NA
      }
    }
  }

  tax_matrix
}
kraken_tax_matrix <- taxa_level_prefix_addition(kraken_tax_matrix)

# Handling Unclassified Taxa
unclassified_label_progigation <- function(tax_matrix) {
  is_unassigned <- function(x) is.na(x) || x == "" || x == "Unclassified"
  for (i in seq_len(nrow(tax_matrix))) {
    row <- tax_matrix[i, , drop = TRUE]
    assigned_index <- which(!vapply(row, is_unassigned, logical(1)))
    if (length(assigned_index) == 0) next
    lowest_assigned_index <- max(assigned_index)
    lowest_assigned_value <- row[lowest_assigned_index]
    fill_value <- paste0("UC_", lowest_assigned_value)
    if (lowest_assigned_index < ncol(tax_matrix)) {
      for (j in (lowest_assigned_index + 1):ncol(tax_matrix)) {
        j
        if (is_unassigned(row[j])) {
          tax_matrix[i, j] <- fill_value
        }
      }
    }
  }
  tax_matrix
}

if (isTRUE(args$add_unclassified_prefix)) {
  kraken_tax_matrix <- unclassified_label_progigation(kraken_tax_matrix)
}


#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
raw_kraken_phyloseq <- phyloseq(otu_table(raw_seq_table, taxa_are_rows = FALSE),
                                sample_data(sample_meta_data_df),
                                tax_table(kraken_tax_matrix))
save(raw_kraken_phyloseq, file = paste0(args$out, "/", args$trialID, "_raw_kraken_phyloseq.RData"))

#--------------------------------------
# Read Counts Normalization
#--------------------------------------
# TODO: fix factor division normalization to account for zero divisor errors
# optionally could do away with pre-normalization altogther as it complicates limma/voom assumptions
counts_normalization <- function(physeq, 
                                 norm_method, 
                                 pseudocount, 
                                 tax_agg_level) {
  message("Made it to normalization function somehow")
  if (!is.null(args$tax_agg_level)) {
    physeq <- tax_glom(physeq, taxrank = args$tax_agg_level)
    taxa_names(physeq) <- as.character(tax_table(physeq)[, args$tax_agg_level])
  }

    otu_divide_by_sample_factor <- function(physeq, factor_column) {
      sample_factors <- sample_data(physeq)[[factor_column]]

      otu_mat <- as(otu_table(physeq), "matrix")

      if (taxa_are_rows(physeq)) {
          otu_mat <- sweep(otu_mat, 2, sample_factors, FUN = "/")
      } else {
          otu_mat <- sweep(otu_mat, 1, sample_factors, FUN = "/")
      }
      otu_mat <- otu_mat * 1e6  # scaling factor to bring values back to a more interpretable range
      otu_table(physeq) <- otu_table(otu_mat, taxa_are_rows = taxa_are_rows(physeq))
      return(physeq)
  }

  if (norm_method == "noNorm") {
      physeq <- physeq
  } else if (norm_method == "log2") {
      physeq <- transform_sample_counts(physeq, function(x) log2(x + pseudocount))
  } else if (norm_method == "RelAbund") {
      physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
  } else if (norm_method == "RawTSS") {
      physeq <- otu_divide_by_sample_factor(physeq, "Raw_reads")
  } else if (norm_method == "HostMapped") {
      physeq <- otu_divide_by_sample_factor(physeq, "Host_mapped_reads")
  } else {
      message("Unknown normalization method provided for limma voom pre-normlaization, no normalization used.")
  }

  return(t(as(otu_table(physeq), "matrix")))
}

#---------------------------------------
# Grouping Samples by Control Type
#---------------------------------------
#' @export
SampleGrouping <- function(Ungrouped_phyloseq, 
                           GroupingType,
                           tax_agg_level) {
  if (GroupingType == "AllControl") {
    keep_samples <- !is.na(get_variable(Ungrouped_phyloseq, "SampleType"))
  } else if (GroupingType == "NegativeControl") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("NegativeControl", "Tumor", "NormalTissue")
  } else if (GroupingType == "PatientSample") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("Tumor", "NormalTissue")
  }

  #Protecting against comparisons where one or more comparison groups are missing
  if (!any(keep_samples, na.rm = TRUE)) {
    message(paste(GroupingType, "comparison is invalid due to lack of any relevant samples"))
    return(NULL)
  }

  Grouped_phyloseq <- prune_samples(keep_samples, Ungrouped_phyloseq)

  #Remove any taxa that now have total zero count after subsetting
  Grouped_phyloseq <- prune_taxa(taxa_sums(Grouped_phyloseq) > 0, Grouped_phyloseq)
  if (!is.null(tax_agg_level)) {
    Glommed_phyloseq <- tax_glom(Grouped_phyloseq,
                                 taxrank = tax_agg_level)
    taxa_names(Glommed_phyloseq) <- as.character(tax_table(Glommed_phyloseq)[, tax_agg_level])
  } else {
    Glommed_phyloseq <- Grouped_phyloseq
  }

  return(Glommed_phyloseq)
}


#----------------------------------------------
# ANCOMBC Differential Abundance Analysis
#----------------------------------------------
#' @export
ANCOMBC_DA <- function(Grouped_phyloseq, 
                       GroupingType, 
                       tax_agg_level, 
                       alpha, lfc_cutoff) {

  #STEP 1: Set formula for ANCOM-BC run, desired columns from ANCOM-BC output, and fomatted output file label
  #         based on comparison
  # TODO: simplify handling of control comparisons versus patient sample comparison (single grouping variable)
  if (GroupingType %in% c("NegativeControl", "AllControl")) {
    formula <- "ControlStatus"
    grouping_variable <- "ControlStatus"
    results_table_groups <- c("taxon", "ControlStatusPatientSample")
    default_struc0_groups <- c("taxon", 
                               "structural_zero (ControlStatus = Control)", 
                               "structural_zero (ControlStatus = PatientSample)")
    comp_file_label <- paste0(GroupingType, "toPS")
    
  } else if (GroupingType == "PatientSample") {
    formula <- "SampleType + PatientID"
    grouping_variable <- "SampleType"
    results_table_groups <- c("taxon", "SampleTypeTumor")
    default_struc0_groups <- c("taxon", 
                               "structural_zero (SampleType = NormalTissue)", 
                               "structural_zero (SampleType = Tumor)")
    comp_file_label <- "NTtoT"
  }

  # STEP 2: Run ANCOM-BC analysis
  ancombc_output <- ancombc(data = Grouped_phyloseq, 
                            tax_level = tax_agg_level,
                            formula = formula,
                            p_adj_method = "holm",
                            group = grouping_variable,
                            struc_zero = TRUE,
                            alpha = alpha)

  
  # STEP 3: Construct formatted ANCOM-BC results data frame
  ancombc_lfc_df <- ancombc_output[["res"]][["lfc"]][, results_table_groups]
  colnames(ancombc_lfc_df) <- c("taxon", "log2FoldChange")

  ancombc_p_df <- ancombc_output[["res"]][["p_val"]][, results_table_groups]
  colnames(ancombc_p_df) <- c("taxon", "p")

  ancombc_padj_df <- ancombc_output[["res"]][["q_val"]][, results_table_groups]
  colnames(ancombc_padj_df) <- c("taxon", "padj")

  ancombc_se_df <- ancombc_output[["res"]][["se"]][, results_table_groups]
  colnames(ancombc_se_df) <- c("taxon", "se")
  
  ancombc_struc0_df <- ancombc_output[["zero_ind"]][, default_struc0_groups]
  colnames(ancombc_struc0_df) <- c("taxon", "struc0_group1", "struc0_group2")

  ancombc_struc0_df <- ancombc_struc0_df |>
    mutate(
      struc0_group1 = as.logical(struc0_group1),
      struc0_group2 = as.logical(struc0_group2),
      struc0 = case_when(
        struc0_group1 & !struc0_group2 ~ "group1",
        !struc0_group1 & struc0_group2 ~ "group2",
        TRUE ~ NA
      )
    ) |> select(taxon, struc0)
  
  ancombc_results_df_list <- list(
    ancombc_lfc_df,
    ancombc_p_df,
    ancombc_padj_df,
    ancombc_se_df,
    ancombc_struc0_df
  )
  
  ancombc_results_df <- ancombc_results_df_list |>
    reduce(left_join, by = "taxon")


 
  
  # STEP 4: Add significance labels
  ancombc_results_df <- ancombc_results_df |>
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      p = as.numeric(p),
      padj = as.numeric(padj), 
      se = as.numeric(se),
      significance = if_else(
        padj < alpha &
          (abs(log2FoldChange) > lfc_cutoff |
             !is.na(struc0)),
        "Sig",
        "NotSig"),
      direction = case_when(
        log2FoldChange > 0 | struc0 == "group1" ~ "pos",
        log2FoldChange < 0 | struc0 == "group2" ~ "neg",
        log2FoldChange == 0 ~ "none")) |>
    select(taxon, log2FoldChange, p, padj, struc0, se, significance, direction)


  # STEP 5: Export results
  # TODO: shift this responsibility to snakemake
  if (!dir.exists(paste0(args$out, "/ANCOMBC/", comp_file_label))) {
    dir.create(paste0(args$out, "/ANCOMBC/", comp_file_label),
               recursive = TRUE)
  }
  
  write.table(
    ancombc_results_df,
    file = paste0(args$out, "/ANCOMBC/", comp_file_label, "/", 
                  args$trialID, "_", comp_file_label, "_ANCOMBCResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  # STEP 7: Return all ANCOM-BC results
  ancombc_results_df
  
}

#' @export
LIMMA_VOOM_DA <- function(Grouped_phyloseq,
                          GroupingType,
                          norm_method, psuedocount, tax_agg_level,
                          alpha, lfc_cutoff,
                          lib_size_cutoff = 1) {
  
  metadata <- as.data.frame(as.matrix(sample_data(Grouped_phyloseq)))
  counts <- as.matrix(t(as(otu_table(Grouped_phyloseq), "matrix")))

  # Filtering out samples that are low count to avoid divisor normalization issues
  fltrd_lib_sizes <- colSums(counts)
  low_count_samples <- names(fltrd_lib_sizes[fltrd_lib_sizes < lib_size_cutoff])
  
  if (length(low_count_samples) > 0) {
    warning("Removing extremely low-count samples: ",
            paste(low_count_samples, collapse = ", "))
    counts <- counts[, !(colnames(counts) %in% low_count_samples), drop = FALSE]
    metadata <- metadata[!(rownames(metadata) %in% low_count_samples), , drop = FALSE]
  }

  if (GroupingType %in% c("NegativeControl", "AllControl")) {
    design <- model.matrix(~ ControlStatus,
                           data = metadata)
    results_table_header <- "ControlStatusPatientSample"
    comp_file_label <- paste0(GroupingType, "toPS")
  } else if (GroupingType == "PatientSample") {
    design <- model.matrix(~ SampleType + PatientID,
                          data = metadata)
    results_table_header <- "SampleTypeTumor"
    comp_file_label <-  "NTtoT"
  }

  if (norm_method != "noNorm") {
    # Currently not functioning properly, need to figure out how to properly filter out 
    # samples that result in zero divisor issues
    # Plus looking to eleminate pre-normalization all together if possible
    counts <- as.matrix(counts_normalization(
      physeq = Grouped_phyloseq,
      norm_method = norm_method,
      pseudocount = pseudocount,
      tax_agg_level = tax_agg_level
    ))
    # Filtering any samples that were removed due to zero divisor problems from metadata
    metadata <- metadata[colnames(counts), ]
  } else {
  # Default to limma/voom's recommended normalization protocol 
  # when not trying to force out own normalization method through limma/voom
  # fits better with their assumptions for DA analysis
    counts <- DGEList(counts = counts)
    keep <- filterByExpr(counts, design)
    counts <- counts[keep, ,keep.lib.size=FALSE]
    counts <- calcNormFactors(counts,
                              method = "TMM")
  }

      # TODO: shift this responsibility to snakemake
  if (!dir.exists(paste0(args$out, "/LIMMA_VOOM/", comp_file_label))) {
    dir.create(paste0(args$out, "/LIMMA_VOOM/", comp_file_label),
               recursive = TRUE)
  }

  png(paste0(args$out, "/LIMMA_VOOM/", comp_file_label, "/",
             args$trialID, "_", comp_file_label, "_", args$norm_method, "_VoomVariance.png"), width = 800, height = 600)
  voom <- voom(counts, design, plot=TRUE)
  dev.off()

  fit <-lmFit(voom, design)
  fit <- eBayes(fit)

  limma_voom_results_df <- topTable(fit,
                                    coef= results_table_header,
                                    number= nrow(counts))
  
  limma_voom_results_df <- limma_voom_results_df |>
    rownames_to_column(var="taxon") |>
    mutate(
      log2FoldChange = as.numeric(logFC),
      padj = as.numeric(adj.P.Val), 
      significance = if_else(
        (padj < alpha &
          abs(log2FoldChange) > lfc_cutoff),
        "Sig",
        "NotSig"),
      direction = case_when(
        log2FoldChange > 0 ~ "pos",
        log2FoldChange < 0 ~ "neg",
        log2FoldChange == 0 ~ "none"
      )) |> 
      select(taxon, log2FoldChange, padj, significance, direction)


  
  write.table(
    limma_voom_results_df,
    file = paste0(args$out, "/LIMMA_VOOM/", comp_file_label, "/", 
                  args$trialID, "_", comp_file_label, "_LIMMA_VOOMResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )

  limma_voom_results_df
}

#' @export
label_formatting <- function(DA_results_df,
                             Grouped_phyloseq,
                             tax_label_level) {
  
  DA_results_df$ASVid <- DA_results_df$taxon

  tax_df <- as.data.frame(as(tax_table(Grouped_phyloseq), "matrix")) |> 
    rownames_to_column(var = "ASVid")

  DA_results_df <- left_join(DA_results_df, tax_df, by = "ASVid")

  DA_results_df$label <- ifelse(is.na(DA_results_df[[tax_label_level]]) | 
                                        DA_results_df[[tax_label_level]] == "",
                                      paste0("UC (", DA_results_df$ASVid, ")"),
                                      paste0(DA_results_df[[tax_label_level]], 
                                            "(", DA_results_df$ASVid, ")"))
  DA_results_df
}

#--------------------------------------------
# Volcano Plotting
#--------------------------------------------
#' @export
DA_volcano_plotting <- function(DA_results_df, 
                                GroupingType, 
                                alpha, lfc_cutoff, 
                                DA_method) {
  
  #Subset to significant ASVs
  sig_DA_results_df <- subset(DA_results_df,
                              significance == "Sig" & !is.na(padj))
  pos_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange > 0 & abs(log2FoldChange) > lfc_cutoff)
  neg_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange < 0 & abs(log2FoldChange) > lfc_cutoff)


  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("CellControl", "NegativeControl", "AllControl")) {
    comp_file_label <- paste0(GroupingType, "toPS")
    plot_title_label <- paste(GroupingType, "vs Patient Samples")
  } else if (GroupingType == "PatientSample"){
    comp_file_label <- "NTtoT"
    plot_title_label <- paste("Normal Tissue vs Tumor")
  }
  
  # Build ggplot object
  p_volcano <- ggplot(DA_results_df,
                      aes(x = log2FoldChange,
                          y = -log10(padj))) +
    theme(
      axis.title.x  = element_text(size = 21),
      axis.title.y  = element_text(size = 21),
      plot.title    = element_text(size = 25)
    ) +
    geom_point(alpha = 0.6, size = 4, color = "grey40") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkred") +
    labs(title = paste("Volcano Plot:", plot_title_label),
         x = "Effect size: log2(Fold Change)",
         y = "-log10(adjusted p-value)")
  
  if (nrow(pos_sig_DA_results_df) > 0) {
    p_volcano <- p_volcano +
      geom_point(data = pos_sig_DA_results_df,
                aes(x = log2FoldChange, y = -log10(padj)),
                color = "firebrick1",
                size = 5) +
      geom_text_repel(data = pos_sig_DA_results_df,
                      aes(x = log2FoldChange, y = -log10(padj), label = label),
                      size = 6, 
                      color = "firebrick1",
                      force = 3,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      point.padding = 0.4,
                      # segment(leader) line style
                      segment.size = 0.5,
                      segment.color = "grey40",
                      segment.alpha = 0.9,
                      min.segment.length = 0.5)
  }

  if (nrow(neg_sig_DA_results_df) > 0) {
    p_volcano <- p_volcano +
      geom_point(data = neg_sig_DA_results_df,
                aes(x = log2FoldChange, y = -log10(padj)),
                color = "dodgerblue1",
                size = 5) +
      geom_text_repel(data = neg_sig_DA_results_df,
                      aes(x = log2FoldChange, y = -log10(padj), label = label),
                      size = 6,
                      color = "dodgerblue1", 
                      force = 5,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      point.padding = 0.4,
                      # segment(leader) line style
                      segment.size = 0.5,
                      segment.color = "grey40",
                      segment.alpha = 0.9,
                      min.segment.length = 0.5)
  }

  
  # Save plot
  ggsave(
    filename = paste0(args$out, "/", DA_method, "/", comp_file_label, "/", 
                      args$trialID, "_", comp_file_label, "_", DA_method, "Volcano.png"),
    plot = p_volcano,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
  )
}



#--------------------------------------------
# Heatmap Plotting
#--------------------------------------------
# TODO: Change to standardize the normalization we use for visualization
#' @export
DA_heatmap_plotting <- function(Grouped_phyloseq,
                                DA_results_df,
                                GroupingType,
                                tax_agg_level, tax_label_level,
                                alpha, lfc_cutoff,
                                DA_method) {

  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("CellControl", "NegativeControl", "AllControl")) {
    comp_file_label <- paste0(GroupingType, "toPS")
    plot_title_label <- paste(GroupingType, "vs Patient Samples")
  } else if (GroupingType == "PatientSample"){
    comp_file_label <- "NTtoT"
    plot_title_label <- paste("Normal Tissue vs Tumor")
  }

  sig_DA_results_df <- subset(DA_results_df, significance == "Sig")                              
  sig_taxa <- sig_DA_results_df[, "taxon"]

  if (length(sig_taxa) == 0) {
    return(NULL)
  }

  Pruned_phyloseq <- prune_taxa(sig_taxa, Grouped_phyloseq)

  if (DA_method == "ANCOMBC") {
    row_order <- sig_DA_results_df$taxon[order(sig_DA_results_df$log2FoldChange, sig_DA_results_df$struc0)]
  } else {
    row_order <- sig_DA_results_df$taxon[order(sig_DA_results_df$log2FoldChange)]
  }

  column_order <- rownames(sample_data(Pruned_phyloseq))[order(sample_data(Pruned_phyloseq)$SampleType,
                                                               sample_data(Pruned_phyloseq)$PatientID)]

  p_heatmap <- plot_heatmap(Pruned_phyloseq,
                            method = NULL,
                            sample.order = column_order,
                            sample.label = "SampleID",
                            taxa.order = row_order,
                            taxa.label = tax_label_level,
                            low="#000033", high="#FF3300", na.value = "black")
  
  # Adding vertical lines to heatmap to separate sample types
  type_ordered <- sample_data(Pruned_phyloseq)[column_order, "SampleType"]
  type_block_sizes <- table(type_ordered)
  type_block_cuts <- cumsum(type_block_sizes)

  p_heatmap <- p_heatmap + geom_vline(xintercept = type_block_cuts+0.5, linewidth = 1, color = "white")

  # Adding horizontal lines to heatmap to separate taxa by direction of differential abundance
  neg_lfc_block_size <- sum(sig_DA_results_df$direction == "neg")
  
  if (neg_lfc_block_size > 0 & neg_lfc_block_size < length(row_order)) {
    p_heatmap <- p_heatmap + geom_hline(yintercept = neg_lfc_block_size + 0.5, linewidth = 1, color = "white", linetype = "dashed")
  }

  ggsave(
    filename = paste0(args$out, "/", DA_method, "/", comp_file_label, "/", 
                  args$trialID, "_", comp_file_label, "_", DA_method, "Heatmap.png"),
    plot = p_heatmap,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
  )
}




#-----------------------------------------
# Differential Analysis Wrapper
#-----------------------------------------
#' @export
DAxPlottingWrapper <- function(Grouped_phyloseq, 
                               DA_method,
                               GroupingType,
                               norm_method,
                               tax_agg_level = "Genus", tax_label_level = "Genus",
                               alpha = 0.01, lfc_cutoff = 1) {
  
  if (DA_method == "ANCOMBC") {
    DA_results_df <- ANCOMBC_DA(Grouped_phyloseq = Grouped_phyloseq,
                                GroupingType = GroupingType,
                                tax_agg_level = tax_agg_level,
                                alpha = alpha, lfc_cutoff = lfc_cutoff)
  } else if (DA_method == "LIMMA_VOOM"){
    DA_results_df <- LIMMA_VOOM_DA(Grouped_phyloseq = Grouped_phyloseq,
                                   GroupingType = GroupingType,
                                   norm_method = norm_method,
                                   tax_agg_level = tax_agg_level,
                                   alpha = alpha, lfc_cutoff = lfc_cutoff)
  }

  # Adding desired taxa level label to have inaddition to ASV's 
  # (i.e. if tax_agg_level is NULL)
  if (is.null(tax_agg_level)) {
    DA_results_df <- label_formatting(
      DA_results_df = DA_results_df,
      Grouped_phyloseq = Grouped_phyloseq,
      tax_label_level = tax_label_level
    )
  } else {
    DA_results_df$label <- DA_results_df$taxon
  }

  DA_volcano_plotting(DA_results_df = DA_results_df,
                      GroupingType = GroupingType,
                      alpha = alpha, lfc_cutoff = lfc_cutoff,
                      DA_method = DA_method)
  
  DA_heatmap_plotting(Grouped_phyloseq = Grouped_phyloseq,
                      GroupingType = GroupingType,
                      DA_results_df = DA_results_df,
                      tax_agg_level = tax_agg_level, tax_label_level = tax_label_level,
                      alpha = alpha, lfc_cutoff = lfc_cutoff,
                      DA_method = DA_method)

}


#----------------------------------
# Differential Abundance Execution
#----------------------------------
#' @export
all_DA_analysis <- function(phyloseq,
                            DA_methods,
                            norm_method,
                            Comparisons,
                            tax_agg_level = NULL, tax_label_level = "Genus") {
  
  for (DA_method in DA_methods) {
    for (Comparison in Comparisons) {
      Grouped_phyloseq <- SampleGrouping(Ungrouped_phyloseq = phyloseq, 
                                         GroupingType = Comparison,
                                         tax_agg_level = tax_agg_level)
      if (nsamples(Grouped_phyloseq) == 0 ||
          is.null(sample_data(Grouped_phyloseq, errorIfNULL = FALSE)) ||
          length(unique(sample_data(Grouped_phyloseq)[,"SampleType"])) < 2) {
        message(paste("Skipping", DA_method, Comparison, "due to lack of one or more sample-type group(s)"))
        next
      }

      DAxPlottingWrapper(
        Grouped_phyloseq = Grouped_phyloseq,
        DA_method = DA_method,
        norm_method = norm_method,
        tax_agg_level = tax_agg_level, tax_label_level = tax_label_level,
        GroupingType = Comparison)
    }
  }
}


all_DA_analysis(phyloseq = raw_kraken_phyloseq,
                DA_methods = c("ANCOMBC", "LIMMA_VOOM"),
                norm_method = args$norm_method,
                Comparisons = c("NegativeControl", "PatientSample"),
                tax_agg_level = "Genus")
#-------------------------------
# R Session Cataloging
#-------------------------------
# allows for more accessible downstream exploratory data analysis
save.image(file = paste0(args$out, "/", args$trialID, "_PhyloSeqSession.RData"))


#--------------------------------------
# Abundance Table and Plot
#--------------------------------------
# GenusAbundance_tableXplot <- function(norm_phyloseq) {
#   # STEP 1: Genus Glomming
#   genusGlom_phyloseq <- tax_glom(norm_phyloseq,
#                                  taxrank = "Genus",
#                                  NArm = TRUE)
#   top10genus_phyloseq <- prune_taxa(names(sort(taxa_sums(genusGlom_phyloseq),
#                                                decreasing = TRUE))[1:10],
#                                     genusGlom_phyloseq)

#   # STEP 2: Genus Table Construction
#   # constructing genus name matrix from tax table of genus glommed phyloseq object
#   genus_name_matrix <- as(tax_table(genusGlom_phyloseq), "matrix")[, "Genus"]

#   # setting taxa names in phyloseq object to genus level name matrix
#   taxa_names(genusGlom_phyloseq) <- genus_name_matrix

#   # Constructing genus count table from glommed OTU table
#   genus_matrix <- as.matrix(otu_table(genusGlom_phyloseq))

#   if (taxa_are_rows(genusGlom_phyloseq)) {
#     genus_matrix <- t(genus_matrix)
#   }

#   # Add total row
#   genus_matrix <- rbind(genus_matrix, Total = colSums(genus_matrix))

#   write.table(
#     genus_matrix,
#     file = paste0(args$out, "/", args$trialID, "_kraken_genus_table.tsv"),
#     sep = "\t",
#     quote = FALSE,
#     row.names = TRUE,
#     col.names = NA
#   )

#   # STEP 3: Count Plot
#   theme_set(theme_bw())

#   abunXtypeXsample_plot <- plot_bar(top10genus_phyloseq, "SampleID", fill = "Genus")
#   abunXtypeXsample_plot +
#     theme(
#       legend.position = "bottom",
#       legend.text = element_text(size = 6)
#     ) +
#     facet_wrap(~SampleType, scales = "free_x") +
#     labs(title = "Read Counts per Sample by Sample Type",
#          x = "Sample",
#          y = "Read Count")
#   ggsave(paste0(args$out, "/", args$trialID, "_kraken_countXtypeXsample.png"), width = 14, height = 12, units = "in")
# }

# GenusAbundance_tableXplot(norm_kraken_phyloseq)