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

parser <- ArgumentParser()

parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number to attach to output files")
parser$add_argument("--B1-physeq",
                    type = "character",
                    help = ".RData file contatining phyloseq object for required first batch")
parser$add_argument("--B2-physeq",
                    type = "character",
                    default = NULL,
                    help = ".RData file contatining phyloseq object for optional second batch")
parser$add_argument("--DA-methods",
                    type = "character",
                    nargs = "+",
                    help = "Choice of differential abunance tools to use [ANCOMBC, LIMMA_VOOM]")
parser$add_argument("--DA-comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Choice of differential abunance comprisons to use [CellLineControltoTumor, CellLineControltoNontumor, NegativeControl, PatientSample, TumorType]")
parser$add_argument("--norm-method",
                    type = "character",
                    default = "noNorm",
                    help = "Determine which normalized sequence table to feed into limma voom and scale heatmaps values (default = noNorm)")
parser$add_argument("--pseudocount",
                    type = "double",
                    default = 1.0,
                    help = "pseudocount value to replace zero counts (default: 1.0)")
parser$add_argument("--tax-agg-level",
                    type = "character",
                    default = "Genus",
                    help = "taxonomic level to agglomerate to for DA analysis (e.g. Genus, Family, etc.)")
parser$add_argument("--tax-label-level",
                    type = "character",
                    default = "Genus",
                    help = "taxonomic level to use for taxon labels in plots (e.g. Genus, Family, etc.)")
parser$add_argument("--select-taxa-names",
                    type = "character",
                    nargs = "+",
                    help = "optional list of files containing lists of taxa to subset differential abundance analysis to (union of all unique names provided)")
parser$add_argument("--alpha",
                    type = "double",
                    default = 0.05,
                    help = "Alpha cutoff for differential abundance significance (default = 0.05)")
parser$add_argument("--lfc-cutoff",
                    type = "double",
                    default = 1.0,
                    help = "Log-fold change cutoff for differential abundance significance (default = 1.0)")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()

#--------------------------------------
# Read Counts Normalization
#--------------------------------------
# TODO: fix factor division normalization to account for zero divisor errors
# optionally could do away with pre-normalization altogther as it complicates limma/voom assumptions
counts_normalization <- function(physeq, 
                                 norm_method, 
                                 pseudocount) {


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
        return(physeq)
    } else if (norm_method == "log2") {
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else if (norm_method == "RelAbund") {
        return(transform_sample_counts(physeq, function(x) x / sum(x)))
    } else if (norm_method == "RawTSS") {
        return(otu_divide_by_sample_factor(physeq, "Raw_reads"))
    } else if (norm_method == "HostMapped") {
        return(otu_divide_by_sample_factor(physeq, "Host_mapped_reads"))
    } else if (norm_method == "log2HostMapped") {
        physeq <- otu_divide_by_sample_factor(physeq, "Host_mapped_reads")
        return(transform_sample_counts(physeq, function(x) log2(x + pseudocount)))
    } else {
        message("Unknown normalization method provided for limma voom pre-normlaization, no normalization used.")
    }

}

#---------------------------------------
# Grouping Samples by Control Type
#---------------------------------------
SampleGrouping <- function(Ungrouped_phyloseq, 
                           GroupingType,
                           tax_agg_level) {
  if (GroupingType == "AllControl") {
    keep_samples <- !is.na(get_variable(Ungrouped_phyloseq, "SampleType"))
  } else if (GroupingType == "TumorType") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") == "Tumor"
  } else if (GroupingType == "CellLineControltoTumor") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("CellLineControl", "Tumor")
  } else if (GroupingType == "CellLineControltoNontumor") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("CellLineControl", "Nontumor", "NormalTissue")
  } else if (GroupingType == "NegativeControl") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("NegativeControl", "Tumor", "Nontumor", "NormalTissue")
  } else if (GroupingType == "PatientSample") {
    keep_samples <- get_variable(Ungrouped_phyloseq, "SampleType") %in% c("Tumor", "Nontumor", "NormalTissue")
  }

  #Protecting against comparisons where one or more comparison groups are missing
  if (!any(keep_samples)) {
    return(NULL)
  }

  Grouped_phyloseq <- prune_samples(keep_samples, Ungrouped_phyloseq)
  Grouped_phyloseq <- tax_glom(Grouped_phyloseq, tax_agg_level)
  taxa_names(Grouped_phyloseq) <- as.character(tax_table(Grouped_phyloseq)[, tax_agg_level])

  #Remove any taxa that now have total zero count after subsetting and glomming
  Grouped_phyloseq <- prune_taxa(taxa_sums(Grouped_phyloseq) > 0, Grouped_phyloseq)

  return(Grouped_phyloseq)
}


#----------------------------------------------
# ANCOMBC Differential Abundance Analysis
#----------------------------------------------
ANCOMBC_DA <- function(Grouped_phyloseq, 
                       GroupingType, 
                       tax_agg_level, 
                       alpha, lfc_cutoff) {

  #STEP 1: Set formula for ANCOM-BC run, desired columns from ANCOM-BC output, and fomatted output file label
  #         based on comparison
  # TODO: simplify handling of control comparisons versus patient sample comparison (single grouping variable)
  if (GroupingType %in% c("CellLineControltoTumor", "CellLineControltoNontumor", "NegativeControl", "AllControl")) {
    formula <- "ControlStatus"
    grouping_variable <- "ControlStatus"
    results_table_groups <- c("taxon", "ControlStatusPatientSample")
    default_struc0_groups <- c("taxon", 
                               "structural_zero (ControlStatus = Control)", 
                               "structural_zero (ControlStatus = PatientSample)")
    comp_file_label <- paste0(GroupingType, "toPS")

  } else if (GroupingType == "TumorType") {
    formula <- "TumorType"
    grouping_variable <- "TumorType"
    results_table_groups <- c("taxon", "TumorTypeiCC")
    default_struc0_groups <- c("taxon",
                               "structural_zero (TumorType = HCC)",
                               "structural_zero (TumorType = iCC)")
    comp_file_label <- "TumorType"

  } else if (GroupingType == "PatientSample") {
    formula <- "SampleType + PatientID"
    grouping_variable <- "SampleType"
    results_table_groups <- c("taxon", "SampleTypeTumor")

    if ("NormalTissue" %in% levels(as(sample_data(Grouped_phyloseq), "data.frame")$SampleType)) {
      default_struc0_groups <- c("taxon", 
                                "structural_zero (SampleType = NormalTissue)", 
                                "structural_zero (SampleType = Tumor)")
    } else if ("Nontumor" %in% levels(as(sample_data(Grouped_phyloseq), "data.frame")$SampleType)) {
            default_struc0_groups <- c("taxon", 
                                "structural_zero (SampleType = Nontumor)", 
                                "structural_zero (SampleType = Tumor)")
    }

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
  if (!dir.exists(paste0(args$out, "/ANCOMBC/", GroupingType))) {
    dir.create(paste0(args$out, "/ANCOMBC/", GroupingType),
               recursive = TRUE)
  }
  
  write.table(
    ancombc_results_df,
    file = paste0(args$out, "/ANCOMBC/", GroupingType, "/", 
                  args$trialID, "_", GroupingType, "_ANCOMBCResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  # STEP 7: Return all ANCOM-BC results
  ancombc_results_df
  
}

LIMMA_VOOM_DA <- function(Grouped_phyloseq,
                          GroupingType,
                          norm_method, pseudocount, 
                          tax_agg_level,
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

  if (GroupingType %in% c("CellLineControltoTumor", "CellLineControltoNontumor", "NegativeControl", "AllControl")) {
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
  if (!dir.exists(paste0(args$out, "/LIMMA_VOOM/", GroupingType))) {
    dir.create(paste0(args$out, "/LIMMA_VOOM/", GroupingType),
               recursive = TRUE)
  }

  png(paste0(args$out, "/LIMMA_VOOM/", GroupingType, "/",
             args$trialID, "_", GroupingType, "_", args$norm_method, "_VoomVariance.png"), width = 800, height = 600)
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
    file = paste0(args$out, "/LIMMA_VOOM/", GroupingType, "/", 
                  args$trialID, "_", GroupingType, "_LIMMA_VOOMResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )

  limma_voom_results_df
}

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
DA_volcano_plotting <- function(DA_results_df, 
                                GroupingType, 
                                alpha, lfc_cutoff, 
                                DA_method) {
  #Subset to non-structural zeros
  nostruc0_DA_results_df <- DA_results_df |> subset(is.na(struc0))
  #Subset to significant ASVs
  sig_DA_results_df <- subset(nostruc0_DA_results_df,
                              significance == "Sig" & !is.na(padj))
  pos_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange > 0 & abs(log2FoldChange) > lfc_cutoff)
  neg_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange < 0 & abs(log2FoldChange) > lfc_cutoff)


  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("TumorType", "CellLineControltoTumor", "CellLineControltoNontumor")) {
    comp_file_label <- GroupingType
    plot_title_label <- GroupingType
  } else if (GroupingType %in% c("NegativeControl", "AllControl")) {
    comp_file_label <- paste0(GroupingType, "toPS")
    plot_title_label <- paste(GroupingType, "vs Patient Samples")
  } else if (GroupingType == "PatientSample"){
    comp_file_label <- "NTtoT"
    plot_title_label <- paste("Normal Tissue vs Tumor")
  }
  
  # Build ggplot object
  p_volcano <- ggplot(nostruc0_DA_results_df, #exclude structural zeros (make no sense in volcano plot)
                      aes(x = log2FoldChange,
                          y = -log10(padj))) +
    theme(
      axis.title.x  = element_text(size = 21),
      axis.title.y  = element_text(size = 21),
      plot.title    = element_text(size = 25),
      axis.text     = element_text(size = 18)
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
    filename = paste0(args$out, "/", DA_method, "/", GroupingType, "/", 
                      args$trialID, "_", GroupingType, "_", DA_method, "Volcano.png"),
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
DA_heatmap_plotting <- function(Grouped_phyloseq,
                                DA_results_df,
                                GroupingType,
                                norm_method, pseudocount,
                                tax_agg_level, tax_label_level,
                                alpha, lfc_cutoff,
                                DA_method) {

  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("TumorType", "CellLineControltoTumor", "CellLineControltoNontumor")) {
    comp_file_label <- GroupingType
    plot_title_label <- GroupingType
  } else if (GroupingType %in% c("NegativeControl", "AllControl")) {
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

  
  Normalized_phyloseq <-   counts_normalization(physeq = Grouped_phyloseq, 
                                                norm_method = norm_method, 
                                                pseudocount = pseudocount)

  Pruned_phyloseq <- prune_taxa(sig_taxa, Normalized_phyloseq)
  Pruned_phyloseq <- prune_samples(sample_sums(Pruned_phyloseq) > 0, Pruned_phyloseq)

  if (DA_method == "ANCOMBC" & length(sig_taxa) > 1) {
    row_order <- sig_DA_results_df$taxon[order(sig_DA_results_df$log2FoldChange, sig_DA_results_df$struc0)]
  } else if (DA_method == "LIMMA_VOOM" & length(sig_taxa) > 1) {
    row_order <- sig_DA_results_df$taxon[order(sig_DA_results_df$log2FoldChange)]
  } else {
    row_order = NULL
  }

  column_order <- rownames(sample_data(Pruned_phyloseq))[order(sample_data(Pruned_phyloseq)$SampleType,
                                                               sample_data(Pruned_phyloseq)$PatientID)]

  p_heatmap <- plot_heatmap(Pruned_phyloseq,
                            method = "PCoA",
                            distance = "bray",
                            sample.order = column_order,
                            sample.label = "SampleID",
                            #taxa.order = row_order,
                            taxa.label = tax_label_level,
                            low="#000033", high="#FF3300", na.value = "black")
  
  # Adding vertical lines to heatmap to separate sample types
  type_ordered <- sample_data(Pruned_phyloseq)[column_order, "SampleType"]
  type_block_sizes <- table(type_ordered)
  type_block_cuts <- cumsum(type_block_sizes)

  p_heatmap <- p_heatmap + geom_vline(xintercept = type_block_cuts+0.5, linewidth = 1, color = "white")

  # Adding horizontal lines to heatmap to separate taxa by direction of differential abundance
  neg_lfc_block_size <- sum(sig_DA_results_df$direction == "neg")
  
  # if (neg_lfc_block_size > 0 & neg_lfc_block_size < length(row_order)) {
  #   p_heatmap <- p_heatmap + geom_hline(yintercept = neg_lfc_block_size + 0.5, linewidth = 1, color = "white", linetype = "dashed")
  # }

  ggsave(
    filename = paste0(args$out, "/", DA_method, "/", GroupingType, "/", 
                  args$trialID, "_", GroupingType, "_", DA_method, "Heatmap.png"),
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
DAxPlottingWrapper <- function(Grouped_phyloseq, 
                               DA_method,
                               GroupingType,
                               norm_method, pseudocount,
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
                                   norm_method = norm_method, pseudocount = pseudocount,
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
                      norm_method = norm_method, pseudocount = pseudocount,
                      tax_agg_level = tax_agg_level, tax_label_level = tax_label_level,
                      alpha = alpha, lfc_cutoff = lfc_cutoff,
                      DA_method = DA_method)

}


#----------------------------------
# Differential Abundance Execution
#----------------------------------
all_DA_analysis <- function(phyloseq,
                            DA_methods,
                            Comparisons,
                            norm_method, pseudocount,
                            tax_agg_level = NULL, tax_label_level = "Genus",
                            select_taxa = NULL,
                            alpha, lfc_cutoff) {
  
  for (DA_method in DA_methods) {
    for (Comparison in Comparisons) {
      Grouped_phyloseq <- SampleGrouping(Ungrouped_phyloseq = phyloseq, 
                                         GroupingType = Comparison,
                                         tax_agg_level = tax_agg_level)
      if (!is.null(select_taxa)) {
        Grouped_phyloseq <- prune_taxa(select_taxa, Grouped_phyloseq)
      }
      
      if (Comparison == "TumorType") {
        Grouping_col <- "TumorType"
      } else {
        Grouping_col <- "SampleType"
      }

      if (nsamples(Grouped_phyloseq) == 0 ||
          is.null(Grouped_phyloseq) ||
          length(unique(as.matrix(sample_data(Grouped_phyloseq))[, Grouping_col])) < 2) {
        message_string <- paste("Skipped", DA_method, Comparison, "due to lack of one or more sample-type group(s)")
        message(message_string)

        if (!dir.exists(paste0(args$out, "/", DA_method, "/", Comparison))) {
          dir.create(paste0(args$out, "/", DA_method, "/", Comparison),
                    recursive = TRUE)
        }

        writeLines(message_string,
                   paste0(args$out, "/", DA_method, "/", Comparison, "/", 
                   args$trialID, "_", Comparison, "_", DA_method, "Results.tsv"))
        next
      }

      DAxPlottingWrapper(
        Grouped_phyloseq = Grouped_phyloseq,
        DA_method = DA_method,
        norm_method = norm_method, pseudocount = pseudocount,
        tax_agg_level = tax_agg_level, tax_label_level = tax_label_level,
        GroupingType = Comparison,
        alpha = alpha, lfc_cutoff = lfc_cutoff)
    }
  }
}

load(args$B1_physeq)
B1physeq <- physeq

if (!is.null(args$B2_physeq)) {
  load(args$B2_physeq)
  B2physeq <- physeq
  if (!is.null(args$tax_agg_level)) {
    B1physeq <- tax_glom(B1physeq, args$tax_agg_level)
    B2physeq <- tax_glom(B2physeq, args$tax_agg_level)
    taxa_names(B1physeq) <- as.character(tax_table(B1physeq)[, args$tax_agg_level])
    taxa_names(B2physeq) <- as.character(tax_table(B2physeq)[, args$tax_agg_level])
  }
  CompPhyseq <- merge_phyloseq(B1physeq, B2physeq)
} else {
  CompPhyseq <- B1physeq
}
# Read in selected taxa names if chosen to subset
if (!is.null(args$select_taxa_names)) {
  select_taxa <- unique(unlist(lapply(args$select_taxa_names, function(name_file) {
    read.csv(name_file, stringsAsFactors = FALSE)[[1]]
  })))
} else {
  select_taxa = NULL
}


all_DA_analysis(phyloseq = CompPhyseq,
                DA_methods = args$DA_methods,
                Comparisons = args$DA_comparisons,
                norm_method = args$norm_method,
                pseudocount = args$pseudocount,
                tax_agg_level = args$tax_agg_level,
                select_taxa = select_taxa,
                alpha = args$alpha, 
                lfc_cutoff = args$lfc_cutoff)