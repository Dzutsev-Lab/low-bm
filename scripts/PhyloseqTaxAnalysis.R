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
library(readxl)
library(stringr)
library(tibble)

library(ComplexHeatmap)
library(circlize)


# Keep ggplot from producing Rplots.pdf
if (!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--mothur-file",
                    type = "character",
                    help = "taxonomy classification file path")
parser$add_argument("--kraken-file",
                    type = "character",
                    help = "kraken taxonomy classification file path")
parser$add_argument("--dump-dir",
                    type = "character",
                    help = "directory containing nodes.dmp and names.dmp 
                            to be used in database construction")
parser$add_argument("--norm-method",
                    type = "character",
                    help = "normalization method used to produce norm SeqTable")
parser$add_argument("--norm-seq-table",
                    type = "character",
                    help = "normalized seq table tsv file path")
parser$add_argument("--raw-seq-table",
                    type = "character",
                    help = "un-normalized seq table tsv file path")
parser$add_argument("--trialID",
                    type = "character",
                    help = "ID number for the trial")
parser$add_argument("--add-unclassified-prefix",
                    help = "Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action = "store_true",
                    default = FALSE)
parser$add_argument("--metadata",
                    type = "character",
                    help = "standardized metadata sheet as .xslx file")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store output abundance plots")


args <- parser$parse_args()


#-----------------------------------------
# Sequence Table Construction
#-----------------------------------------
# NORMALIZED
norm_seq_table <- read.delim(args$norm_seq_table,
                             header = TRUE,
                             row.names = 1)

# UN-NORMALIZED
raw_seq_table <- read.delim(args$raw_seq_table,
                            header = TRUE,
                            row.names = 1)


#--------------------------------
# Sample Data Table Construction
#--------------------------------
# sample_names <- rownames(norm_seq_table)
# sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)

# # Sample Type
# sample_type <- sub("\\d+[A-Za-z]*$", "", sample_info)

# # Technical Rep (Will also denote patient ID for patient samples)
# tech_rep <- sub(".*?(\\d+[A-Za-z]*)$", "\\1", sample_info)

# sample_metadata_df <- data.frame(
#   SampleType = factor(sample_type),
#   Replicate = factor(tech_rep),
#   row.names = rownames(norm_seq_table),
#   stringsAsFactors = FALSE
# )

# sample_metadata_df <- sample_metadata_df |>
#   mutate(
#     SampleType = case_when(
#                            (is.na(SampleType) | SampleType == "") &
#                              (grepl("NT$", Replicate) | grepl("N$", Replicate))
#                            ~ "NormalTissue", #check before tumor otherwise all would be labeled tumor (ending with T)
#                            (is.na(SampleType) | SampleType == "") &
#                              grepl("T$", Replicate) ~ "Tumor",
#                            TRUE ~ SampleType),
#     SampleType = case_when(
#                            SampleType == "expcontrol" ~ "CellControl",
#                            SampleType == "NEGATIVECONTROL" ~ "NegativeControl",
#                            TRUE ~ as.character(SampleType)),
#     SampleType = factor(SampleType),
#     Replicate = as.character(Replicate),
#     PatientID = case_when(
#                           SampleType %in% c("Tumor", "NormalTissue") ~ parse_number(Replicate),
#                           TRUE ~ NA),
#     PatientID = factor(PatientID)
#   )
sample_metadata_df <- read_excel(args$metadata)
sample_metadata_df <- sample_metadata_df |>
  as.data.frame() |>
  mutate(
    SampleID = as.character(SampleID),
    PatientID = ifelse(!str_detect(SampleType, "Control"),
                       parse_number(SampleID),
                       NA)) |>
  column_to_rownames("SampleName")

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

kraken_tax_matrix <- if (isTRUE(args$add_unclassified_prefix)) {
  unclassified_label_progigation(kraken_tax_matrix)
}

save.image(file = paste0(args$out, "/PhyloseqTaxAnalysis.RData"))
#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
norm_kraken_phyloseq <- phyloseq(otu_table(norm_seq_table, taxa_are_rows = FALSE),
                                 sample_data(sample_metadata_df),
                                 tax_table(kraken_tax_matrix))

raw_kraken_phyloseq <- phyloseq(otu_table(raw_seq_table, taxa_are_rows = FALSE),
                                sample_data(sample_metadata_df),
                                tax_table(kraken_tax_matrix))

#--------------------------------------
# Abundance Table and Plot
#--------------------------------------
GenusAbundance_tableXplot <- function(norm_phyloseq) {
  # STEP 1: Genus Glomming
  genusGlom_phyloseq <- tax_glom(norm_phyloseq,
                                 taxrank = "Genus",
                                 NArm = TRUE)
  top10genus_phyloseq <- prune_taxa(names(sort(taxa_sums(genusGlom_phyloseq),
                                               decreasing = TRUE))[1:10],
                                    genusGlom_phyloseq)

  # STEP 2: Genus Table Construction
  # constructing genus name matrix from tax table of genus glommed phyloseq object
  genus_name_matrix <- as(tax_table(genusGlom_phyloseq), "matrix")[, "Genus"]

  # setting taxa names in phyloseq object to genus level name matrix
  taxa_names(genusGlom_phyloseq) <- genus_name_matrix

  # Constructing genus count table from glommed OTU table
  genus_matrix <- as.matrix(otu_table(genusGlom_phyloseq))

  if (taxa_are_rows(genusGlom_phyloseq)) {
    genus_matrix <- t(genus_matrix)
  }

  # Add total row
  genus_matrix <- rbind(genus_matrix, Total = colSums(genus_matrix))

  write.table(
    genus_matrix,
    file = paste0(args$out, "/", args$trialID, "_kraken_genus_table.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  # STEP 3: Count Plot
  theme_set(theme_bw())

  abunXtypeXsample_plot <- plot_bar(top10genus_phyloseq, "SampleID", fill = "Genus")
  abunXtypeXsample_plot +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 6)
    ) +
    facet_wrap(~SampleType, scales = "free_x") +
    labs(title = "Read Counts per Sample by Sample Type",
         x = "Sample",
         y = "Read Count")
  ggsave(paste0(args$out, "/", args$trialID, "_kraken_countXtypeXsample.png"), width = 14, height = 12, units = "in")
}





#---------------------------------------
# Grouping Samples by Control Type
#---------------------------------------
SampleGrouping <- function(Ungrouped_phyloseq, 
                           GroupingType) {
  
  sample_df <- as(sample_data(Ungrouped_phyloseq), "data.frame")
  
  if (GroupingType %in% c("AllControl", "NegativeControl", "CellControl")) {
    sample_df <- sample_df |>
      mutate(SampleType = ifelse(SampleType %in% c("NormalTissue", "Tumor"),
                                 "PatientSample",
                                 as.character(SampleType)))

    if (GroupingType == "AllControl") {
      sample_df <- sample_df |>
        mutate(SampleType = ifelse(SampleType %in% c("NegativeControl", "CellControl"),
                                   "AllControl",
                                   as.character(SampleType)))
    }

    GroupingLevels <- c(GroupingType, "PatientSample")

  } else if (GroupingType == "PatientSample") {
    GroupingLevels <- c("NormalTissue", "Tumor")
  }
  #ensure SampleType is correctly factored
  sample_df$SampleType <- factor(sample_df$SampleType)

  #prune samples not included in current comparison
  sample_df <- sample_df[sample_df$SampleType %in% GroupingLevels, ]
  select_samples <- rownames(sample_df)
  Grouped_phyloseq <- prune_samples(select_samples, Ungrouped_phyloseq)
  
  #update phyloseq object with refactored sample meta data
  sample_data(Grouped_phyloseq) <- sample_data(sample_df)

  #Remove any taxa that now have total zero count after subsetting
  Grouped_phyloseq <- prune_taxa(taxa_sums(Grouped_phyloseq) > 0, Grouped_phyloseq)

  Grouped_phyloseq
}


#----------------------------------------------
# ANCOMBC Differential Abundance Analysis
#----------------------------------------------
ANCOMBC_DA <- function(Grouped_phyloseq, 
                       GroupingType, 
                       tax_agg_level, 
                       tax_label_level, 
                       alpha, 
                       lfc_cutoff) {

  #STEP 1: Set formula for ANCOM-BC run, desired columns from ANCOM-BC output, and fomatted output file label
  #         based on comparison
  if (GroupingType %in% c("CellControl", "NegativeControl", "AllControl")) {
    formula <- "SampleType"
    lfc_adjp_groups <- c("taxon", "SampleTypePatientSample")
    default_struc0_groups <- c("taxon", 
                               paste0("structural_zero (SampleType = ", GroupingType, ")"), 
                               "structural_zero (SampleType = PatientSample)")
    struc0_group1_label <- GroupingType
    struc0_group2_label <- "PatientSample"
    comp_file_label <- paste0(GroupingType, "vsPatientSample")
    
  } else if (GroupingType == "PatientSample") {
    formula <- "SampleType + PatientID"
    lfc_adjp_groups <- c("taxon", "SampleTypeTumor")
    default_struc0_groups <- c("taxon", 
                               "structural_zero (SampleType = NormalTissue)", 
                               "structural_zero (SampleType = Tumor)")
    struc0_group1_label <- "NormalTissue"
    struc0_group2_label <- "Tumor"
    comp_file_label <- "NormalvsTumor"
  }

  # STEP 2: Run ANCOM-BC analysis
  ancombc_output <- ancombc(data = Grouped_phyloseq, 
                            tax_level = tax_agg_level,
                            formula = formula,
                            p_adj_method = "holm",
                            group = "SampleType",
                            struc_zero = TRUE,
                            alpha = alpha)

  
  # STEP 3: Construct formatted ANCOM-BC results data frame
  ancombc_lfc_df <- ancombc_output[["res"]][["lfc"]][, lfc_adjp_groups]
  colnames(ancombc_lfc_df) <- c("taxon", "log2FoldChange")
  
  ancombc_padj_df <- ancombc_output[["res"]][["q_val"]][, lfc_adjp_groups]
  colnames(ancombc_padj_df) <- c("taxon", "padj")
  
  ancombc_struc0_df <- ancombc_output[["zero_ind"]][, default_struc0_groups]
  colnames(ancombc_struc0_df) <- c("taxon", "struc0_group1", "struc0_group2")
  ancombc_struc0_df <- ancombc_struc0_df |>
    mutate(
      struc0_group1 = as.logical(struc0_group1),
      struc0_group2 = as.logical(struc0_group2),
      struc0 = case_when(
        struc0_group1 & !struc0_group2 ~ struc0_group1_label,
        !struc0_group1 & struc0_group2 ~ struc0_group2_label,
        TRUE ~ NA
      )
    )
  

  
  ancombc_results_df <- ancombc_lfc_df |>
    left_join(ancombc_padj_df, by = "taxon") |>
    left_join(ancombc_struc0_df[, c("taxon", "struc0")], by = "taxon")
 
  
  # STEP 4: Add significance labels
  ancombc_results_df <- ancombc_results_df |>
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj), 
      Significance = if_else(
        padj < alpha &
          (abs(log2FoldChange) > lfc_cutoff |
             !is.na(struc0)),
        "Sig",
        "NotSig"
      )
    )

  # STEP 5: Add taxon assignments labels to ANCOM-BC results data frame 
  #         only needed if no taxa level aggregation done in ANCOM-BC (otherwise just tax_level used)
  if (is.null(tax_agg_level)) {
    ancombc_results_df$ASVid <- ancombc_results_df$taxon

    tax_df <- as(tax_table(Grouped_phyloseq), "data.frame")
    if (!taxa_are_rows(Grouped_phyloseq)) {
      tax_df <- t(tax_df)
    }
    tax_df$ASVid <- rownames(tax_df)

    ancombc_results_df <- left_join(ancombc_results_df, tax_df, by = "ASVid")

    ancombc_results_df$Label <- ifelse(is.na(ancombc_results_df[[tax_label_level]]) | 
                                         ancombc_results_df[[tax_label_level]] == "",
                                       paste0("UC (", ancombc_results_df$ASVid, ")"),
                                       paste0(ancombc_results_df[[tax_label_level]], 
                                              "(", ancombc_results_df$ASVid, ")"))
  } else {
    ancombc_results_df$Label <- ancombc_results_df$taxon
  }

  # STEP 6: Export significant results
  sig_ancombc_results_df <- subset(ancombc_results_df,
                                   Significance == "Sig")

  # TODO: shift this responsibility to snakemake
  if (!dir.exists(paste0(args$out, "/ANCOMBC/", comp_file_label))) {
    dir.create(paste0(args$out, "/ANCOMBC/", comp_file_label),
               recursive = TRUE)
  }
  
  write.table(
    sig_ancombc_results_df,
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

#--------------------------------------------
# Volcano Plotting
#--------------------------------------------
DA_volcano_plotting <- function(DA_results_df, 
                                GroupingType, 
                                alpha, lfc_cutoff, 
                                DA_method) {
  
  #Subset to significant ASVs
  sig_DA_results_df <- subset(DA_results_df,
                              Significance == "Sig" & !is.na(padj))
 
  pos_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange > 0 & abs(log2FoldChange) > lfc_cutoff)
  neg_sig_DA_results_df <- subset(sig_DA_results_df,
                                  log2FoldChange < 0 & abs(log2FoldChange) > lfc_cutoff)



  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("CellControl", "NegativeControl", "AllControl")) {
    comp_file_label <- paste0(GroupingType, "vsPatientSample")
    plot_title_label <- paste(GroupingType, "vs Patient Samples")
  } else if (GroupingType == "PatientSample"){
    comp_file_label <- "NormalvsTumor"
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
                      aes(x = log2FoldChange, y = -log10(padj), label = Label),
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
                      aes(x = log2FoldChange, y = -log10(padj), label = Label),
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
DA_heatmap_plotting <- function(Grouped_phyloseq,
                                DA_results_df,
                                GroupingType,
                                tax_agg_level,
                                alpha,
                                lfc_cutoff,
                                pseudocount = 1) {


  #Create plot title and file labels depending on comparison
  if (GroupingType %in% c("CellControl", "NegativeControl", "AllControl")) {
    sample_type_order <- c(GroupingType, "PatientSample")
    comp_file_label <- paste0(GroupingType, "vsPatientSample")
    plot_title_label <- paste(GroupingType, "vs Patient Samples")
    group_var = "CtrlStatus"
  } else if (GroupingType == "PatientSample") {
    sample_type_order <- c("NormalTissue", "Tumor")
    comp_file_label <- "NormalvsTumor"
    plot_title_label <- "Normal Tissue vs Tumor"
    group_var = "SampleType"
  }
  plot_title_label <- paste0(plot_title_label,
                             " (alpha = ", alpha,
                             ", |LFC| > ", lfc_cutoff, ")")

  #Subset to significant ASVs
  sig_DA_results_df <- DA_results_df[DA_results_df$Significance == "Sig", , drop = FALSE]
  if (nrow(sig_DA_results_df) == 0) {
    warning("No siginifcant differnetially abundant taxa found, skipping heatmap plotting.")
    return(NULL)
  }
  sig_DA_taxa <- unique(as.character(sig_DA_results_df$taxon))
  
  # STEP ??: Get OTU table glommed at the specificed taxa aggregation level
  if (!(is.null(tax_agg_level))) {
    Glommed_phyloseq <- tax_glom(Grouped_phyloseq,
                                 taxrank = tax_agg_level, 
                                 NArm = FALSE)

  } else {
    Glommed_phyloseq <- Grouped_phyloseq
  }
  OTU_mat <- t(as(otu_table(Glommed_phyloseq), "matrix"))
  genus_names <- as.character(tax_table(Glommed_phyloseq)[, "Genus"])
  rownames(OTU_mat) <- make.unique(genus_names)

  # STEP ??: Subset OTU table to significant 
  OTU_mat <- OTU_mat[sig_DA_taxa, , drop = FALSE]
  
  # STEP ??: Normalize pseudocount values to column sums with log transfrom
  libsizes <- sample_sums(Glommed_phyloseq)
  relative_OTU_mat <- sweep(OTU_mat, 2, libsizes, FUN = "/")
  log_relative_OTU_mat <- log2(relative_OTU_mat * 1e6 + pseudocount)


  # STEP ??: Get sample metadata 
  META_df <- as(sample_data(Glommed_phyloseq), "data.frame")
  META_df$SampleName <- rownames(META_df)
  META_df$SampleType <- factor(META_df$SampleType, levels = sample_type_order)
  META_df$PatientID <- factor(META_df$PatientID)


  # STEP ??: Order Samples by SampleType followed by PatientID
  META_df <- META_df[order(META_df$SampleType,
                           META_df$PatientID), , drop = FALSE]


  
  # STEP ??: Get row labels
  da_label <- if_else(
    sig_DA_results_df$log2FoldChange >= 0,
    paste0("DA Up (LFC = ", sprintf("%.2f", sig_DA_results_df$log2FoldChange), ")"),
    paste0("DA Down (LFC = ", sprintf("%.2f", sig_DA_results_df$log2FoldChange), ")")
  )


  struc0_label <- case_when(
    is.na(sig_DA_results_df$struc0) ~ "",
    sig_DA_results_df$struc0 == sample_type_order[1] ~ paste0("Structural Zero: ", sample_type_order[1]),
    TRUE ~ paste0("Structural Zero: ", sample_type_order[2])
  )


  row_label <- ifelse(struc0_label != "", struc0_label, da_label)


  row_class_levels <- c(
    paste0("Struc0 ", sample_type_order[1]),
    paste0("Struc0 ", sample_type_order[2]),
    paste0("DA ", sample_type_order[1]),
    paste0("DA ", sample_type_order[2])
  )
  row_class <- ifelse(struc0_label != "",
                      ifelse(sig_DA_results_df$struc0 == sample_type_order[1], 
                             row_class_levels[1], 
                             row_class_levels[2]),
                      ifelse(sig_DA_results_df$log2FoldChange < 0, 
                             row_class_levels[3], 
                             row_class_levels[4]))

  # STEP ??: Ordering rows by DA
  row_order <- order(row_class, -abs(sig_DA_results_df$log2FoldChange), rownames(log_relative_OTU_mat))
  log_relative_OTU_mat <- log_relative_OTU_mat[row_order, , drop = FALSE]
  sig_DA_results_df <- sig_DA_results_df[row_order, , drop = FALSE]
  row_label <- row_label[row_order]
  row_class <- row_class[row_order]

  
  # STEP ??: Setting Colors
  # TODO: adjust to be able to handle more than two sample types
  sample_type_colors <- setNames(c("orchid", 
                                   "darkorange"), 
                                 sample_type_order)
  row_class_colors <- setNames(c("steelblue4",
                                 "firebrick4",
                                 "steelblue", 
                                 "firebrick"),
                               row_class_levels)

  # STEP ??: Heatmap Cell Coloring
  quantile_stat <- quantile(log_relative_OTU_mat,
                            probs = c(0.05, 0.5, 0.95),
                            na.rm = TRUE,
                            names = FALSE)
  if (diff(range(quantile_stat)) == 0) {
    quantile_stat <- c(min(log_relative_OTU_mat, na.rm = TRUE) - 1, 0, max(log_relative_OTU_mat) + 1)
  }

  color_function <- colorRamp2(quantile_stat,
                               c("steelblue", "white", "firebrick"))


  top_heatmap_anno <- ComplexHeatmap::HeatmapAnnotation(
    SampleType = META_df$SampleType,
    col = list(SampleType = sample_type_colors),
    annotation_name_gp = grid::gpar(fontsize = 9)
  )

  left_heatmap_anno <- ComplexHeatmap::HeatmapAnnotation(
    Class = row_class,
    Status = row_label,
    col = list(Class = row_class_colors),
    annotation_name_gp = grid::gpar(fontsize = 9),
    which = "row"
  )

  htmap <- ComplexHeatmap::Heatmap(
    log_relative_OTU_mat,
    name = "log2 Reads per Million",
    col = color_function,
    na_col = "grey90",
    top_annotation = top_heatmap_anno,
    left_annotation = left_heatmap_anno,
    column_split = META_df$SampleType,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_names = TRUE,
    column_names_rot = 90,
    row_names_side = "left",
    column_title = plot_title_label,
    use_raster = TRUE
  )

  png(
    file = paste0(args$out, "/ANCOMBC/", comp_file_label, "/", 
                  args$trialID, "_", comp_file_label, "_", "ANCOMBCHeatmap.png"),
    width = 10,
    height = 8,
    units = "in",
    res = 300
  )

  ComplexHeatmap::draw(htmap,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
  
  dev.off()
}




#-----------------------------------------
# Differential Analysis Wrapper
#-----------------------------------------
DAxPlottingWrapper <- function(Raw_phyloseq,
                               GroupingType, 
                               tax_agg_level = NULL, tax_label_level = "Genus",
                               alpha = 0.01, lfc_cutoff = 1) {

  Grouped_phyloseq <- SampleGrouping(Ungrouped_phyloseq = Raw_phyloseq, 
                                     GroupingType = GroupingType)
      
  DA_results_df <- ANCOMBC_DA(Grouped_phyloseq = Grouped_phyloseq,
                              GroupingType = GroupingType,
                              tax_agg_level = tax_agg_level,
                              tax_label_level = tax_label_level,
                              alpha = alpha,
                              lfc_cutoff = lfc_cutoff)
    
  DA_heatmap_plotting(Grouped_phyloseq = Grouped_phyloseq,
                      DA_results_df = DA_results_df,
                      GroupingType = GroupingType,
                      tax_agg_level = tax_agg_level,
                      alpha = alpha, lfc_cutoff = lfc_cutoff)

  DA_volcano_plotting(DA_results_df = DA_results_df,
                      GroupingType = GroupingType,
                      alpha = alpha,
                      lfc_cutoff = lfc_cutoff,
                      DA_method = "ANCOMBC")
  

}


#----------------------------------
# Differential Abundance Execution
#----------------------------------
all_DA_analysis <- function(DA_method, tax_agg_level = NULL, tax_label_level = "Genus") {

  GenusAbundance_tableXplot(norm_kraken_phyloseq)
  
  CtrlTypes <- c(
    "AllControl",
    "CellControl",
    "NegativeControl"
  )

  for (CtrlType in CtrlTypes) {
    DAxPlottingWrapper(Raw_phyloseq = raw_kraken_phyloseq,
                       tax_agg_level = tax_agg_level,
                       tax_label_level = tax_label_level,
                       GroupingType = CtrlType)
  }

  DAxPlottingWrapper(Raw_phyloseq = raw_kraken_phyloseq,
                     tax_agg_level = tax_agg_level,
                     tax_label_level = tax_label_level,
                     GroupingType = "PatientSample")
}

all_DA_analysis(DA_method = "ANCOMBC",
                tax_agg_level = "Genus")
#-------------------------------
# R Session Cataloging
#-------------------------------
# allows for more accessible downstream exploratory data analysis
save.image(file = paste0(args$out, "/", args$trialID, "_PhyloSeqSession.RData"))


# #--------------------------------------------
# # Violin Plotting
# #--------------------------------------------
# DA_violin_plotting <- function (Norm_phyloseq, DA_results_df, CtrlType, DA_method) {
#   # DA_results_df: data.frame of DA results (rownames = ASV ids)
#   # PStoCtrl_phyloseq: phyloseq object used for DA (contains same taxa)
  
#   # STEP 1: Subset DA-results to significant ASVs (extract significant ASVids)
#   sig_DA_results_df <- subset(DA_results_df,
#                               Significance == "Sig" & !is.na(padj))
#   if (nrow(sig_DA_results_df) == 0) {
#     message("No plotting rows (plot_df is empty) for CtrlType=", CtrlType, " — skipping violin output.")
#     return(invisible(NULL))
#   }
#   sig_ASVids <- rownames(sig_DA_results_df)
  
#   # STEP 2: Group Normalized Phyloseq by Control Type
#   norm_PStoCtrl_phyloseq <- SampleGrouping(Ungrouped_phyloseq = Norm_phyloseq,
#                                            CtrlType = CtrlType)

#   # STEP 3: Get Normalized Read Counts from Normalized Phyloseq (subset to significant ASVs)
#   otu_df <- as.data.frame(otu_table(norm_PStoCtrl_phyloseq)) # have to first extract as data frame
#   otu_mat <- as.matrix(otu_df) # otu_table() does not accept immediate matrix transformation
  
#   if (taxa_are_rows(norm_PStoCtrl_phyloseq)) {
#     otu_mat <- t(otu_mat)  # ASV's are columns
#   }

#   common_ASVids <- intersect(sig_ASVids, as.character((colnames(otu_mat))))
#   missing_ASVids <- setdiff(sig_ASVids, common_ASVids)
#   if (length(missing_ASVids) > 0) {
#     message("Warning: ", length(missing_ASVids),
#             " significant ASVs not found in normalized OTU matrix and will be skipped. Examples: ",
#             paste(head(missing_ASVids, 10), collapse = ", "))
#   }
  
#   sig_otu_mat <- otu_mat[, sig_ASVids, drop = FALSE]   # samples x sig_taxa
  
#   # STEP 4: Extract Sample Metadata from Normalized Phyloseq Object
#   sample_data_df <- as(sample_data(norm_PStoCtrl_phyloseq), "data.frame")
#   sample_data_df$SampleID <- rownames(sample_data_df)
  
#   # STEP 5: Join Normalized Read Counts, Sample Metadata, and Taxa Classification Labels (from Significant DA Results)
#   # pivot OTU to long form and join to sample data and significant DA results
#   # Resulting Data Frame: one row per read count value, Sample Meta Data and Taxa Classification duplicated as appropriate
#   plot_df <- as.data.frame(sig_otu_mat, stringsAsFactors = FALSE) %>%
#     mutate(SampleID = rownames(.)) %>%
#     pivot_longer(cols = -SampleID, names_to = "ASVid", values_to = "Count") %>%
#     left_join(sample_data_df %>% select(SampleID, CtrlStatus, SampleType, Replicate), by = "SampleID") %>%
#     left_join(sig_DA_results_df %>% select(ASVid, Genus, Label), by = "ASVid")
  
#   # STEP 6: Format dataframe for plotting
#   # ensure CtrlStatus is factor
#   plot_df$CtrlStatus <- factor(plot_df$CtrlStatus, levels = unique(plot_df$CtrlStatus))
  
#   # order ASVs by max mean count (faceting ends up in interpretable order)
#   asv_order <- plot_df %>%
#     group_by(ASVid) %>%
#     summarize(meanCount = mean(Count, na.rm = TRUE)) %>% 
#     arrange(desc(meanCount)) %>% 
#     pull(ASVid)
#   plot_df$ASVid <- factor(plot_df$ASVid, levels = asv_order)
  
#   # In case where there a many significant ASVids, view only top N:
#   top_n <- 10
#   if (length(asv_order) > top_n) {
#     chosen <- asv_order[1:top_n]
#     plot_df <- filter(plot_df, ASVid %in% chosen)
#   }
  
#   # Normalization Lab
#   norm_label_map <- c(
#     noNorm = "Un-normalized",
#     log2 = "Log2 Transformed",
#     rawTSS = "per Raw Read",
#     hostTSS = "per Host Read",
#     rawTSSlog2 = "log2(per Raw Reads)",
#     hostTSSlog2 = "log2(per Host Reads)",
#     rawTSSpM = "per 10^6 Raw Reads",
#     hostTSSpM = "per 10^6 Host Reads",
#     rawTSSpMlog2 = "log2(per 10^6 Raw Read)",
#     hostTSSpMlog2 = "log2(per 10^6 Host Reads)"              
#   )
#   norm_label <- norm_label_map[[args$norm_method]]
  
#   # Building plot
#   p_violin <- ggplot(plot_df, aes(x = CtrlStatus, y = Count, fill = CtrlStatus)) +
#     geom_violin(trim = FALSE, alpha = 0.6, width = 0.9) +
#     stat_summary(fun = mean, geom = "point", color = "black", size = 2) +
#     facet_wrap(~ Label, scales = "free_y", ncol = 2) +
#     theme_bw() +
#     theme(strip.text = element_text(size = 8),
#           axis.text.x = element_text(angle = 0, vjust = 0.5),
#           axis.title.x = element_blank()) +
#     labs(title = paste0("Per-sample counts for significant ASVs (", length(sig_ASVids), " taxa)"),
#          y = paste0("Count (", norm_label, ")")) +
#     guides(fill = guide_legend(title = "Group"))
  
#   # save
#   ggsave(
#     filename = paste0(args$out, '/', DA_method, '/', "PSto", CtrlType, '/', 
#                       args$trialID, "_PSto", CtrlType, "_", DA_method,"Violin.png"),
#     plot = p_violin,
#     width = 16, height = 10, units = "in", dpi = 300
#   )
# }