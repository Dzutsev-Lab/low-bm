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

# Keep ggplot from producing Rplots.pdf
if(!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--mothur-file", type="character", help="taxonomy classification file path")
parser$add_argument("--kraken-file", type="character", help="kraken taxonomy classification file path")
parser$add_argument("--dump-dir", type="character", help="directory containing nodes.dmp and names.dmp to be used in database construction")
parser$add_argument("--norm-method", type="character", help="normalization method used to produce normalized seq table")
parser$add_argument("--norm-seq-table", type="character", help="normalized seq table tsv file path")
parser$add_argument("--raw-seq-table", type="character", help="un-normalized seq table tsv file path")
parser$add_argument("--bacterial-names", type="character", help="names of bacterial ASVs")
parser$add_argument("--trialID", type="character", help="ID number for the trial")
parser$add_argument("--add-unclassified-prefix", help="Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action="store_true",
                    default=FALSE)
parser$add_argument("--read-count-file", type="character", help="combined read count file path")
parser$add_argument("--out", type="character", help="directory to store output abundance plots")


args <- parser$parse_args()


#----------------------------------------
# Positive Bacterial ASVID Extraction
#----------------------------------------
bacterial_IDs <- readLines(args$bacterial_names) # split IDs by line
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines

#-----------------------------------
# Read Count Extraction
#-----------------------------------
read_counts_df <- read.delim(args$read_count_file, header = TRUE, row.names = 1)

#-----------------------------------------
# Sequence Table Construction
#-----------------------------------------
# NORMALIZED
norm_seq_table <- read.delim(args$norm_seq_table, 
                             header = TRUE, 
                             row.names = 1)
norm_seq_table <- norm_seq_table[, bacterial_IDs]

# UN-NORMALIZED
raw_seq_table <- read.delim(args$raw_seq_table, 
                            header = TRUE, 
                            row.names = 1)
raw_seq_table <- raw_seq_table[, bacterial_IDs]


#--------------------------------
# Sample Data Table Construction
#--------------------------------
sample_names <- rownames(norm_seq_table)
sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)

# Sample Type
sample_type <- sub("\\d+[A-Za-z]*$", "", sample_info)

# Technical Rep (Will also denote patient ID for patient samples)
tech_rep <- sub(".*?(\\d+[A-Za-z]*)$", "\\1", sample_info)

sample_meta_data_df <- data.frame(
  SampleType = factor(sample_type),
  Replicate = factor(tech_rep),
  row.names = rownames(norm_seq_table),
  stringsAsFactors = FALSE
)

sample_meta_data_df <- sample_meta_data_df %>% 
  mutate(
    SampleType = case_when(
      (is.na(SampleType) | SampleType == "") & (grepl("NT$", Replicate) | grepl("N$", Replicate)) ~ "NormalTissue", #check before tumor otherwise all would be labeled tumor (ending with T)
      (is.na(SampleType) | SampleType == "") & grepl("T$", Replicate) ~ "Tumor",
      TRUE ~ SampleType
    )
  )

sample_meta_data_df <- sample_meta_data_df %>%
  mutate(
    CtrlStatus = case_when(
      SampleType %in% c("Tumor", "NormalTissue") ~ "PatientSample",
      SampleType == "expcontrol" ~ "ExpCtrl",
      SampleType == "NEGATIVECONTROL" ~ "NegCtrl",
      TRUE ~ as.character(SampleType)
    ),
    CtrlStatus = factor(CtrlStatus)
  )

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

desiredTaxa <- c("superkingdom","phylum","class","order","family","genus","species")

kraken_tax_df <- getTaxonomy(ids = kraken_info$taxid, 
                      sqlFile = sql_db, 
                      desiredTaxa = desiredTaxa)
kraken_tax_df <- as.data.frame(kraken_tax_df, stringsAsFactors = FALSE)
kraken_tax_df$ASVid <- kraken_info$ASVid

kraken_tax_df <- rename( kraken_tax_df, 
                  Domain  = superkingdom,
                  Phylum  = phylum,
                  Class   = class,
                  Order   = order,
                  Family  = family,
                  Genus   = genus,
                  Species = species)


kraken_tax_matrix <- as.matrix(kraken_tax_df[, c("Domain","Phylum","Class","Order","Family","Genus","Species")])
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
        tax_matrix[j, i] <- paste0(tax_prefixes[colnames(tax_matrix)[i]], "__", tax_matrix[j, i])
      } else {
        tax_matrix[j, i] <- NA
      }
    }
  }

  return(tax_matrix)
}
kraken_tax_matrix <- taxa_level_prefix_addition(kraken_tax_matrix)

# Handling Unclassified Taxa
unclassified_label_progigation <- function(tax_matrix) {
  
  is_unassigned <- function(x) is.na(x) || x == "" || x == "Unclassified"
  
  for (i in 1:nrow(tax_matrix)) {
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
  return(tax_matrix)
}

kraken_tax_matrix <- if (isTRUE(args$add_unclassified_prefix)) {
  unclassified_label_progigation(kraken_tax_matrix)
}


#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
norm_kraken_phyloseq <- phyloseq(otu_table(norm_seq_table, taxa_are_rows = FALSE),
                                       sample_data(sample_meta_data_df),
                                       tax_table(kraken_tax_matrix))

raw_kraken_phyloseq <- phyloseq(otu_table(raw_seq_table, taxa_are_rows = FALSE),
                                sample_data(sample_meta_data_df),
                                tax_table(kraken_tax_matrix))

#--------------------------------------
# Abundance Table and Plot
#--------------------------------------
GenusAbundance_tableXplot <- function(norm_phyloseq) {
  
  # STEP 1: Genus Glomming
  genusGlom_phyloseq <- tax_glom(norm_phyloseq, 
                                 taxrank = "Genus",
                                 NArm = TRUE)
  
  top10genus_phyloseq <- prune_taxa(
    names(sort(taxa_sums(genusGlom_phyloseq), decreasing = TRUE))[1:10],
    genusGlom_phyloseq
  )
  
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
    file = paste0(args$out, '/', args$trialID, "_kraken_genus_table.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  # STEP 3: Count Plot
  theme_set(theme_bw())
  
  abunXtypeXsample_plot <- plot_bar(top10genus_phyloseq, "Replicate", fill = "Genus")
  abunXtypeXsample_plot +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 6)
    ) +
    facet_wrap(~SampleType, scales = "free_x") +
    labs(title = "Read Counts per Sample by Sample Type",
         x = "Sample",
         y = "Read Count")
  ggsave(paste0(args$out, '/', args$trialID, "_kraken_countXtypeXsample.png"), width = 14, height = 12, units = "in")
}

GenusAbundance_tableXplot(norm_kraken_phyloseq)



#---------------------------------------
# Grouping Samples by Control Type
#---------------------------------------
CtrlGrouping <- function(Ungrouped_phyloseq, CtrlType) {
  
  sample_df_orig <- as(sample_data(Ungrouped_phyloseq), "data.frame")
  
  if (CtrlType == "AllCtrl") {
    sample_df <- sample_df_orig %>%
      mutate(
        CtrlStatus = if_else(
          CtrlStatus %in% c("NegCtrl", "ExpCtrl"),
          "AllCtrl",
          as.character(CtrlStatus)
        )
      )
    #ensure sample_df is correctly formatted with rownames and factored CtrlStatus
    rownames(sample_df) <- rownames(sample_df_orig)
    sample_df$CtrlStatus <- factor(sample_df$CtrlStatus)
    
    CtrlGrouped_phyloseq <- Ungrouped_phyloseq
    sample_data(CtrlGrouped_phyloseq) <- sample_data(sample_df[ sample_names(CtrlGrouped_phyloseq), , drop = FALSE ])
    
  } else {
    sample_df <- sample_df_orig
    select_samples <- rownames(sample_df)[ sample_df$CtrlStatus %in% c("PatientSample", CtrlType) ]
    
    CtrlGrouped_phyloseq <- prune_samples(select_samples, Ungrouped_phyloseq)
    
    # Checking if comparison is valid (i.e. we have at least one of each sample category)
    count_check <- as(sample_data(CtrlGrouped_phyloseq), "data.frame")
    group_counts <- table(count_check$CtrlStatus)
    if (length(group_counts) < 2 || any(group_counts == 0)) {
      message("Skipping comparison PatientSample vs ", CtrlType,
              " because one or more groups have zero samples.")
      message("Group counts: ", paste(names(group_counts), group_counts, collapse = " | "))
      return(invisible(NULL))
    }
    
    #remove any taxa that now have zero count after filtering by sample type
    CtrlGrouped_phyloseq <- prune_taxa(taxa_sums(CtrlGrouped_phyloseq) > 0, CtrlGrouped_phyloseq)
  }
  
  return(CtrlGrouped_phyloseq)
}


#----------------------------------------------
# ANCOMBC Differential Abundance Anlysis
#----------------------------------------------
ANCOMBC_PStoCtrl_DA <- function(Raw_phyloseq, CtrlType, alpha, lfc_cutoff) {
  
  #STEP 1: Group by Control Status
  Grouped_phyloseq <- CtrlGrouping(Ungrouped_phyloseq = Raw_phyloseq, 
                                   CtrlType = CtrlType)
  
  #STEP 2: Run ANCOM-BC Analysis
  # TODO: parameterize tax_level
  ancombc_output <- ancombc(data = Grouped_phyloseq, tax_level = "Genus",
                            formula = "CtrlStatus",
                            p_adj_method = "holm",
                            group = "CtrlStatus",
                            alpha = alpha)
  
  # STEP 3: Construct formatted ANCOM-BC results data frame
  ancombc_select_cols <- c("taxon", "CtrlStatusPatientSample")
  ancombc_lfc_df <- ancombc_output[["res"]][["lfc"]][, ancombc_select_cols]
  colnames(ancombc_lfc_df) <- c("Genus", "log2FoldChange")
  
  ancombc_padj_df <- ancombc_output[["res"]][["q_val"]][, ancombc_select_cols]
  colnames(ancombc_padj_df) <- c("Genus", "padj")
  
  ancombc_results_df <- left_join(ancombc_lfc_df, ancombc_padj_df, by = "Genus")
  
  # STEP 4: Add significant and plotting labels
  ancombc_results_df$Significance <- "NotSig"
  ancombc_results_df$Significance[
    ancombc_results_df$padj < alpha &
    abs(ancombc_results_df$log2FoldChange) > lfc_cutoff
  ] <- "Sig"
  
  ancombc_results_df$Label <- ancombc_results_df$Genus
  
  # STEP 5: Export significant results
  sig_ancombc_results_df <- subset(ancombc_results_df,
                                   Significance == "Sig" & !is.na(padj))
  
  # TODO: shift this responsibility to snakemake
  if (!dir.exists(paste0(args$out, '/ANCOMBC/', "PSto", CtrlType))) {
    dir.create(paste0(args$out, '/ANCOMBC/', "PSto", CtrlType),
               recursive = TRUE)
  }
  
  write.table(
    sig_ancombc_results_df,
    file = paste0(args$out, '/ANCOMBC/', "PSto", CtrlType,'/', 
                  args$trialID, "_PSto", CtrlType,"_ANCOMBCResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  return(ancombc_results_df)
  
}

#----------------------------------------------
# DESeq2 Differential Abundance Analysis
#----------------------------------------------
DESeq_PStoCtrl_DA <- function(Raw_phyloseq, CtrlType, alpha, lfc_cutoff) {
  
  # STEP 1: Group by Control Status and convert phyloseq object to DESeqDataSet
  Grouped_phyloseq <- CtrlGrouping(Ungrouped_phyloseq = Raw_phyloseq, 
                                   CtrlType = CtrlType)
  PStoCtrl_DDS <- phyloseq_to_deseq2(Grouped_phyloseq, ~ CtrlStatus)
  
  # STEP 2: Prepare normMatrix with host read normalization factors
  sample_order <- sample_names(Grouped_phyloseq)
  
  chimera_filtered_vec <- as.numeric(read_counts_df[sample_order, "chimera.filtered"])
  host_unmapped_vec <- as.numeric(read_counts_df[sample_order, "HostUnmapped_reads"])
  host_reads_per_sample <- chimera_filtered_vec - host_unmapped_vec
  
  names(host_reads_per_sample) <- sample_order
  host_reads_per_sample[host_reads_per_sample <= 0] <- 1 # replacing any zeros with small number to avoid zero divisor problem
  
  n_features <- ntaxa(Grouped_phyloseq) # number of taxa (ASVids) will dictate dimensions of normMatrix
  normMatrix <- matrix(rep(host_reads_per_sample, each = n_features),
                       nrow= n_features, ncol = length(host_reads_per_sample),
                       dimnames = list(taxa_names(Grouped_phyloseq),
                                       names(host_reads_per_sample)))
  row_geom_mean <- apply(normMatrix, 1, function(row) exp(mean(log(row))))
  normMatrix <- normMatrix / row_geom_mean
  
  #check that normMatrix rows match DDS counts table
  dds_colnames <- colnames(counts(PStoCtrl_DDS))
  if (!identical(colnames(normMatrix), dds_colnames)) {
    normMatrix <- normMatrix[, dds_colnames, drop=FALSE]
  }
  # Final explicit matrix formatting (all should be redundant with previous steps)
  normMatrix <- apply(normMatrix, 2, as.numeric)
  rownames(normMatrix) <- taxa_names(Grouped_phyloseq)
  colnames(normMatrix) <- dds_colnames
  
  # STEP 3: estimate size factors 
  # using Positive Counts Method (addresses sparse data set issues, avoiding log 0 errors) 
  # and Host Reac Count Normalization Matrix
  PStoCtrl_DDS <- estimateSizeFactors(PStoCtrl_DDS, type="poscounts", normMatrix = normMatrix)
  
  # STEP 4: Run DA Analysis
  PStoCtrl_DDS <- DESeq(PStoCtrl_DDS, test="Wald", fitType="parametric")
  DA_results_df <- as.data.frame(results(PStoCtrl_DDS, cooksCutoff = FALSE), stringAsFactor = FALSE)
  DA_results_df$ASVid <- rownames(DA_results_df)
  
  # STEP 5: Attach Taxa Classifications
  tax_df <- as.data.frame(tax_table(Grouped_phyloseq))
  tax_df$ASVid <- rownames(tax_df)
  
  # combine to final table
  DA_results_df <- left_join(DA_results_df, tax_df, by="ASVid")
  
  # STEP 6: Add additional label fields 
  # Taxa label with ASVid (with just ASVid fallback when no classification available)
  DA_results_df$Label <- ifelse(is.na(DA_results_df$Genus) | DA_results_df$Genus == "",
                                DA_results_df$ASVid,
                                paste0(DA_results_df$Genus, '(', DA_results_df$ASVid, ')'))
  
  # Significance label by alpha and lfc cutoff
  DA_results_df$Significance <- "NotSig"
  DA_results_df$Significance[
    DA_results_df$padj < alpha &
    abs(DA_results_df$log2FoldChange) > lfc_cutoff
  ] <- "Sig"
  
  rownames(DA_results_df) <- DA_results_df$ASVid
  
  # STEP 7: Export Significant Results
  sig_idx <- which(DA_results_df$Significance == "Sig")
  if (length(sig_idx) == 0) {
    sig_DA_results_df <- DA_results_df[0, , drop = FALSE] #empty data.frame with the same columns as DA results
  } else {
    sig_DA_results_df <- DA_results_df[sig_idx, , drop = FALSE]
  }
  
  # Create directory for specific comparison
  # TODO: shift this responsibility to snakemake
  if (!dir.exists(paste0(args$out, '/DESeq/', "PSto", CtrlType))) {
    dir.create(paste0(args$out, '/DESeq/', "PSto", CtrlType),
               recursive = TRUE)
  }
  
  write.table(
    sig_DA_results_df,
    file = paste0(args$out, '/DESeq/', "PSto", CtrlType,'/', 
                  args$trialID, "_PSto", CtrlType,"_DESeqResults.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  # STEP 8: Export DESeq Normalized Counts with per-sample sums
  DESeq_norm_counts <- counts(PStoCtrl_DDS, normalized = TRUE)
  sample_totals_vec <- colSums(DESeq_norm_counts)
  sample_totals_df <- as.data.frame(sample_totals_vec)
  write.table(sample_totals_df,
              file = paste0(args$out, '/DESeq/', "PSto", CtrlType,'/', 
                            args$trialID, "_PSto", CtrlType,"_DESeqNormCountTotals.tsv"),
              sep = "\t",
              quote = FALSE,
              col.names = NA 
  )
  
  return(DA_results_df)
}

#--------------------------------------------
# Volcano Plotting
#--------------------------------------------
DA_volcano_plotting <- function(DA_results_df, CtrlType, alpha, lfc_cutoff, DA_method) {
  
  #Subset to significant ASVs
  sig_DA_results_df <- subset(DA_results_df,
                              Significance == "Sig" & !is.na(padj))
  
  negLFC_DA_results_df <- subset(DA_results_df,
                                 (log2FoldChange < 0))
  
  # Build ggplot object
  p_volcano <- ggplot(DA_results_df,
                      aes(x = log2FoldChange,
                          y = -log10(padj),
                          color = Significance)) +
    theme_bw() +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "darkred") +
    scale_color_manual(values = c("NotSig" = "grey70",
                                  "Sig" = "deeppink4")) +
    labs(title = paste("Volcano Plot: PatientSample vs", CtrlType),
         x = "Effect size: log2(Fold Change)",
         y = "-log10(adjusted p-value)")
  
  if (nrow(sig_DA_results_df) > 0) {
    p_volcano <- p_volcano +
      geom_point(data = sig_DA_results_df,
                aes(x = log2FoldChange, y = -log10(padj)),
                color = "deeppink4",
                size = 2) +
      geom_text_repel(data = sig_DA_results_df,
                      aes(x = log2FoldChange, y = -log10(padj), label = Label),
                      size = 3.2, 
                      force = 3,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      point.padding = 0.4,
                      # segment(leader) line style
                      segment.size = 0.3,
                      segment.color = "grey40",
                      segment.alpha = 0.9,
                      min.segment.length = 0)
  }
  if (nrow(negLFC_DA_results_df) > 0) {
    p_volcano <- p_volcano +
      geom_point(data = negLFC_DA_results_df,
                 aes(x = log2FoldChange, y = -log10(padj)),
                 color = "darkolivegreen",
                 size = 2) +
      geom_text_repel(data = negLFC_DA_results_df,
                      aes(x = log2FoldChange, y = -log10(padj), label = Label),
                      color = "darkolivegreen",
                      size = 3.2, 
                      force = 3,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      point.padding = 0.4,
                      # segment(leader) line style
                      segment.size = 0.3,
                      segment.color = "grey40",
                      segment.alpha = 0.9,
                      min.segment.length = 0)
  }

  
  # Save plot
  ggsave(
    filename = paste0(args$out, '/', DA_method, '/', "PSto", CtrlType, '/', 
                      args$trialID, "_PSto", CtrlType,"_", DA_method, "Volcano.png"),
    plot = p_volcano,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
  )
}

#--------------------------------------------
# Violin Plotting
#--------------------------------------------
DA_violin_plotting <- function (Norm_phyloseq, DA_results_df, CtrlType, DA_method) {
  # DA_results_df: data.frame of DA results (rownames = ASV ids)
  # PStoCtrl_phyloseq: phyloseq object used for DA (contains same taxa)
  
  # STEP 1: Subset DA-results to significant ASVs (extract significant ASVids)
  sig_DA_results_df <- subset(DA_results_df,
                              Significance == "Sig" & !is.na(padj))
  if (nrow(sig_DA_results_df) == 0) {
    message("No plotting rows (plot_df is empty) for CtrlType=", CtrlType, " — skipping violin output.")
    return(invisible(NULL))
  }
  sig_ASVids <- rownames(sig_DA_results_df)
  
  # STEP 2: Group Normalized Phyloseq by Control Type
  norm_PStoCtrl_phyloseq <- CtrlGrouping(Ungrouped_phyloseq = Norm_phyloseq,
                                         CtrlType = CtrlType)

  # STEP 3: Get Normalized Read Counts from Normalized Phyloseq (subset to significant ASVs)
  otu_df <- as.data.frame(otu_table(norm_PStoCtrl_phyloseq)) # have to first extract as data frame
  otu_mat <- as.matrix(otu_df) # otu_table() does not accept immediate matrix transformation
  
  if (taxa_are_rows(norm_PStoCtrl_phyloseq)) {
    otu_mat <- t(otu_mat)  # ASV's are columns
  }

  common_ASVids <- intersect(sig_ASVids, as.character((colnames(otu_mat))))
  missing_ASVids <- setdiff(sig_ASVids, common_ASVids)
  if (length(missing_ASVids) > 0) {
    message("Warning: ", length(missing_ASVids),
            " significant ASVs not found in normalized OTU matrix and will be skipped. Examples: ",
            paste(head(missing_ASVids, 10), collapse = ", "))
  }
  
  sig_otu_mat <- otu_mat[, sig_ASVids, drop = FALSE]   # samples x sig_taxa
  
  # STEP 4: Extract Sample Metadata from Normalized Phyloseq Object
  sample_data_df <- as(sample_data(norm_PStoCtrl_phyloseq), "data.frame")
  sample_data_df$SampleID <- rownames(sample_data_df)
  
  # STEP 5: Join Normalized Read Counts, Sample Metadata, and Taxa Classification Labels (from Significant DA Results)
  # pivot OTU to long form and join to sample data and significant DA results
  # Resulting Data Frame: one row per read count value, Sample Meta Data and Taxa Classification duplicated as appropriate
  plot_df <- as.data.frame(sig_otu_mat, stringsAsFactors = FALSE) %>%
    mutate(SampleID = rownames(.)) %>%
    pivot_longer(cols = -SampleID, names_to = "ASVid", values_to = "Count") %>%
    left_join(sample_data_df %>% select(SampleID, CtrlStatus, SampleType, Replicate), by = "SampleID") %>%
    left_join(sig_DA_results_df %>% select(ASVid, Genus, Label), by = "ASVid")
  
  # STEP 6: Format dataframe for plotting
  # ensure CtrlStatus is factor
  plot_df$CtrlStatus <- factor(plot_df$CtrlStatus, levels = unique(plot_df$CtrlStatus))
  
  # order ASVs by max mean count (faceting ends up in interpretable order)
  asv_order <- plot_df %>%
    group_by(ASVid) %>%
    summarize(meanCount = mean(Count, na.rm = TRUE)) %>% 
    arrange(desc(meanCount)) %>% 
    pull(ASVid)
  plot_df$ASVid <- factor(plot_df$ASVid, levels = asv_order)
  
  # In case where there a many significant ASVids, view only top N:
  top_n <- 10
  if (length(asv_order) > top_n) {
    chosen <- asv_order[1:top_n]
    plot_df <- filter(plot_df, ASVid %in% chosen)
  }
  
  # Normalization Lab
  norm_label_map <- c(
    noNorm = "Un-normalized",
    log2 = "Log2 Transformed",
    rawTSS = "per Raw Read",
    hostTSS = "per Host Read",
    rawTSSlog2 = "log2(per Raw Reads)",
    hostTSSlog2 = "log2(per Host Reads)",
    rawTSSpM = "per 10^6 Raw Reads",
    hostTSSpM = "per 10^6 Host Reads",
    rawTSSpMlog2 = "log2(per 10^6 Raw Read)",
    hostTSSpMlog2 = "log2(per 10^6 Host Reads)"              
  )
  norm_label <- norm_label_map[[args$norm_method]]
  
  # Building plot
  p_violin <- ggplot(plot_df, aes(x = CtrlStatus, y = Count, fill = CtrlStatus)) +
    geom_violin(trim = FALSE, alpha = 0.6, width = 0.9) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 2) +
    facet_wrap(~ Label, scales = "free_y", ncol = 2) +
    theme_bw() +
    theme(strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 0, vjust = 0.5),
          axis.title.x = element_blank()) +
    labs(title = paste0("Per-sample counts for significant ASVs (", length(sig_ASVids), " taxa)"),
         y = paste0("Count (", norm_label, ")")) +
    guides(fill = guide_legend(title = "Group"))
  
  # save
  ggsave(
    filename = paste0(args$out, '/', DA_method, '/', "PSto", CtrlType, '/', 
                      args$trialID, "_PSto", CtrlType, "_", DA_method,"Violin.png"),
    plot = p_violin,
    width = 16, height = 10, units = "in", dpi = 300
  )
}

DAxPlottingWrapper <- function(Raw_phyloseq, Norm_phyloseq, CtrlType,
                               alpha = 0.01, lfc_cutoff = 2.5,
                               DA_method) {
  if (DA_method == "DESeq") {
    DA_results_df <- DESeq_PStoCtrl_DA(Raw_phyloseq = Raw_phyloseq,
                                             CtrlType = CtrlType,
                                             alpha = alpha,
                                             lfc_cutoff = lfc_cutoff)
  } else if (DA_method == "ANCOMBC") {
    DA_results_df <- ANCOMBC_PStoCtrl_DA(Raw_phyloseq = Raw_phyloseq,
                                         CtrlType = CtrlType,
                                         alpha = alpha,
                                         lfc_cutoff = lfc_cutoff)
  }

  DA_volcano_plotting(DA_results_df = DA_results_df,
                      CtrlType = CtrlType,
                      alpha = alpha,
                      lfc_cutoff = lfc_cutoff,
                      DA_method = DA_method)
  
  # Violin plotting currently only works with DESeq's results data frame because it contains ASVids
  # ANCOM-BC collapses at the desired taxa level so need to consider how to reintegrate ASVid
  # Alternative: make separate violin plot capability that facets by genus not asvid
  if (DA_method == "DESeq") {
    DA_violin_plotting(Norm_phyloseq = Norm_phyloseq,
                       DA_results_df = DA_results_df,
                       CtrlType = CtrlType,
                       DA_method = DA_method)
  }
}


#----------------------------------
# Differential Abundance Execution
#----------------------------------
CtrlTypes <- c(
  'AllCtrl',
  'ExpCtrl',
  'NegCtrl'
)

for (CtrlType in CtrlTypes) {
  DAxPlottingWrapper(Raw_phyloseq = raw_kraken_phyloseq,
                     Norm_phyloseq = norm_kraken_phyloseq,
                     CtrlType = CtrlType,
                     DA_method = "ANCOMBC")
}


#-------------------------------
# R Session Cataloging
#-------------------------------
# allows for more accessible downstream exploratory data analysis
save.image(file = paste0(args$out, '/', args$trialID, "_PhyloSeqSession.RData"))



#------------------------------
# Beta Diversity Plots
#------------------------------
# distance_methods <- unlist(distanceMethodList)
# # Remove distance metrics that require trees as a part of the phyloseq object
# distance_methods <- distance_methods[-(1:3)]

# # Remove user-defined method (we havent defined it)
# distance_methods <- distance_methods[-which(distance_methods=='ANY')]

# head(distance_methods)
# #construct plot list for distance panel
# plist <- vector("list", length(distance_methods))
# names(plist) = distance_methods
 # for( i in distance_methods ){
 #   dist_matrix <- distance(kraken_phyloseq, method = i)
 #   ordination_matrix <- ordinate(kraken_phyloseq, "MDS", distance = dist_matrix)
 #  
 #   #PLOTTING
 #   # clear previous plot
 #   p <- NULL
 #  
 #   p <- plot_ordination(kraken_phyloseq, ordination_matrix, color="SampleType")
 #   p <- p + ggtitle(paste("MDS using ditance method ", i, sep=""))
 #   plist[[i]] = p
 # }
# p_df = ldply(plist, function(x) x$data)
# names(p_df)[1] <- "distance"
# p = ggplot(p_df, aes(Axis.1, Axis.2, color=SampleType))
# p = p + geom_point(size=3, alpha=0.5)
# p = p + facet_wrap(~distance, scales = "free")
# p = p + ggtitle("MDS on various distance metrics for Kraken Classification")
# p

# ggsave(paste0(args$out, '/', args$trialID, "_kraken_beta.png"), width = 8, height = 6, units = "in")



# #------------------------------------
# # Mothur Taxonomy Table Construction
# #------------------------------------
# #read in taxonomy file
# mothur_file <- args$mothur_file
# mothur_tax_lines <- readLines(mothur_file) # split records by line
# mothur_tax_lines <- mothur_tax_lines[nzchar(mothur_tax_lines)] # removes empty lines
# 
# 
# # convert to simple data frame
# mothur_tax_df <- do.call(rbind, strsplit(unlist(mothur_tax_lines), "\t")) # splits ASV ID and taxonomy record
# colnames(mothur_tax_df) <- c("ASVid", "taxonomy_record")
# mothur_tax_df <- as.data.frame(mothur_tax_df, stringsAsFactors = FALSE)
# str(mothur_tax_df)
# 
# # split taxonomical levels
# mothur_split_tax_df <- strsplit(sub(";+$", "", mothur_tax_df$taxonomy_record), ";") # removes trailing and ";" and splits record field per level at ";"
# 
# maxranks <- 9 #hardset maximum number of rank assignment (down to sub-strain)
# 
# # function for removing confidence score from level entry
# clean_confidence_score <- function(x) sub("\\(.*\\)$", "", x) 
# 
# 
# # construct taxonomy matrix via matrix transposition of split taxonomy records
# mothur_tax_matrix <- t(vapply(mothur_split_tax_df, function(tax_entry) {
#   tax_entry <- sapply(tax_entry, clean_confidence_score)
#   # pad unassigned lower levels with NA
#   if (length(tax_entry) < maxranks) {
#     tax_entry <- c(tax_entry, rep(NA, maxranks - length(tax_entry)))
#   }
#   tax_entry
# }, FUN.VALUE = character(maxranks)))
# 
# rownames(mothur_tax_matrix) <- mothur_tax_df$ASVid
# colnames(mothur_tax_matrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species","Strain","Substrain")
# mothur_tax_matrix <- as.matrix(mothur_tax_matrix) #convert to matrix for ease of use with phyloseq
# 
# 
# str(mothur_tax_matrix)