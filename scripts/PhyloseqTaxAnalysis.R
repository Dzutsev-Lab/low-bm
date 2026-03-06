library(phyloseq)
library(Biostrings)
library(ggplot2)
library(argparse)
library(dplyr)
library(taxonomizr)
library(DESeq2)

# Keep ggplot from producing Rplots.pdf
if(!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--mothur-file", type="character", help="taxonomy classification file path")
parser$add_argument("--kraken-file", type="character", help="kraken taxonomy classification file path")
parser$add_argument("--dump-dir", type="character", help="directory containing nodes.dmp and names.dmp to be used in database construction")
parser$add_argument("--norm-seq-table", type="character", help="normalized seq table tsv file path")
parser$add_argument("--raw-seq-table", type="character", help="un-normalized seq table tsv file path")
parser$add_argument("--bacterial-names", type="character", help="names of bacterial ASVs")
parser$add_argument("--trialID", type="character", help="ID number for the trial")
parser$add_argument("--add-unclassified-prefix", help="Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action="store_true",
                    default=FALSE)
parser$add_argument("--out", type="character", help="directory to store output abundance plots")


args <- parser$parse_args()



#-----------------------------------------
# Normalized Sequence Table Construction
#-----------------------------------------
#read in all ASV sequence table
norm_seq_table <- read.delim(args$norm_seq_table, header = TRUE, row.names = 1)
'Sequence Table'
str(norm_seq_table)

# read in positive bacterial ASV IDs
bacterial_IDs <- readLines(args$bacterial_names) # split IDs by line
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines
'Bacterial ASV IDs'
str(bacterial_IDs)

# select columns from sequence table to positive bacterial ASVs
norm_seq_table <- norm_seq_table[, bacterial_IDs]


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

# add taxonomic level prefixes to kraken taxonomy matrix for consistancy with mothur taxonomy
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
str(kraken_tax_matrix)

#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
kraken_phyloseq <- phyloseq(otu_table(norm_seq_table, taxa_are_rows = FALSE),
                                       sample_data(sample_meta_data_df),
                                       tax_table(kraken_tax_matrix))
str(kraken_phyloseq)

#------------------------------------
# Unclassified Taxa Handling Function
#------------------------------------
kraken_forabund_phyloseq <- if (isTRUE(args$add_unclassified_prefix)) {
  org_tax_matrix <- as(tax_table(kraken_phyloseq), "matrix")
  filled_tax_matrix <- org_tax_matrix
  dim(org_tax_matrix)
  dim(filled_tax_matrix)
  is_unassigned <- function(x) is.na(x) || x == "" || x == "Unclassified"
  
  for (i in 1:nrow(filled_tax_matrix)) {
    i
    row <- filled_tax_matrix[i, , drop = TRUE]
    assigned_index <- which(!vapply(row, is_unassigned, logical(1)))
    if (length(assigned_index) == 0) next
    lowest_assigned_index <- max(assigned_index)
    lowest_assigned_value <- row[lowest_assigned_index]
    fill_value <- paste0("Unclassified_", lowest_assigned_value)
    if (lowest_assigned_index < ncol(filled_tax_matrix)) {
      for (j in (lowest_assigned_index + 1):ncol(filled_tax_matrix)) {
        j
        if (is_unassigned(row[j])) {
          filled_tax_matrix[i, j] <- fill_value
        }
      }
    }
  }
  dim(org_tax_matrix)
  dim(filled_tax_matrix)
  
  kraken_filled_phyloseq <- kraken_phyloseq
  tax_table(kraken_filled_phyloseq) <- filled_tax_matrix
  kraken_filled_phyloseq
} else {
  kraken_phyloseq
}

kraken_genusGlom_phyloseq <- tax_glom(kraken_forabund_phyloseq, 
                                       taxrank = "Genus",
                                       NArm = TRUE
                                       )

kraken_top10genus_phyloseq <- prune_taxa(
  names(sort(taxa_sums(kraken_genusGlom_phyloseq), decreasing = TRUE))[1:10],
  kraken_genusGlom_phyloseq
)


#--------------------------------------
# Kraken Genus Table Construction
#--------------------------------------
# constructing genus name matrix from tax table of genus glommed phyloseq object
kraken_genus_name_matrix <- as(tax_table(kraken_genusGlom_phyloseq), "matrix")[, "Genus"]

# setting taxa names in phyloseq object to genus level name matrix
taxa_names(kraken_genusGlom_phyloseq) <- kraken_genus_name_matrix

# Constructing genus count table from glommed OTU table
kraken_genus_matrix <- as.matrix(otu_table(kraken_genusGlom_phyloseq))
str(kraken_genus_matrix)

#ensure proper orientation of genus table
if (taxa_are_rows(kraken_genusGlom_phyloseq)) {
  kraken_genus_matrix <- t(kraken_genus_matrix)
}

# Add total row
kraken_genus_matrix <- rbind(kraken_genus_matrix, Total = colSums(kraken_genus_matrix))

write.table(
  kraken_genus_matrix,
  file = paste0(args$out, '/', args$trialID, "_kraken_genus_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA 
)

#------------------------------
# Count Plot
#------------------------------
theme_set(theme_bw())

kraken_abunXtypeXsample_plot <- plot_bar(kraken_top10genus_phyloseq, "Replicate", fill = "Genus")
kraken_abunXtypeXsample_plot +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  facet_wrap(~SampleType, scales = "free_x") +
  labs(title = "Read Counts per Sample by Sample Type",
       x = "Sample",
       y = "Read Count")
ggsave(paste0(args$out, '/', args$trialID, "_kraken_countXtypeXsample.png"), width = 14, height = 12, units = "in")

#----------------------------------------------------
# Differential Abundance Phyloseq Object Contruction
#----------------------------------------------------
#copy normalized phyloseq object
raw_kraken_phyloseq <- kraken_phyloseq

#import raw read count sequence table
raw_seq_table <- read.delim(args$raw_seq_table, 
														header = TRUE, 
														row.names = 1)
raw_seq_table <- raw_seq_table[, bacterial_IDs]

otu_table(raw_kraken_phyloseq) <- otu_table(raw_seq_table, taxa_are_rows = FALSE)

#------------------------------
# Differential Analysis
#------------------------------
PStoCtrlStatusDA <- function(CtrlType) {
  
  sample_df_orig <- as(sample_data(raw_kraken_phyloseq), "data.frame")
  
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
    
    PStoCtrl_phyloseq <- raw_kraken_phyloseq
    sample_data(PStoCtrl_phyloseq) <- sample_data(sample_df[ sample_names(PStoCtrl_phyloseq), , drop = FALSE ])
                                                  
  } else {
    sample_df <- sample_df_orig
    select_samples <- rownames(sample_df)[ sample_df$CtrlStatus %in% c("PatientSample", CtrlType) ]

    PStoCtrl_phyloseq <- prune_samples(select_samples, raw_kraken_phyloseq)
    
    # Checking if comparison is valid (i.e. we have at least one of each sample category)
    count_check <- as(sample_data(PStoCtrl_phyloseq), "data.frame")
    group_counts <- table(count_check$CtrlStatus)
    if (length(group_counts) < 2 || any(group_counts == 0)) {
      message("Skipping comparison PatientSample vs ", CtrlType,
              " because one or more groups have zero samples.")
      message("Group counts: ", paste(names(group_counts), group_counts, collapse = " | "))
      return(invisible(NULL))
    }
    
    #remove any taxa that now have zero count after filtering by sample type
    PStoCtrl_phyloseq <- prune_taxa(taxa_sums(PStoCtrl_phyloseq) > 0, PStoCtrl_phyloseq)
  }
  

  # convert phyloseq object to DESeqDataSet split by CtrlStatus label (i.e. PS, selected control type)
  PStoCtrl_DDS <- phyloseq_to_deseq2(PStoCtrl_phyloseq, ~ CtrlStatus)
  
  # estimate size factors using positive counts (addresses sparse dataset issues, avoiding log 0 errors)
  PStoCtrl_DDS <- estimateSizeFactors(PStoCtrl_DDS, type="poscounts")
  
  # run DA analysis 
  PStoCtrl_DDS <- DESeq(PStoCtrl_DDS, test="Wald", fitType="parametric")
  
  # Get DA results and turn into data.frame
  DA_results <- results(PStoCtrl_DDS, cooksCutoff = FALSE)
  DA_df <- as.data.frame(DA_results, stringAsFactor = FALSE)
  
  # Filter significant rows
  alpha = 0.01
  sig_idx <- which(!is.na(DA_df$padj) & (DA_df$padj < alpha))
  if (length(sig_idx) == 0) {
    sig_DA_results <- DA_df[0, , drop = FALSE] #empty data.frame with the same columns as DA results
  } else {
    sig_DA_results <- DA_df[sig_idx, , drop = FALSE]
    
    #retrieve tax table
    tax_tab <- as.matrix(tax_table(PStoCtrl_phyloseq))
    
    # get canonical list of significant taxa (intersection of taxa that exist in taxa table and those found to be significant in DA)
    # should be the same with perfect overlap (simple safety measure for ordering purposes)
    sig_taxa <- intersect(rownames(sig_DA_results), rownames(tax_tab))
    
    # subset and reorder for perfect match between significant DA results and taxa table (subset to siginificant taxa)
    sig_DA_results <- sig_DA_results[sig_taxa, , drop = FALSE]
    sig_tax_tab <- tax_tab[sig_taxa, , drop = FALSE]
    
    # combine to final table
    sig_DA_results <- cbind(sig_DA_results, as.data.frame(sig_tax_tab, stringsAsFactors = FALSE))
  }
  
  write.table(
    sig_DA_results,
    file = paste0(args$out, '/', args$trialID, "_PSto", CtrlType,"_DA_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA 
  )
  
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {}
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  x = tapply(sig_DA_results$log2FoldChange, sig_DA_results$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sig_DA_results$Phylum = factor(as.character(sig_DA_results$Phylum), levels=names(x))
  ggplot(sig_DA_results, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + theme(axis.test.x = element_text(angle= -90, hjust = 0, vjust=0.5))
  ggsave(paste0(args$out, '/', args$trialID, "_PSto", CtrlType,"_DA_results.png"), width = 14, height = 12, units = "in")
}

PStoCtrlStatusDA("ExpCtrl")
PStoCtrlStatusDA("NegCtrl")
PStoCtrlStatusDA("AllCtrl")

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