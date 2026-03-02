library(phyloseq)
library(Biostrings)
library(ggplot2)
library(argparse)
library(dplyr)
library(taxonomizr)

# Keep ggplot from producing Rplots.pdf
if(!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--mothur-file", type="character", help="taxonomy classification file path")
parser$add_argument("--kraken-file", type="character", help="kraken taxonomy classification file path")
parser$add_argument("--dump-dir", type="character", help="directory containing nodes.dmp and names.dmp to be used in database construction")
parser$add_argument("--norm-seq-table", type="character", help="normalized seq table tsv file path")
parser$add_argument("--bacterial-names", type="character", help="names of bacterial ASVs")
parser$add_argument("--trialID", type="character", help="ID number for the trial")
parser$add_argument("--add-unclassified-prefix", help="Whether to add prefix to unclassified taxa based on lowest assigned taxonomic level",
                    action="store_true",
                    default=FALSE)
parser$add_argument("--out", type="character", help="directory to store output abundance plots")


args <- parser$parse_args()



#-----------------------------
# Sequence Table Construction
#-----------------------------
#read in all ASV sequence table
seq_table <- read.delim(args$norm_seq_table, header = TRUE, row.names = 1)
'Sequence Table'
str(seq_table)

# read in positive bacterial ASV IDs
bacterial_IDs <- readLines(args$bacterial_names) # split IDs by line
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines
'Bacterial ASV IDs'
str(bacterial_IDs)

# select columns from sequence table to positive bacterial ASVs
seq_table <- seq_table[, bacterial_IDs]


#--------------------------------
# Sample Data Table Construction
#--------------------------------
sample_names <- rownames(seq_table)
sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)

# Sample Type
sample_type <- sub("\\d+[A-Za-z]*$", "", sample_info)

# Technical Rep
tech_rep <- sub(".*?(\\d+[A-Za-z]*)$", "\\1", sample_info)

sample_meta_data_df <- data.frame(
  SampleType = factor(sample_type),
  Replicate = factor(tech_rep),
  row.names = rownames(seq_table),
  stringsAsFactors = FALSE
)

#------------------------------------
# Mothur Taxonomy Table Construction
#------------------------------------
#read in taxonomy file
mothur_file <- args$mothur_file
mothur_tax_lines <- readLines(mothur_file) # split records by line
mothur_tax_lines <- mothur_tax_lines[nzchar(mothur_tax_lines)] # removes empty lines


# convert to simple data frame
mothur_tax_df <- do.call(rbind, strsplit(unlist(mothur_tax_lines), "\t")) # splits ASV ID and taxonomy record
colnames(mothur_tax_df) <- c("ASVid", "taxonomy_record")
mothur_tax_df <- as.data.frame(mothur_tax_df, stringsAsFactors = FALSE)
str(mothur_tax_df)

# split taxonomical levels
mothur_split_tax_df <- strsplit(sub(";+$", "", mothur_tax_df$taxonomy_record), ";") # removes trailing and ";" and splits record field per level at ";"

maxranks <- 9 #hardset maximum number of rank assignment (down to sub-strain)

# function for removing confidence score from level entry
clean_confidence_score <- function(x) sub("\\(.*\\)$", "", x) 


# construct taxonomy matrix via matrix transposition of split taxonomy records
mothur_tax_matrix <- t(vapply(mothur_split_tax_df, function(tax_entry) {
  tax_entry <- sapply(tax_entry, clean_confidence_score)
  # pad unassigned lower levels with NA
  if (length(tax_entry) < maxranks) {
    tax_entry <- c(tax_entry, rep(NA, maxranks - length(tax_entry)))
  }
  tax_entry
}, FUN.VALUE = character(maxranks)))

rownames(mothur_tax_matrix) <- mothur_tax_df$ASVid
colnames(mothur_tax_matrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species","Strain","Substrain")
mothur_tax_matrix <- as.matrix(mothur_tax_matrix) #convert to matrix for ease of use with phyloseq


str(mothur_tax_matrix)

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
      }
    }
  }

  return(tax_matrix)
}
kraken_tax_matrix <- taxa_level_prefix_addition(kraken_tax_matrix)
str(kraken_tax_matrix)

#------------------------------------
# Unclassified Taxa Handling
#------------------------------------
if (args$add_unclassified_prefix) {

  fill_unassigned_by_lowest <- function(tax_matrix, prefix = "Unclassified_") {
    tax_matrix <- as.matrix(tax_matrix)
    is_unassigned <- function(x) is.na(x) || x == "Unclassified"

    for (i in seq_len(nrow(tax_matrix))) {
      row <- tax_matrix[i, ]

      assigned_index <- which(!vapply(row, is_unassigned, logical(1)))
      if (length(assigned_index) == 0) next

      lowest_assigned_index <- max(assigned_index)
      lowest_assigned_rank <- colnames(tax_matrix)[lowest_assigned_index]
      lowest_assigned_value <- row[lowest_assigned_index]

      fill_value <- paste0(prefix, lowest_assigned_value)

      if (lowest_assigned_index < length(row)) {
        for (j in (lowest_assigned_index + 1):length(row)) {
          if (is_unassigned(row[j])) {
            tax_matrix[i, j] <- fill_value
          }
        }
      }
    }

    return(tax_matrix)
  }

  mothur_tax_matrix <- fill_unassigned_by_lowest(mothur_tax_matrix)
  kraken_tax_matrix <- fill_unassigned_by_lowest(kraken_tax_matrix)
}
#--------------------------------------
# Mothur Phyloseq Objects Construction
#--------------------------------------
mothur_phyloseq <- phyloseq(otu_table(seq_table, taxa_are_rows = FALSE),
                            sample_data(sample_meta_data_df),
                            tax_table(mothur_tax_matrix))
str(mothur_phyloseq)

# glom all samples based on common genus (sets all low assignments to NA)
mothur_genusGlom_phyloseq <- tax_glom(mothur_phyloseq,
                                      taxrank = "Genus",
                                      NArm = TRUE)  
                                      # generates NA column containing count of all ASVs unclassified at genus level

str(mothur_genusGlom_phyloseq)

# subset to only top 10 most prevelant genus
mothur_top10genus_phyloseq <- prune_taxa(
  names(sort(taxa_sums(mothur_genusGlom_phyloseq), decreasing = TRUE))[1:10],
  mothur_genusGlom_phyloseq
)


#--------------------------------------
# Kraken Phyloseq Objects Construction
#--------------------------------------
kraken_phyloseq <- phyloseq(otu_table(seq_table, taxa_are_rows = FALSE),
                                       sample_data(sample_meta_data_df),
                                       tax_table(kraken_tax_matrix))
str(kraken_phyloseq)

kraken_genusGlom_phyloseq <- tax_glom(kraken_phyloseq, 
                                       taxrank = "Genus",
                                       NArm = TRUE
                                       )

kraken_top10genus_phyloseq <- prune_taxa(
  names(sort(taxa_sums(kraken_genusGlom_phyloseq), decreasing = TRUE))[1:10],
  kraken_genusGlom_phyloseq
)




#--------------------------------------
# Mothur Genus Table Construction
#--------------------------------------
# constructing genus name matrix from tax table of genus glommed phyloseq object
mothur_genus_name_matrix <- as(tax_table(mothur_genusGlom_phyloseq), "matrix")[, "Genus"]

# setting taxa names in phyloseq object to genus level name matrix
taxa_names(mothur_genusGlom_phyloseq) <- mothur_genus_name_matrix

# Constructing genus count table from glommed OTU table
mothur_genus_matrix <- as.matrix(otu_table(mothur_genusGlom_phyloseq))
str(mothur_genus_matrix)

#ensure proper orientation of genus table
if (taxa_are_rows(mothur_genusGlom_phyloseq)) {
  mothur_genus_matrix <- t(mothur_genus_matrix)
}

# Add total row
mothur_genus_matrix <- rbind(mothur_genus_matrix, Total = colSums(mothur_genus_matrix))

write.table(
  mothur_genus_matrix,
  file = paste0(args$out, '/', args$trialID, "_mothur_genus_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA 
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
# Count Plots
#------------------------------
theme_set(theme_bw())

abunXtypeXsample_plot <- plot_bar(mothur_top10genus_phyloseq, "Replicate", fill = "Genus")
abunXtypeXsample_plot +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  facet_wrap(~SampleType, scales = "free_x") +
  labs(title = "Read Counts per Sample by Sample Type",
       x = "Sample",
       y = "Read Count")
ggsave(paste0(args$out, '/', args$trialID, "_mothur_countXtypeXsample.png"), width = 14, height = 12, units = "in")

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
  