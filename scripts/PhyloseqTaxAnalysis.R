# #----------------------------
# # Log Construction
# #----------------------------  
# logfile <- snakemake@log[[1]]
# if (!is.null(logfile) && nzchar(logfile)) {
#   logcon <- file(logfile, open = "wt")
#   sink(logcon, type = "output")
#   sink(logcon, type = "message")
#   cat("----- R script started -----\n")
#   cat("Working dir:", getwd(), "\n")
#   cat("Sys.getenv PATH:", Sys.getenv("PATH"), "\n")
#   flush.console()
# }

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(argparse)
library(dplyr)

# Keep ggplot from producing Rplots.pdf
if(!interactive()) pdf(NULL)

parser <- ArgumentParser()

parser$add_argument("--tax-file", type="character", help="taxonomy classification file path")
parser$add_argument("--kraken-file", type="character", help="kraken taxonomy classification file path")
parser$add_argument("--dump-dir", type="character", help="directory containing nodes.dmp and names.dmp to be used in database construction")
parser$add_argument("--norm-seq-table", type="character", help="normalized seq table tsv file path")
parser$add_argument("--bacterial-names", type="character", help="names of bacterial ASVs")
parser$add_argument("--abund-plot-dir", type="character", help="directory to store output abundance plots")

args <- parser$parse_args()

#-----------------------------
# Taxonomy Table Construction
#-----------------------------

#read in taxonomy file
tax_file <- args$tax_file
tax_lines <- readLines(tax_file) # split records by line
tax_lines <- tax_lines[nzchar(tax_lines)] # removes empty lines


# convert to simple data frame
tax_df <- do.call(rbind, strsplit(unlist(tax_lines), "\t")) # splits ASV ID and taxonomy record
colnames(tax_df) <- c("ASVid", "taxonomy_record")
tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
str(tax_df)

# split taxonomical levels
split_tax_df <- strsplit(sub(";+$", "", tax_df$taxonomy_record), ";") # removes trailing and ";" and splits record field per level at ";"

maxranks <- 9 #hardset maximum number of rank assignment (down to sub-strain)

# function for removing confidence score from level entry
clean_confidence_score <- function(x) sub("\\(.*\\)$", "", x) 


# construct taxonomy matrix via matrix transposition of split taxonomy records
tax_matrix <- t(vapply(split_tax_df, function(tax_entry) {
  tax_entry <- sapply(tax_entry, clean_confidence_score)

  
  # pad unassigned lower levels with NA
  if (length(tax_entry) < maxranks) {
    tax_entry <- c(tax_entry, rep(NA, maxranks - length(tax_entry)))
  }
  tax_entry
}, FUN.VALUE = character(maxranks)))

rownames(tax_matrix) <- tax_df$ASVid
colnames(tax_matrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species","Strain","Substrain")
tax_matrix <- as.matrix(tax_matrix) #convert to matrix for ease of use with phyloseq
str(tax_matrix)

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


#--------------------------------
# Phyloseq Objects Construction
#--------------------------------
all_sample_phyloseq <- phyloseq(otu_table(seq_table, taxa_are_rows = FALSE),
                                sample_data(sample_meta_data_df),
                                tax_table(tax_matrix))
str(all_sample_phyloseq)



# glom all samples based on common genus (sets all low assignments to NA)
all_sample_phyloseq_genus_glom <- tax_glom(all_sample_phyloseq,
                                           taxrank = "Genus",
                                           NArm = TRUE) # generates NA column containing count of all ASVs unclassified at genus level
                                                        # shouldn't do anything with mothur confidence threshold at 0 (all taxa level given assignment)
str(all_sample_phyloseq_genus_glom)

# constructing genus name matrix from tax table of genus glommed phyloseq object
genus_name_matrix <- as(tax_table(all_sample_phyloseq_genus_glom), "matrix")[, "Genus"]
genus_name_matrix[is.na(genus_name_matrix) | genus_name_matrix == ""] <- "Unassigned" # dealing with Unassigned column label

# setting taxa names in phyloseq object to genus level name matrix
taxa_names(all_sample_phyloseq_genus_glom) <- genus_name_matrix

# Constructing genus count table from glommed OTU table
genus_matrix <- as.matrix(otu_table(all_sample_phyloseq_genus_glom))
str(genus_matrix)

#ensure proper orientation of genus table
if (taxa_are_rows(all_sample_phyloseq_genus_glom)) {
  genus_matrix <- t(genus_matrix)
}

write.table(
  genus_matrix,
  file = file.path(args$abund_plot_dir, "genus_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA 
)



# subset to only top 10 most prevelant genus
all_sample_phyloseq_top10g <- prune_taxa(
  names(sort(taxa_sums(all_sample_phyloseq_genus_glom), decreasing = TRUE))[1:10],
  all_sample_phyloseq_genus_glom
)

#------------------------------
# General Abundance Plots
#------------------------------
theme_set(theme_bw())

abunXtype_plot <- plot_bar(all_sample_phyloseq_top10g, x="SampleType", fill = "Genus")
abunXtype_plot + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  )
ggsave(file.path(args$abund_plot_dir, "abunXtype.png"))

abunXsample_plot <- plot_bar(all_sample_phyloseq_top10g, fill = "Genus")
abunXsample_plot + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  )
ggsave(file.path(args$abund_plot_dir, "abunXsample.png"))


abunXtypeXsample_plot <- plot_bar(all_sample_phyloseq_top10g, "Replicate", fill = "Genus")
abunXtypeXsample_plot +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  facet_wrap(~SampleType, scales = "free_x") +
  labs(title = "Read Counts per Sample by Sample Type",
       x = "Sample",
       y = "Read Count")
ggsave(file.path(args$abund_plot_dir, "abunXtypeXsample.png"), width = 10, height = 12, units = "in")




#-----------------------------
# Kraken Phyloseq Analysis
#-----------------------------
library(taxonomizr)
# Taxonomy Database File Paths
names_dmp <- file.path(args$dump_dir, "names.dmp")
nodes_dmp <- file.path(args$dump_dir, "nodes.dmp")
sql_db <- file.path(args$abund_plot_dir, paste0("tax_db_", basename(args$dump_dir), ".sqlite"))

# Taxonomy Database Construction
read.names.sql(names_dmp, sql_db, overwrite = TRUE)
read.nodes.sql(nodes_dmp, sql_db, overwrite = TRUE)

# Import ASVid -> TaxID info from .kraken2 file
kraken_info <- read.delim(args$kraken_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(kraken_info)[1:3] <- c("status", "ASVid", "taxid")

# trimming kraken information to ASVid and TaxID for all with positive classifications
kraken_info <- kraken_info[kraken_info$status == "C" & kraken_info$taxid != 0, c("ASVid", "taxid")]
kraken_info$taxid <- as.integer(kraken_info$taxid)

desiredTaxa <- c("superkingdom","phylum","class","order","family","genus","species")

tax_df <- getTaxonomy(ids = kraken_info$taxid, 
                      sqlFile = sql_db, 
                      desiredTaxa = desiredTaxa)
tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
tax_df$ASVid <- kraken_info$ASVid

tax_df <- rename( tax_df, 
                  Kingdom = superkingdom,
                  Phylum  = phylum,
                  Class   = class,
                  Order   = order,
                  Family  = family,
                  Genus   = genus,
                  Species = species)


tax_matrix <- as.matrix(tax_df[, c("Kingdom","Phylum","Class","Order","Family","Genus","Species")])
rownames(tax_matrix) <- tax_df$ASVid
str(tax_matrix)


kraken_all_sample_phyloseq <- phyloseq(otu_table(seq_table, taxa_are_rows = FALSE),
                                       sample_data(sample_meta_data_df),
                                       tax_table(tax_matrix))
str(kraken_all_sample_phyloseq)

kraken_all_sample_phyloseq_genus_glom <- tax_glom(kraken_all_sample_phyloseq, 
                                                  taxrank = "Genus",
                                                  NArm = TRUE)

kraken_all_sample_phyloseq_top10g <- prune_taxa(
  names(sort(taxa_sums(kraken_all_sample_phyloseq_genus_glom), decreasing = TRUE))[1:10],
  kraken_all_sample_phyloseq_genus_glom
)

kraken_abunXtypeXsample_plot <- plot_bar(kraken_all_sample_phyloseq_top10g, "Replicate", fill = "Genus")
kraken_abunXtypeXsample_plot +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  facet_wrap(~SampleType, scales = "free_x") +
  labs(title = "Read Counts per Sample by Sample Type",
       x = "Sample",
       y = "Read Count")
ggsave(file.path(args$abund_plot_dir, "kraken_abunXtypeXsample.png"), width = 10, height = 12, units = "in")

# #----------------------------
# # Log Close
# #----------------------------
# if (!is.null(logfile) && nzchar(logfile)) {
#   sink(type = "message"); sink(type = "output")
#   close(logcon)
# }