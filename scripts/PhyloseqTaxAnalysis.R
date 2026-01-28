
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(readr)
library(patchwork)

#----------------------------
# Log Construction
#----------------------------  
logfile <- snakemake@log[[1]]
if (!is.null(logfile) && nzchar(logfile)) {
  logcon <- file(logfile, open = "wt")
  sink(logcon, type = "output")
  sink(logcon, type = "message")
  cat("----- R script started -----\n")
  cat("Working dir:", getwd(), "\n")
  cat("Sys.getenv PATH:", Sys.getenv("PATH"), "\n")
  flush.console()
}


#-----------------------------
# Taxonomy Table Construction
#-----------------------------

#read in taxonomy file
tax_file <- snakemake@input$taxfile
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

# function for removing level demarcation (redundant with column names)
clean_level_marker <- function(x) sub("^[^_]+__", "", x)


# construct taxonomy matrix via matrix transposition of split taxonomy records
tax_matrix <- t(vapply(split_tax_df, function(tax_entry) {
  tax_entry <- sapply(tax_entry, clean_confidence_score)
  tax_entry <- sapply(tax_entry, clean_level_marker)
  
  # pad un assigned lower levels with NA
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
seq_table <- read.delim(snakemake@input$seq_table, header = TRUE, row.names = 1)
'Sequence Table'
str(seq_table)

# read in positive bacterial ASV IDs
bacterial_IDs <- readLines(snakemake@input$bacterial_names) # split IDs by line
bacterial_IDs <- bacterial_IDs[nzchar(bacterial_IDs)] # removes empty lines
'Bacterial ASV IDs'
str(bacterial_IDs)

# select columns from sequence table to positive bacterial ASVs
seq_table <- seq_table[, bacterial_IDs]


#--------------------------------
# Sample Data Table Construction
#--------------------------------
sample_names <- rownames(seq_table)
cleaned_sample_names <- sapply(strsplit(sample_names, "\\."), `[`, 2)
rownames(seq_table) <- cleaned_sample_names

# Bacterial Spike Concentration
sample_info <- sapply(strsplit(sample_names, "_"), `[`, 3)
bact_conc <- sub("b.*$", "", sample_info)

# Number of technical replicate
tech_rep <- sub(".*bac([0-9]+)$*", "\\1", sample_info)

sample_meta_data_df <- data.frame(
  Concentration = factor(bact_conc),
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
all_sample_phyloseq_genus_glom <- tax_glom(all_sample_phyloseq, taxrank = "Genus")

# subset to only top 10 most prevelant genus
all_sample_phyloseq_top10g <- prune_taxa(
  names(sort(taxa_sums(all_sample_phyloseq_genus_glom), decreasing = TRUE))[1:10],
  all_sample_phyloseq_genus_glom
)

#------------------------------
# General Abundance Plots
#------------------------------
theme_set(theme_bw())

abunXconc_plot <- plot_bar(all_sample_phyloseq_top10g, x="Concentration", fill = "Genus")
abunXconc_plot + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  )
ggsave(snakemake@output$abunXconc_plot)

abunXsample_plot <- plot_bar(all_sample_phyloseq_top10g, fill = "Genus")
abunXsample_plot + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  )
ggsave(snakemake@output$abunXsample_plot)


abunXconcXsample_plot <- plot_bar(all_sample_phyloseq_top10g, "Replicate", fill = "Genus", facet_grid = ~Concentration)
abunXconcXsample_plot +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  )

ggsave(snakemake@output$abunXconcXsample_plot)

#------------------------------
# Top 20 Analysis
#------------------------------
# top20_names <- names(sort(taxa_sums(all_sample_phyloseq), decreasing=TRUE))[1:20]
# all_sample_phyloseq.top20 <- transform_sample_counts(all_sample_phyloseq, function(OTU) OTU/sum(OTU))
# all_sample_phyloseq.top20 <- prune_taxa(top20_names, all_sample_phyloseq.top20)
# plot_bar(all_sample_phyloseq.top20, x="Concentration", fill="Family")
# ggsave(snakemake@output$abundance_plot)


#------------------------------
# Bray-Curtis Distances Analysis
#------------------------------
# sample.phyloseq.prop <- transform_sample_counts(sample.phyloseq, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(sample.phyloseq.prop, method="NMDS", distance="bray")


#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}