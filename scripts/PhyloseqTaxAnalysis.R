
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(readr)

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


#--------------------------
# Taxonomy File Parsing
#--------------------------

#read in taxonomy file
tax_file <- snakemake@input$taxfile

# split tax file by line
tax_lines <- readLines(tax_file)
tax_lines <- tax_lines[nzchar(tax_lines)] # removes empty lines

# convert to simple data frame
tax_df <- do.call(rbind, strsplit(unlist(tax_lines), "\t")) # splits ASV ID and taxonomy record
colnames(tax_df) <- c("ASVid", "taxonomy_record")
tax_df <- as.data.frame(tax_df, stringAsFactors = FALSE)

split_tax <- strsplit(tax_df$taxonomy_record, ";") # splits record filed per level at ";"

maxranks <- max(sapply(split_tax, length))

# function for removing confidence score from level entry
clean_confidence_score <- function(x) sub("\\(.*\\)$", "", x) 

tax_matrix <- t(sapply(split_tax, function(tax_entry) {
  clean_tax_levels <- sapply(tax_entry, clean_confidence_score)
  
  # pad un assigned lower levels with NA
  if (length(clean_tax_levels) < maxranks) {
    clean_tax_levels <- c(clean_tax_levels, rep(NA, maxranks - length(clean_tax_levels)))
  }
}))

rownames(tax_matrix) <- tax_df$ASVid
colnames(tax_matrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species","Strain","Substrain")
tax_matrix <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
head(tax_matrix)




# seqtab.nochim <- read.delim("../Exp_Output/Subset/clean.fastq/Filtered/SeqTable.tsv", header = TRUE, row.names = 1)
# taxadf <- read.delim("../Exp_Output/Subset/clean.fastq/Filtered/ASVTaxonomy.tsv", header = TRUE, row.names = 1)
# taxa <- as.matrix(taxadf)
# 
# theme_set(theme_bw())
# 
# samples.out <- rownames(seqtab.nochim)
# bact_conc <- sapply(strsplit(samples.out, "_"), `[`, 3)
# bact_conc <- sub("b.*$", "", bact_conc)
# 
# sampledf <- data.frame(Concentration=bact_conc)
# rownames(sampledf) <- samples.out
# 
# 
# #construct phyloseq object from dada2 output
# sample.phyloseq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
#                              sample_data(sampledf),
#                              tax_table(taxa))
# 
# # Renaming sequences as numbered ASV ID (instead of whole sequence)
# dna <- Biostrings::DNAStringSet(taxa_names(sample.phyloseq))
# names(dna) <- taxa_names(sample.phyloseq)
# sample.phyloseq <- merge_phyloseq(sample.phyloseq, dna)
# taxa_names(sample.phyloseq) <- paste0("ASV", seq(ntaxa(sample.phyloseq)))
# sample.phyloseq
# 
# # Transform data to proportions as appropriate for Bray-Curtis distances
# sample.phyloseq.prop <- transform_sample_counts(sample.phyloseq, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(sample.phyloseq.prop, method="NMDS", distance="bray")
# 
# top20 <- names(sort(taxa_sums(sample.phyloseq), decreasing=TRUE))[1:175]
# sample.phyloseq.top20 <- transform_sample_counts(sample.phyloseq, function(OTU) OTU/sum(OTU))
# sample.phyloseq.top20 <- prune_taxa(top20, sample.phyloseq.top20)
# plot_bar(sample.phyloseq.top20, x="Concentration", fill="Family")
# 
# # #----------------------------
# # # Phyloseq Abundance Analysis
# # #----------------------------
# # theme_set(theme_bw())
# # 
# # samples.out <- rownames(seqtab.nochim)
# # bact_conc <- sapply(strsplit(samples.out, "_"), `[`, 3)
# # bact_conc <- sub("b.*$", "", bact_conc)
# # 
# # sampledf <- data.frame(Concentration=bact_conc)
# # rownames(sampledf) <- samples.out
# # 
# # 
# # #construct phyloseq object from dada2 output
# # sample.phyloseq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
# #                             sample_data(sampledf),
# #                             tax_table(taxa))
# # 
# # # Renaming sequences as numbered ASV ID (instead of whole sequence)
# # dna <- Biostrings::DNAStringSet(taxa_names(sample.phyloseq))
# # names(dna) <- taxa_names(sample.phyloseq)
# # sample.phyloseq <- merge_phyloseq(sample.phyloseq, dna)
# # taxa_names(sample.phyloseq) <- paste0("ASV", seq(ntaxa(sample.phyloseq)))
# # sample.phyloseq
# # 
# # # Top 20 most abundant 
# # top20 <- names(sort(taxa_sums(sample.phyloseq), decreasing=TRUE))[1:20]
# # sample.phyloseq.top20 <- transform_sample_counts(sample.phyloseq, function(OTU) OTU/sum(OTU))
# # sample.phyloseq.top20 <- prune_taxa(top20, sample.phyloseq.top20)
# # plot_bar(sample.phyloseq.top20, x="Concentration", fill="Family")
# # ggsave(file.path(out_dir,"abundanceXconc.png"))


#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}