library(dada2)
library(readr)
library(Biostrings)
library(DECIPHER)

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


# Seq Table tsv read-in
seqtab.nochim <- as.matrix(read.delim(snakemake@input$seq_table, 
                                      header = TRUE, 
                                      sep = "\t",
                                      row.names = 1,
                                      check.names = FALSE)) # ensures that DNA sequence names are not mangled by R on import
# ensures count values in sequence table are stored as integers
storage.mode(seqtab.nochim) <- "integer"

# TAXONOMY REFERENCE
training_ref <- snakemake@input$tax_ref

training_ref_dna <- readDNAStringSet(training_ref)

# parse the headers to obtain a taxonomy for training set
s <- strsplit(names(training_ref_dna), " ")
genus <- sapply(s, `[`, 1)
species <- sapply(s, `[`, 2)
training_taxonomy <- paste("Root", genus, species, sep="; ")


#-----------------------
# Classifier Training
#-----------------------
# TODO: find pre-trained model that can be plugged in instead of retraining each time
trainingSet <- LearnTaxa(train = training_ref_dna,
                         taxonomy = training_taxonomy)
png(filename = "/Exp_Output/Subset_Samples/training_plot.png", width=800, height=600)
plot(trainingSet)
dev.off()


#----------------------------
# Taxa Assignment
#----------------------------
# dna <- DNAStringSet(getSequences(seqtab.nochim))
# # struggling to download training set from DECIPHER website
# # using training set provided for illustrative purposes only as a part of DECIPHER package
# 
# ids <- IdTaxa(dna, trainingSet, strand = "top", processors = NULL, verbose = FALSE)
# ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# taxa <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# 
# colnames(taxa) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
# 
# write.table(taxa, file = snakemake@output$taxonomy_table, sep = '\t', quote = FALSE, col.names = TRUE)
# 


#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}