#!/usr/bin/env Rscript

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

library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(DECIPHER)

# ------ params from snakemake ----
# DIRECTORIES
in_dir <- snakemake@params$bacterial_dir
denoise_out_dir <- snakemake@params$denoised_dir
filterAndTrim_out_dir <- file.path(denoise_out_dir, "filteredAndTrimmed")
count_out_dir <- snakemake@params$count_dir

# filterAndTrim PARAMETERS
chunk_size <- snakemake@params$chunk_size
truncLen <-  snakemake@params$truncLen
primerLen <- snakemake@params$primerLen
maxN <- snakemake@params$maxN
maxEE <- snakemake@params$maxEE
truncQ <- snakemake@params$truncQ

# TAXONOMY REFERENCE
taxRef <- snakemake@params$taxRef



#----------------------------
# File Accounting
#----------------------------  
#internal to R script, snakemake should be standing sentinal
# ---- discover fastqs ----
fq.files <- sort(snakemake@input$bacterial_reads)

# ---- extract sample names ----
sample.names <- sapply(strsplit(basename(fq.files), "\\."), `[`, 2)
names(fq.files) <- sample.names

# --- create filtered fastqs ----
filtered.fq.files <- snakemake@output$filtered_reads

#----------------------------
# Filtering and Trimming
#----------------------------  
chunk_idx <- split(seq_along(fq.files), ceiling(seq_along(fq.files) / chunk_size))
outs_list <- vector("list", length(chunk_idx))

for (i in seq_along(chunk_idx)) {
  j <- chunk_idx[[i]]
  message(sprintf('[chunk %d/%d] processing %d samples', i, length(chunk_idx), length(j)))
  filter.result <- filterAndTrim( fwd = fq.files[j], 
                                  filt = filtered.fq.files[j],
                                  truncLen = truncLen,
                                  trimLeft = primerLen,
                                  maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                                  rm.phix = TRUE,
                                  compress = FALSE, 
                                  multithread = FALSE)
  outs_list[[i]] <- filter.result

}
out_mat <- do.call(rbind, outs_list)
rownames(out_mat) <- names(fq.files)

# ---- find filtered, non-empty fastqs ----
nonzero.filt <- out_mat[, "reads.out"] > 0
nonzero.fq.files <- filtered.fq.files[nonzero.filt]
message(sprintf("Filtered+Trimmed, non-empty files: %d / %d", length(nonzero.fq.files), length(filtered.fq.files)))


#----------------------------
# Error Modeling
#----------------------------

# ---- Insufficient FASTQ check ----
if (length(nonzero.fq.files) < 2) stop("Not enough filtered+trimmed FASTQs to learn error model.")

error.model <- learnErrors(nonzero.fq.files, multithread = TRUE)

# ---- plot error models ----
png(filename = snakemake@output$seq_err_plot, width=800, height=600)
plotErrors(error.model, nominalQ = TRUE)
dev.off()

#----------------------------
# Unique Sequence Denoising
#----------------------------
denoised.reads <- dada(nonzero.fq.files, err=error.model, multithread=TRUE)

# ---- Construct Sequence Table ----
seqtab.denoise <- makeSequenceTable(denoised.reads)

# ---- Remove Chimeras ----
# TODO: ensure that sequence table is beig exported as tsv properly
seqtab.nochim <- removeBimeraDenovo(seqtab.denoise, method = "consensus", multithread = TRUE)
write.table(seqtab.nochim, 
            file = snakemake@output$seq_table, 
            sep = '\t', 
            quote = FALSE, 
            col.names = NA) #check to see if this is not mutating the matrix

#----------------------------
# Sequence Tracking
#----------------------------
getN <- function(x) sum(getUniques(x))
tracker <- cbind(out_mat, rowSums(seqtab.denoise), rowSums(seqtab.nochim))
colnames(tracker) <- c("ID", "filteredAndTrimmed", "denoised", "chimera-filtered")
rownames(tracker) <- sample.names
write.table(tracker, 
            file = snakemake@output$filter_stage_counts, 
            sep = '\t', 
            quote = FALSE, 
            col.names = TRUE)

#----------------------------
# Classification FASTA Prep
#----------------------------
#stable ASV IDs (not just DNA sequence)
ASV_IDs <- paste0("ASV", seq_len(ncol(seqtab.nochim)))
dna_strings <- DNAStringSet(colnames((seqtab.nochim)))
names(dna_strings) <- ASV_IDs
writeXStringSet(dna_strings, 
                filepath = snakemake@output$rep_asv_fasta,
                format = "fasta")


#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}