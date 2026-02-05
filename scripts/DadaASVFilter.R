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
library(Biostrings)
library(ggplot2)
library(DECIPHER)

# ------ params from snakemake ----

# filterAndTrim PARAMETERS
chunk_size <- snakemake@params$chunk_size
truncLen <-  snakemake@params$truncLen
primerLen <- snakemake@params$primerLen
maxN <- snakemake@params$maxN
maxEE <- snakemake@params$maxEE
truncQ <- snakemake@params$truncQ



#----------------------------
# File Accounting
#----------------------------  
#internal to R script, snakemake should be standing sentinal
# ---- discover fastqs ----
fq.files <- snakemake@input$umi_dedup_reads

# ---- import sample names ----
sample.names <- readLines(snakemake@input$sample_names)
sample.names <- sample.names[nzchar(sample.names)] # removes empty lines


# --- import outpaths for filtered fastqs ----
filtered.fq.files <- snakemake@output$filtered_reads

# --- set names for fastq vectors to sample names ---
names(fq.files) <- sample.names
names(filtered.fq.files) <- sample.names

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
                                  trimLeft = primerLen,
                                  truncLen = truncLen,
                                  maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                                  rm.phix = TRUE,
                                  compress = FALSE, 
                                  multithread = TRUE)
  outs_list[[i]] <- filter.result

}
out_mat <- do.call(rbind, outs_list)
rownames(out_mat) <- names(fq.files)

# ---- find filtered, non-empty fastqs ----
nonzero.filt <- out_mat[, "reads.out"] > 0
nonzero.fq.files <- filtered.fq.files[nonzero.filt]
message(sprintf("Filtered+Trimmed, non-empty files: %d / %d", length(nonzero.fq.files), length(filtered.fq.files)))

zero.fq.idx <- which(!nonzero.filt)
if (length(zero.fq.idx) > 0) {
  message(sprintf("Creating %d empty placeholder filtered fastq files for samples with 0 reads after filtering", length(zero.fq.idx)))
  for (k in zero.fq.idx) {
    empty.fq.path <- filtered.fq.files[k]
    if (file.exists(empty.fq.path)) {
      # truncate to zero bytes, how does this work???
      con <-file(empty.fq.path, open = "wb")
      close(con)
    } else {
      file.create(empty.fq.path)
    }
  }
}

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
seqtab.nochim <- removeBimeraDenovo(seqtab.denoise, method = "consensus", multithread = TRUE)

# ---- Add Empty Entries for Completely Filetered Samples ----
missing.samples <- setdiff(sample.names, rownames(seqtab.nochim))

if (length(missing.samples) > 0) {
  zero.matrix <- matrix(
    0,
    nrow = length(missing.samples),
    ncol = ncol(seqtab.nochim),
    dimnames = list(missing.samples, colnames(seqtab.nochim))
  )
  
  seqtab.nochim <- rbind(seqtab.nochim, zero.matrix)
}

#reordering to match original sample ordering
seqtab.nochim <- seqtab.nochim[sample.names, , drop = FALSE]

#----------------------------
# Sequence Tracking
#----------------------------
getN <- function(x) sum(getUniques(x))
tracker <- cbind(out_mat, 
                 rowSums(seqtab.denoise), 
                 rowSums(seqtab.nochim))
colnames(tracker) <- c("deduped", "filteredAndTrimmed", "denoised", "chimera-filtered")
rownames(tracker) <- sample.names

tracker_df <- data.frame(
  SampleID = rownames(tracker),
  tracker,
  check.names = FALSE
)

write.table(tracker_df, 
            file = snakemake@output$filter_stage_counts, 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE,
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

#----------------------------------------------
# Export Sequence Table Relabeled with ASV IDs
#----------------------------------------------
colnames(seqtab.nochim) <- ASV_IDs
write.table(seqtab.nochim, 
            file = snakemake@output$seq_table, 
            sep = '\t', 
            quote = FALSE, 
            col.names = NA) #check to see if this is not mutating the matrix

#----------------------------
# Log Close
#----------------------------
if (!is.null(logfile) && nzchar(logfile)) {
  sink(type = "message"); sink(type = "output")
  close(logcon)
}