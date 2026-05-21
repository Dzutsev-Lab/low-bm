library(dada2)
library(Biostrings)
library(digest)
library(ggplot2)
library(argparse)

#---------------------------------
# Execution Arguments
#---------------------------------
parser <- ArgumentParser()

#input files
parser$add_argument("--fqs", type="character", nargs='+', help="List of fastq file pathes to be used in denoising") # nargs='+' tells it to look for one or more arguments
parser$add_argument("--sample-names", type="character", help="File path to .names file with sample names, one per line")

#output files
parser$add_argument("--filtered-fqs", type="character", nargs='+', help="List of fastq file pathes to store filtered reads.")
parser$add_argument("--err-plt", type="character", help="File path to store error modeling plot.")
parser$add_argument("--filt-counts", type="character", help="File path to store filter stage read counts tsv.")
parser$add_argument("--asv-fa", type="character", help="File path to store stable ASV IDs")
parser$add_argument("--seq-table", type="character", help="File path to store sequence table")
parser$add_argument("--asv-map", type="character", help="File path to store asv ID to sequence map")

#trimAndFilerter Paramters
parser$add_argument("--chunk-size", type="integer")
parser$add_argument("--truncLen", type="integer")
parser$add_argument("--primerLen", type="integer")
parser$add_argument("--maxN", type="integer")
parser$add_argument("--maxEE", type="integer")
parser$add_argument("--truncQ", type="integer")

#computational resource parameters
parser$add_argument("--threads", type="integer")

args <- parser$parse_args()

# ------ params from snakemake ----

# filterAndTrim PARAMETERS
chunk_size <- args$chunk_size
truncLen <-  args$truncLen
primerLen <- args$primerLen
maxN <- args$maxN
maxEE <- args$maxEE
truncQ <- args$truncQ



#----------------------------
# File Accounting
#----------------------------  
#internal to R script, snakemake should be standing sentinal
# ---- discover fastqs ----
fq.files <- args$fqs

# ---- import sample names ----
sample.names <- readLines(args$sample_names)
sample.names <- sample.names[nzchar(sample.names)] # removes empty lines


# --- import outpaths for filtered fastqs ----
filtered.fq.files <- args$filtered_fqs

#check correct number of filtered fastq compared to input fastqs
if (length(fq.files) != length(filtered.fq.files)) {
  stop("Mismatch: fqs=", length(fq.files), " filtered_fqs=", length(filtered.fq.files))
}


# --- set names for fastq vectors to sample names ---
names(fq.files) <- sample.names
names(filtered.fq.files) <- sample.names

#----------------------------
# Filtering and Trimming
#----------------------------  
out_mat <- filterAndTrim( 
  fwd = fq.files, 
  filt = filtered.fq.files,
  trimLeft = primerLen,
  truncLen = truncLen,
  maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
  rm.phix = TRUE,
  compress = FALSE, 
  multithread = args$threads,
  n = chunk_size
)
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

error.model <- learnErrors(nonzero.fq.files, multithread = args$threads)

# ---- plot error models ----
png(filename = args$err_plt, width=800, height=600)
plotErrors(error.model, nominalQ = TRUE)
dev.off()

#----------------------------
# Unique Sequence Denoising
#----------------------------
denoised.reads <- dada(nonzero.fq.files, err=error.model, multithread = args$threads)

# ---- Construct Sequence Table ----
seqtab.denoise <- makeSequenceTable(denoised.reads)

# ---- Remove Chimeras ----
seqtab.nochim <- removeBimeraDenovo(seqtab.denoise, method = "consensus", multithread = args$threads)

# ---- Add Empty Entries for Completely Filetered Samples ----
missing.samples <- setdiff(sample.names, rownames(seqtab.nochim))

if (length(missing.samples) > 0) {

  #adding entries to both no-chimera seq table and just denoised seq table to ensure same size for read count tsv construction
  denoised.zero.matrix <- matrix(
    0,
    nrow = length(missing.samples),
    ncol = ncol(seqtab.denoise),
    dimnames = list(missing.samples, colnames(seqtab.denoise))
  )
  seqtab.denoise <- rbind(seqtab.denoise, denoised.zero.matrix)

  nochim.zero.matrix <- matrix(
    0,
    nrow = length(missing.samples),
    ncol = ncol(seqtab.nochim),
    dimnames = list(missing.samples, colnames(seqtab.nochim))
  )
  seqtab.nochim <- rbind(seqtab.nochim, nochim.zero.matrix)

  # Reordering sequence tables to ensure they match row ordering with other sample matrices
  seqtab.denoise <- seqtab.denoise[sample.names, , drop = FALSE]
  seqtab.nochim  <- seqtab.nochim[sample.names, , drop = FALSE]
}



#----------------------------
# Sequence Tracking
#----------------------------

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
            file = args$filt_counts, 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE)

#----------------------------
# Classification FASTA Prep
#----------------------------
#stable ASV IDs (not just DNA sequence)
seqs <- colnames(seqtab.nochim)
ASV_IDs <- paste0(
  "ASV_", 
  vapply(
    seqs,
    digest,
    character(1),
    algo = "md5",
    serialize = FALSE
  )
)
dna_strings <- DNAStringSet(seqs)
names(dna_strings) <- ASV_IDs
writeXStringSet(dna_strings, 
                filepath = args$asv_fa,
                format = "fasta")

#----------------------------------------------
# Export Sequence Table Relabeled with ASV IDs
#----------------------------------------------
colnames(seqtab.nochim) <- ASV_IDs
write.table(seqtab.nochim, 
            file = args$seq_table, 
            sep = '\t', 
            quote = FALSE, 
            col.names = NA) #check to see if this is not mutating the matrix
            
#----------------------------------------------
# Export ASV Sequence Lookup Table
#----------------------------------------------
asv_map <- data.frame(
  ASV_ID = ASV_IDs,
  Sequence = seqs,
  stringsAsFactors = FALSE
)

write.table(
  asv_map,
  file = args$asv_map,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)