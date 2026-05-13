library(argparse)
library(phyloseq)

parser <- ArgumentParser()

parser$add_argument("--trial-list",
                    type = "character",
                    help = "list of trials to compile into single phyloseq object file path")
parser$add_argument("--techrep-avg",
                    action = "store_true",
                    default = FALSE,
                    help = "Optional flag to average OTU counts between technical replicates")
parser$add_argument("--out",
                    type = "character",
                    help = "desired output directory for R session image with complied phyloseq object")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing the trial folders [default: Exp_Output]")

args <- parser$parse_args()

trials <- readLines(args$trial_list)
trials <- trimws(trials)
trials <- trials[nzchar(trials)]

if (length(trials) == 0) {
  stop("No trial names found in --trial-list.")
}

#Helper to construct the file path to physeq image based on trial naming
trial_to_rdata <- function(trial, base_dir = "Exp_Output") {
  trial_dir <- file.path(base_dir, trial)
  trialID <- sub("_.*$", "", trial)
  file.path(trial_dir, paste0(trialID, "_physeq.RData"))
}

#Help to average techincal replicate count values 
average_by_techrep <- function(physeq) {
    meta_df <- as(sample_data(physeq), "data.frame")
    otu_mat <- as(otu_table(physeq), "matrix")

    if (taxa_are_rows(physeq)) {
        otu_mat <- t(otu_mat)
    }

    # Make vector for grouping by SampleID (should be identical between technical replicates)
    meta_df$SampleID <- as.character(meta_df$SampleID)
    group <- meta_df$SampleID

    # Make explicit factor for SampleID grouping
    group_levels <- unique(group)
    group_factor <- factor(group, levels = group_levels)
    group_counts <- table(group_factor)
    print(str(group_counts))

    # avergae rows within each unique Sample ID
    avg_otu_mat <- rowsum(otu_mat, group = group_factor, reorder = FALSE)
    print(str(avg_otu_mat))
    avg_otu_mat <- sweep(
        avg_otu_mat,
        1,
        as.numeric(group_counts[rownames(avg_otu_mat)]),
        FUN = "/"
    )
    

    # filter meta data to representative sample rows (1 for each technical replicate pair)
    avg_meta_df <- meta_df[match(rownames(avg_otu_mat), meta_df$SampleID), , drop = FALSE]
    rownames(avg_meta_df) <- rownames(avg_otu_mat)
    
    phyloseq(
        otu_table(as.matrix(avg_otu_mat), taxa_are_rows = FALSE),
        sample_data(avg_meta_df),
        tax_table(physeq)
    )
}

# list to store physeq objects before merging
physeq_list <- list()

# add loaded phyloseq objects to list
for (i in seq_along(trials)) {
  trial <- trials[i]
  rdata_path <- trial_to_rdata(trial, args$base_dir)

  if (!file.exists(rdata_path)) {
    stop(sprintf("Missing RData file for trial '%s': %s", trial, rdata_path))
  }

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)

  if (!exists("physeq", envir = e, inherits = FALSE)) {
    stop(sprintf("File does not contain an object named 'physeq': %s", rdata_path))
  }

  physeq_obj <- get("physeq", envir = e)

  if (!inherits(physeq_obj, "phyloseq")) {
    stop(sprintf("Object 'physeq' in %s is not a phyloseq object.", rdata_path))
  }

  physeq_list[[i]] <- physeq_obj
}

# Merge all phyloseq objects
if (length(physeq_list) == 1) {
  CompPhyseq <- physeq_list[[1]]
} else {
  CompPhyseq <- Reduce(function(x, y) merge_phyloseq(x, y), physeq_list)
}

physeq <- CompPhyseq

if (args$techrep_avg) {
  physeq <- average_by_techrep(physeq)
}

# Save output
if (!dir.exists(args$out)) {
    dir.create(args$out, recursive = TRUE)
}
out_file <- file.path(args$out, "CompPhyseq.RData")
save(physeq, file = out_file)

message("Merged phyloseq object saved to: ", out_file)