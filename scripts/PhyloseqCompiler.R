library(argparse)
library(phyloseq)

parser <- ArgumentParser()

parser$add_argument("--trial-list",
                    type = "character",
                    help = "list of trials to compile into single phyloseq object file path")
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

# Save output
if (!dir.exists(args$out)) {
    dir.create(args$out, recursive = TRUE)
}
out_file <- file.path(args$out, "CompPhyseq.RData")
save(physeq, file = out_file)

message("Merged phyloseq object saved to: ", out_file)