library(phyloseq)
library(Biostrings)
library(data.table)
library(tibble)
library(stringr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--DA-comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Choice of differential abunance comprisons to pull significant taxa from [CellLineControltoTumor, CellLineControltoNontumor, NegativeControl, PatientSample, TumorType]")
parser$add_argument("--DA-method",
                    type = "character",
                    help = "Differential abundance method used to perform comparisons")
parser$add_argument("--in-trial",
                    type = "character",
                    help = "Trial with differential abundance analyses to pull significant taxa from")
parser$add_argument("--all-asv-fasta",
                    type = "character",
                    help = "FASTA file with all ASV sequences (ASV IDs as headers)")
parser$add_argument("--out-trial",
                    type = "character",
                    help = "Directory to store output manifests and select ASV fasta's within base directory")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing the trial directory")

args <- parser$parse_args()


#-----------------------------
# Helpers
#-----------------------------
get_asv_abundance <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    rowSums(otu)
  } else {
    colSums(otu)
  }
}


#-----------------------------
# Main export function
#-----------------------------
export_significant_genus_asvs <- function(ps,
                                          asv_fasta,
                                          sig_genera_by_comparison,
                                          out_dir,
                                          taxa_level_col = "Genus") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Read All ASV FASTA
    asv_fasta <- readDNAStringSet(asv_fasta)

    # Extract taxonomy (first as matrix then convert to data frame)
    tax_df <- as.data.frame(as(tax_table(ps), "matrix")) |>
        rownames_to_column("ASVid")

    if (!taxa_level_col %in% colnames(tax_df)) {
        stop("taxa_level_col = '", taxa_level_col, "' not found in tax_table(ps).")
    }

    # Extract abundance
    abund <- get_asv_abundance(ps)
    abund <- abund[tax_df$ASVid]
    if (any(is.na(abund))) {
        stop("Some ASV IDs in the tax table were not found in the OTU table.")
    }

    tax_df$TotalCount <- as.numeric(abund)

    all_manifests <- list()

    for (comparison_name in names(sig_genera_by_comparison)) {
        comp_genera <- sig_genera_by_comparison[[comparison_name]]
        comp_dir <- file.path(out_dir, comparison_name)
        dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

        comp_manifest_list <- list()

        for (genus in comp_genera) {

            dir.create(file.path(comp_dir,genus), recursive = TRUE, showWarnings = FALSE)

            genus_asvs <- tax_df$ASVid[!is.na(tax_df$Genus) & tax_df$Genus == genus]

            if (length(genus_asvs) == 0) {
                warning("No ASVs found for genus '", genus, "' in comparison '", comparison_name, "'.")
                next
            }

            # Subset FASTA
            genus_fasta <- asv_fasta[genus_asvs]

            # Build manifest
            m <- tax_df[match(genus_asvs, tax_df$ASVid), , drop = FALSE]
            m$Comparison <- comparison_name
            m$Sequence <- as.character(genus_fasta)
            m$RelativeAbundanceWithinGenus <- m$TotalCount / sum(m$TotalCount)

            # Reorder columns for readability
            m <- m[, c(
                "Comparison", taxa_level_col,
                "ASVid",
                "TotalCount", "RelativeAbundanceWithinGenus",
                "Sequence",
                setdiff(colnames(m), c(
                "Comparison", taxa_level_col,
                "ASVid",
                "TotalCount", "RelativeAbundanceWithinGenus",
                "Sequence"
                ))
            )]

            manifest_file <- file.path(comp_dir, genus, paste0(genus, "_ASV_manifest.tsv"))
            fasta_file <- file.path(comp_dir, genus, paste0(genus, "_ASV.fasta"))

            fwrite(as.data.table(m), file = manifest_file, sep = "\t", quote = FALSE)
            writeXStringSet(genus_fasta, filepath = fasta_file)

            comp_manifest_list[[genus]] <- m
        }

        if (length(comp_manifest_list) > 0) {
            comp_manifest <- rbindlist(lapply(comp_manifest_list, as.data.table), fill = TRUE)
            comp_manifest_file <- file.path(comp_dir, paste0(comparison_name, "_AllGenera_ASV_manifest.tsv"))
            fwrite(comp_manifest, file = comp_manifest_file, sep = "\t", quote = FALSE)
            all_manifests[[comparison_name]] <- comp_manifest
        }
    }

    invisible(all_manifests)
}

sig_genera_by_comparison <- list()
input_trial_id_num <- sub("_.*", "", args$in_trial)
input_dir <- file.path(args$base_dir, args$in_trial)
output_dir <- file.path(args$base_dir, args$out_trial)

for (comparison in args$DA_comparisons) {
    DA_results_df <- read.delim(file.path(input_dir, args$DA_method, comparison, paste0(input_trial_id_num, "_", comparison, "_", args$DA_method, "Results.tsv")))
    sig_genera <- DA_results_df[DA_results_df$significance == "Sig", "taxon"]
    sig_genera_by_comparison[[comparison]] <- sig_genera
}

load(file.path(input_dir, "CompPhyseq.RData"))

export_significant_genus_asvs(
    ps = physeq,
    asv_fasta = args$all_asv_fasta,
    sig_genera_by_comparison = sig_genera_by_comparison,
    out_dir = output_dir,
    taxa_level_col = "Genus"
)