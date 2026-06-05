library(phyloseq)
library(Biostrings)
library(data.table)
library(tibble)
library(stringr)
library(argparse)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    default = NULL,
                    help = "Analysis YAML with differential_abundance comparison specs")

parser$add_argument("--DA-comparisons",
                    type = "character",
                    nargs = "+",
                    help = "Choice of differential abunance comprisons to pull significant taxa from [CellLineControltoTumor, CellLineControltoNontumor, NegativeControl, PatientSample, TumorType]")
parser$add_argument("--DA-method",
                    type = "character",
                    help = "Differential abundance method used to perform comparisons")
parser$add_argument("--compiled-asv-fasta",
                    type = "character",
                    help = "FASTA file with all ASV sequences (ASV IDs as headers)")
parser$add_argument("--io-dir",
                    type = "character",
                    help = "Directory to store output manifests and select ASV fasta's within base directory")
parser$add_argument("--trialID",
                    type = "character")
parser$add_argument("--taxa-level",
                    type = "character")
parser$add_argument("--compiled-physeq",
                    type = "character")
parser$add_argument("--base-dir",
                    type = "character",
                    default = "Exp_Output",
                    help = "Base directory containing the trial directory")

args <- parser$parse_args()

project_config <- list()
ordination_config <- list()
if (!is.null(args$analysis_config)) {
  cfg <- load_yaml_config(args$analysis_config)
  project_config <- cfg$project
  blast_config <- cfg$blast_confirmation
}

load_input_physeq <- function() {
    if (!is.null(args$compiled_physeq)) {
        return(load_physeq(args$compiled_physeq))
    }
    if (!is.null(project_config$compiled_physeq) && file.exists(project_config$compiled_physeq)) {
        return(load_physeq(project_config$compiled_physeq))
    }
    if (!is.null(project_config$batch_table)) {
        batch_table <- project_config$batch_table
        physeq_paths <- resolve_batch_physeqs(batch_table, base_dir = base_dir)
        return(merge_physeqs(load_physeqs(physeq_paths)))
    }
    stop("Provide --compiled-physeq script argument or project.compiled_physeq or project.batch_table in --analysis-config.", call. = FALSE)
}

#-----------------------------
# Main export function
#-----------------------------
export_significant_genus_asvs <- function(ps,
                                          asv_fasta,
                                          sig_genera_by_comparison,
                                          out_dir,
                                          taxa_level_col = "Genus") {

    # Read All ASV FASTA
    asv_fasta <- readDNAStringSet(asv_fasta)

    # Extract taxonomy (first as matrix then convert to data frame)
    tax_df <- as.data.frame(as(tax_table(ps), "matrix")) |>
        rownames_to_column("ASVid")

    if (!taxa_level_col %in% colnames(tax_df)) {
        stop("taxa_level_col = '", taxa_level_col, "' not found in tax_table(ps).")
    }

    # Extract abundance
    abund <- otu_totals_by_taxa(ps)
    abund <- abund[tax_df$ASVid]
    if (any(is.na(abund))) {
        stop("Some ASV IDs in the tax table were not found in the OTU table.")
    }

    tax_df$TotalCount <- as.numeric(abund)

    all_manifests <- list()

    for (comparison_name in names(sig_genera_by_comparison)) {
        comp_genera <- sig_genera_by_comparison[[comparison_name]]
        comp_dir <- file.path(out_dir, "BlastAnalysis", comparison_name)
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

CompPhyseq <- load_input_physeq()

if (!is.null(args$analysis_config)) {
    trialID <- blast_config$trialID
    io_dir <- blast_config$io_dir
    DA_comparisons <- blast_config$DA_comparisons
    DA_method <- blast_config$DA_method
    taxa_level <- blast_config$taxa_level
    compiled_asv_fasta <- project_config$compiled_asv_fasta
} else {
    trialID <- args$trialID
    io_dir <- args$io_dir
    DA_comparisons <- args$DA_comparisons
    DA_method <- args$DA_method
    taxa_level <- args$taxa_level
    compiled_asv_fasta <- args$compiled_asv_fasta
}

#-------------------------
# Input Check
#-------------------------
required_input <- list(trialID, io_dir, DA_comparisons, DA_method, taxa_level, compiled_asv_fasta)
names(required_input) <- c("--trialID", "--io-dir", "--DA-comparisons", "--DA-method", "--taxa-level", "--compiled-asv-fasta")
missing_input <- names(required_input)[sapply(required_input, is.null)]
if (length(missing_input) > 0) {
    stop(sprintf("Missing required input. Provide: %s as script argument(s) or in --analysis-config.", paste(unlist(missing_input), collapse = ", ")))
}




sig_genera_by_comparison <- list()
for (comparison in DA_comparisons) {
    DA_results_df <- read.delim(file.path(io_dir, DA_method, comparison, 
                                          paste0(trialID, "_", comparison, "_", DA_method, "Results.tsv")))
    sig_genera <- DA_results_df[DA_results_df$significance == "Sig", "taxon"]
    sig_genera_by_comparison[[comparison]] <- sig_genera
}

export_significant_genus_asvs(
    ps = CompPhyseq,
    asv_fasta = compiled_asv_fasta,
    sig_genera_by_comparison = sig_genera_by_comparison,
    out_dir = io_dir,
    taxa_level_col = taxa_level
)