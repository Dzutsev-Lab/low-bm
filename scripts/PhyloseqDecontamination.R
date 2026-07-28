library(argparse)
library(micRoclean)
library(phyloseq)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqTransforms.R"))
source(file.path("scripts", "Rhelpers", "MicRocleanDecontamination.R"))

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    required = TRUE,
                    help = "Meta-analysis YAML with a meta_decontamination section")

args <- parser$parse_args()

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
decontam_config <- cfg$meta_decontamination %||% list()

if (length(decontam_config) == 0) {
  stop("Config is missing required section: meta_decontamination", call. = FALSE)
}

input_physeq <- config_value(decontam_config, "input_physeq") %||%
  config_value(project_config, "compiled_physeq")
if (is.null(input_physeq) || is.na(input_physeq) || !nzchar(trimws(as.character(input_physeq)))) {
  stop("Set meta_decontamination.input_physeq to one phyloseq endpoint.", call. = FALSE)
}
input_physeq <- as.character(input_physeq)

out_dir <- analysis_output_dir(
  project_config,
  decontam_config,
  default = dirname(input_physeq)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_physeq <- config_value(decontam_config, "output_physeq") %||% "DecontamPhyseq.RData"
output_filter_report <- config_value(decontam_config, "output_filter_report") %||% "DecontamFilterReport.tsv"
output_decontaminated_names <- config_value(decontam_config, "output_decontaminated_names") %||%
  "decontaminated.ASV.names"
output_contaminant_names <- config_value(decontam_config, "output_contaminant_names") %||%
  "contaminant.ASV.names"

sample_type_column <- config_value(decontam_config, "sample_type_column") %||% "SampleType"
control_sample_types <- config_value(decontam_config, "control_sample_types") %||%
  config_value(decontam_config, "control_sample_type") %||%
  "NegativeControl"
batch_column <- config_value(decontam_config, "batch_column") %||% "ProcessingBatch"
research_goal <- config_value(decontam_config, "research_goal") %||% "biomarker"
control_name <- config_value(decontam_config, "control_name") %||% "Control"
blocklist <- config_value(decontam_config, "blocklist") %||% "FakeBlockedTaxa"

physeq <- load_physeq(input_physeq)
original_taxa <- phyloseq::taxa_names(physeq)
counts_df <- otu_samples_by_taxa(physeq)
microclean_meta <- build_microclean_metadata(
  physeq = physeq,
  sample_type_column = sample_type_column,
  control_sample_types = control_sample_types,
  batch_column = batch_column
)
counts_df <- counts_df[rownames(microclean_meta), , drop = FALSE]

if (any(abs(counts_df - round(counts_df)) > .Machine$double.eps^0.5, na.rm = TRUE)) {
  warning(
    "Input phyloseq counts are not all integer-like. micRoclean is intended for raw ASV counts.",
    call. = FALSE
  )
}

technical_replicates <- data.frame(
  Batch_1 = character(),
  Batch_2 = character()
)

biomarkerID_results <- micRoclean(
  counts = counts_df,
  meta = microclean_meta,
  research_goal = research_goal,
  control_name = control_name,
  blocklist = as.character(unlist(blocklist, use.names = FALSE)),
  technical_replicates = technical_replicates
)

cleaned_physeq <- apply_decontaminated_taxa(
  physeq = physeq,
  decontaminated_count = biomarkerID_results$decontaminated_count
)

physeq <- cleaned_physeq
out_file <- file.path(out_dir, output_physeq)
save_physeq(physeq, out_file)

decontaminated_names <- phyloseq::taxa_names(cleaned_physeq)
contaminant_names <- setdiff(original_taxa, decontaminated_names)

writeLines(decontaminated_names, con = file.path(out_dir, output_decontaminated_names))
writeLines(contaminant_names, con = file.path(out_dir, output_contaminant_names))
write.table(
  biomarkerID_results$contaminant_id,
  file = file.path(out_dir, output_filter_report),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

message("Decontaminated phyloseq object saved to: ", out_file)
