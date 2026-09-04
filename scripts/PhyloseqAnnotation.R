library(argparse)
library(phyloseq)

source(file.path("scripts", "Rhelpers", "PhyloseqIO.R"))
source(file.path("scripts", "Rhelpers", "PhyloseqAnnotation.R"))

parser <- ArgumentParser()

parser$add_argument("--analysis-config",
                    type = "character",
                    required = TRUE,
                    help = "Meta-analysis YAML with a meta_annotation section")

args <- parser$parse_args()

cfg <- load_yaml_config(args$analysis_config)
project_config <- cfg$project %||% list()
annotation_config <- cfg$meta_annotation %||% list()

if (length(annotation_config) == 0) {
  stop("Config is missing required section: meta_annotation", call. = FALSE)
}

input_physeq <- annotation_scalar(
  config_value(annotation_config, "input_physeq") %||%
    config_value(project_config, "compiled_physeq")
)
if (is.null(input_physeq)) {
  stop("Set meta_annotation.input_physeq to one phyloseq endpoint.", call. = FALSE)
}

annotation_table <- annotation_scalar(config_value(annotation_config, "annotation_table"))
if (is.null(annotation_table)) {
  stop("Set meta_annotation.annotation_table to a TSV annotation file.", call. = FALSE)
}
if (!file.exists(annotation_table)) {
  stop("Missing annotation table: ", annotation_table, call. = FALSE)
}

level <- config_value(annotation_config, "level")
annotation_key <- config_value(annotation_config, "annotation_key")
columns <- config_value(annotation_config, "columns")
prefix <- config_value(annotation_config, "prefix")
overwrite <- config_value(annotation_config, "overwrite") %||% "error"

out_dir <- analysis_output_dir(
  project_config,
  annotation_config,
  default = dirname(input_physeq)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

resolve_annotation_output_path <- function(path, out_dir) {
  path <- annotation_scalar(path)
  if (is.null(path)) {
    return(NULL)
  }
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(out_dir, path)
}

output_physeq <- resolve_annotation_output_path(
  config_value(annotation_config, "output_physeq") %||% "AnnotatedPhyseq.RData",
  out_dir
)
output_report <- resolve_annotation_output_path(
  config_value(annotation_config, "output_report") %||% "PhyseqAnnotationReport.tsv",
  out_dir
)

annotation_df <- read.delim(
  annotation_table,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

physeq <- load_physeq(input_physeq)
annotation_result <- annotate_physeq_sample_data(
  physeq = physeq,
  annotation_df = annotation_df,
  level = level,
  annotation_key = annotation_key,
  columns = columns,
  prefix = prefix,
  overwrite = overwrite
)

save_physeq(annotation_result$physeq, output_physeq)
dir.create(dirname(output_report), recursive = TRUE, showWarnings = FALSE)
write.table(
  annotation_result$audit,
  file = output_report,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Annotated phyloseq object saved to: ", output_physeq)
message("Phyloseq annotation report saved to: ", output_report)
