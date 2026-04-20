library(limma)
library(edgeR)
library(phyloseq)
library(tibble)
library(dplyr)
library(stringr)

library(argparse)

parser <- ArgumentParser()

parser$add_argument("--trialName",
                    type = "character",
                    help = "Trial ID to locate output files (full name: numeber_exp_descript)")
parser$add_argument("--norm-method",
                    type = "character",
                    help = "Determine which sequence table is used")
parser$add_argument("--out",
                    type = "character",
                    help = "directory to store differential abundance output")

args <- parser$parse_args()

trial_name_components <- str_split(args$trialName, "_")[[1]]
trialID <- trial_name_components[1]
trialExp <- trial_name_components[2]

load(paste0("Exp_Output/", args$trialName, "/", trialID, "_raw_kraken_phyloseq.RData"))

physeq <- raw_kraken_phyloseq
subset_physeq <- subset_samples(physeq, SampleType %in% c("Tumor", "NormalTissue"))
metadata <- as.data.frame(as.matrix(sample_data(subset_physeq)))
NTvT_samples <- sample_names(subset_physeq)

seqtable <- read.delim(paste0("Exp_Output/", args$trialName, "/CountNormalization/", trialID, "_", args$norm_method, "_SeqTable.tsv"),
                       header = TRUE,
                       sep = "\t",
                       row.names = 1)
seqtable <- seqtable[rownames(seqtable) %in% NTvT_samples, , drop=FALSE]
tseqtable <- t(seqtable)


counts <- as.matrix(tseqtable)
design <- model.matrix(~ SampleType, 
                       data = metadata)

# TODO: determine if should do prior TMM or some other additional normalization
# would then use dge in place of counts for voom and limma fit
# dge <- DGEList(counts=counts)
# dge <- calcNormFactors(dge)
png(paste0("Exp_Output/", args$out, "/", trialID, "_", trialExp, "_", args$norm_method, "_VoomVariance.png"), width = 800, height = 600)
voom <- voom(counts, design, plot=TRUE)
dev.off()

# TODO: consider using the between-array normalization methods used for microarrays
# suggested for very noisy data by limma package docs (pg 71, sect 15.1)
# voom <- voom(counts, design, plot=TRUE, normlaize="quantile")
# dev.copy(png, paste0(args$out, "/VoomVariance.png"))
# dev.off()

fit <- lmFit(voom, design)
fit <- eBayes(fit)

results <- topTable(fit,
                    coef= ncol(design),
                    number= nrow(counts))

write.table(results, 
            file = paste0("Exp_Output/", args$out, "/", trialID, "_", trialExp, "_", args$norm_method, "_limmavoom_results.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)