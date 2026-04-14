library(tibble)
library(dplyr)
library(ggplot2)

Batch1Results <- read.delim("./041326.1_ItalyLungExp1_RerunMicRocleanTrial/ANCOMBC/NTtoT/041326.1_NTtoT_ANCOMBCResults.tsv", sep = "\t", header = TRUE)
Batch2Results <- read.delim("./041026.1_ItalyLungExp2_micRocleanTrial/ANCOMBC/NTtoT/041026.1_NTtoT_ANCOMBCResults.tsv", sep = "\t", header = TRUE)

Batch1Results <- Batch1Results |> column_to_rownames(var = "taxon")
Batch2Results <- Batch2Results |> column_to_rownames(var = "taxon")

common_taxa <- intersect(rownames(Batch1Results), rownames(Batch2Results))
BatchCompResults <- data.frame(
  row.names = common_taxa,
  "LFC_B1" = Batch1Results[common_taxa, "log2FoldChange"],
  "LFC_B2" = Batch2Results[common_taxa, "log2FoldChange"],
  "padj_B1" = Batch1Results[common_taxa, "padj"],
  "padj_B2" = Batch2Results[common_taxa, "padj"]
)


# LFC vs LFC Plotting
ggplot(BatchCompResults, aes(x = LFC_B1, y = LFC_B2)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Log Fold Change Comparison Between Sequencing Batches",
    x = "LFC (Batch 1)",
    y = "LFC (Batch 2)"
  )

# Correlation Calculations
cor_pearson <- cor(BatchCompResults$LFC_B1,
                   BatchCompResults$LFC_B2,
                   method = "pearson")
cor_spearman <- cor(BatchCompResults$LFC_B1,
                   BatchCompResults$LFC_B2,
                   method = "spearman")

# Concordance of Direction Calculations
prop_same_direction <- sum(
  sign(BatchCompResults$LFC_B1) == sign(BatchCompResults$LFC_B2)
) / nrow(BatchCompResults)
