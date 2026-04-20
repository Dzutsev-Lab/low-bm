library(tibble)
library(dplyr)
library(ggplot2)
library(argparse)
library(stringr)
library(metap)

parser <- ArgumentParser()

parser$add_argument("--in-dir",
                    type = "character",
                    help = "Directory where limma voom DA results can be found")
parser$add_argument("--batch1-Name",
                    type = "character",
                    help = "Name for batch 1 to locate limma voom results")
parser$add_argument("--batch2-Name",
                    type = "character",
                    help = "Name for batch 2 to locate limma voom results")
parser$add_argument("--norm-method",
                    type = "character",
                    help = "normalization method used prior to limma voom")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "desired output directory")

args <- parser$parse_args()

readin_DA <- function(file, batch_label) {
    df <- read.delim(
        file,
        header = TRUE,
        sep = "\t"
    )

    names(df)[1] <- "taxon"

    df |>
        mutate(
            batch = batch_label
        )
}

B1_name_components <- str_split(args$batch1_Name, "_")[[1]]
B1_ID <- B1_name_components[1]
B1_ExpID <- B1_name_components[2]
B2_name_components <- str_split(args$batch2_Name, "_")[[1]]
B2_ID <- B2_name_components[1]
B2_ExpID <- B2_name_components[2]


B1_DA_results <- readin_DA(file = paste0("Exp_Output/", args$in_dir, "/", 
                                         B1_ID, "_", B1_ExpID, "_", args$norm_method, "_limmavoom_results.tsv"),
                           batch_label = B1_ExpID) |>
                    rename(
                        logFC_b1 = logFC,
                        AveExpr_b1 = AveExpr,
                        t_b1 = t,
                        p_b1 = P.Value,
                        adj_p_b1 = adj.P.Val,
                        B_b1 = B
                    )
B2_DA_results <- readin_DA(file = paste0("Exp_Output/", args$in_dir, "/", 
                                         B2_ID, "_", B2_ExpID, "_", args$norm_method, "_limmavoom_results.tsv"),
                           batch_label = B2_ExpID) |>
                    rename(
                        logFC_b2 = logFC,
                        AveExpr_b2 = AveExpr,
                        t_b2 = t,
                        p_b2 = P.Value,
                        adj_p_b2 = adj.P.Val,
                        B_b2 = B
                    )

Comp_DA_results <- inner_join(
                        B1_DA_results |> select(taxon, logFC_b1, p_b1, adj_p_b1),
                        B2_DA_results |> select(taxon, logFC_b2, p_b2, adj_p_b2),
                        by = "taxon"
    ) |> rowwise() |>
    mutate(
        fisher_p = sumlog(c(p_b1, p_b2))$p
    ) |> ungroup() |>
    mutate(
        fisher_q = p.adjust(fisher_p, method = "BH"),
        # TODO: determine what this is doing
        fisher_neglog10_q = -log10(pmax(fisher_q, .Machine$double.xmin))
    )

spearman_res <- cor.test(
  Comp_DA_results$logFC_b1,
  Comp_DA_results$logFC_b2,
  method = "spearman",
  exact = FALSE
)

LFCxLFC_plot <- ggplot(Comp_DA_results, aes(x = logFC_b1, y = logFC_b2, color = fisher_neglog10_q)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  scale_color_gradient(
    low = "grey80",
    high = "firebrick",
    name = "-log10(Fisher BH q)"
  ) +
  labs(
    x = paste0(B1_ExpID, " logFC"),
    y = paste0(B2_ExpID, " logFC"),
    title = paste("DA concordance between", B1_ExpID, "and", B2_ExpID),
    subtitle = sprintf(
      "Shared taxa = %d; Spearman rho = %.2f (p = %.2g)",
      nrow(Comp_DA_results),
      unname(spearman_res$estimate),
      spearman_res$p.value
    ),
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  theme_classic(base_size = 12)

ggsave(
    filename = file.path("Exp_Output", args$out_dir, "batch_logFC_concordance_fisher.png"),
    plot = LFCxLFC_plot,
    width = 7,
    height = 6,
    dpi = 300
)

write.table(
    Comp_DA_results,
    file = file.path("Exp_Output", args$out_dir, "shared_taxa_fisher_meta.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)