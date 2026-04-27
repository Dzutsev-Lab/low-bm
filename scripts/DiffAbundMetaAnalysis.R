library(tibble)
library(dplyr)
library(ggplot2)
library(argparse)
library(stringr)
library(metap)
library(ggrepel)

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
parser$add_argument("--comparison",
                    type = "character",
                    help = "Comparison used for differential abundance (Options: NTtoT, NegativeControltoPS)")
parser$add_argument("--out-dir",
                    type = "character",
                    help = "desired output directory")

args <- parser$parse_args()


#------------------------------------------------
# Stardarized Differential Abundance tsv Read-in
#------------------------------------------------
readin_DA <- function(file, batch_label) {
    df <- read.delim(
        file,
        header = TRUE,
        sep = "\t"
    )

    df |>
        mutate(
            batch = batch_label
        )
}

# Parse Batch Name
B1_name_components <- str_split(args$batch1_Name, "_")[[1]]
B1_ID <- B1_name_components[1]
B1_ExpID <- B1_name_components[2]
B2_name_components <- str_split(args$batch2_Name, "_")[[1]]
B2_ID <- B2_name_components[1]
B2_ExpID <- B2_name_components[2]

#---------------------------------------------------
# Read in Limma Voom Differential Abundance Results
#---------------------------------------------------
B1_DA_results <- readin_DA(file = paste0("Exp_Output/", args$in_dir, "/", 
                                         B1_ID, "_", B1_ExpID, "_", args$comparison, "_ANCOMBCResults.tsv"),
                           batch_label = B1_ExpID) |>
                    rename(
                        b1 = batch,
                        logFC_b1 = log2FoldChange,
                        p_b1 = p,
                        adj_p_b1 = padj,
                        se_b1 = se
                    )
print(str(B1_DA_results))
B2_DA_results <- readin_DA(file = paste0("Exp_Output/", args$in_dir, "/", 
                                         B2_ID, "_", B2_ExpID, "_", args$comparison, "_ANCOMBCResults.tsv"),
                           batch_label = B2_ExpID) |>
                    rename(
                        b2 = batch,
                        logFC_b2 = log2FoldChange,
                        p_b2 = p,
                        adj_p_b2 = padj,
                        se_b2 = se
                    )
str(B2_DA_results)


#------------------------------------------------
# Fisher composite q value calculator
#------------------------------------------------
calc_fisher_q <- function(Composite_DA_df) {
    required_cols <- c("p_b1", "p_b2")
    missing_cols <- setdiff(required_cols, names(Composite_DA_df))
    if (length(missing_cols) > 0) {
    stop("Missing required columns for heterogeneity calculation: ",
            paste(missing_cols, collapse = ", "))
    }

    Composite_DA_df |> rowwise() |>
    mutate(
        fisher_p = sumlog(c(p_b1, p_b2))$p
    ) |> ungroup() |>
    mutate(
        fisher_q = p.adjust(fisher_p, method = "BH"),
        fisher_neglog10_p = -log10(pmax(fisher_p, .Machine$double.xmin)),
        fisher_neglog10_q = -log10(pmax(fisher_q, .Machine$double.xmin))
    )
}

calc_random_effect_heterogeneity <- function(Composite_DA_df) {

    required_cols <- c("logFC_b1", "se_b1", "logFC_b2", "se_b2")
    missing_cols <- setdiff(required_cols, names(Composite_DA_df))
    if (length(missing_cols) > 0) {
    stop("Missing required columns for heterogeneity calculation: ",
            paste(missing_cols, collapse = ", "))
    }

    Composite_DA_df |>
        rowwise() |>
        mutate(
            w_b1 = if (is.finite(se_b1) && se_b1 > 0) 1 / (se_b1^2) else NA_real_,
            w_b2 = if (is.finite(se_b2) && se_b2 > 0) 1 / (se_b2^2) else NA_real_,
            fixed_effect = if (is.finite(w_b1) && is.finite(w_b2) && (w_b1 + w_b2) > 0) {
            (w_b1 * logFC_b1 + w_b2 * logFC_b2) / (w_b1 + w_b2)
            } else {
            NA_real_
            },
            het_Q = if (is.finite(fixed_effect) && is.finite(w_b1) && is.finite(w_b2)) {
            w_b1 * (logFC_b1 - fixed_effect)^2 + w_b2 * (logFC_b2 - fixed_effect)^2
            } else {
            NA_real_
            },
            het_Q_df = if (is.finite(het_Q)) 1L else NA_integer_,
            het_Q_p = if (is.finite(het_Q)) pchisq(het_Q, df = 1, lower.tail = FALSE) else NA_real_
        ) |>
        ungroup() |>
        mutate(
            het_Q_q = p.adjust(het_Q_p, method = "BH"),
            # basing heterogeneity on unadjusted het_Q p-value out of caution
            non_heterogeneous = het_Q_p > 0.2,
            het_neglog10_q = -log10(pmax(het_Q_q, .Machine$double.xmin)),
            heterogeneity_class = ifelse(
                non_heterogeneous,
                "Low heterogeneity",
                "High heterogeneity")
        )
}

Comp_DA_results <- inner_join(B1_DA_results, B2_DA_results, by = "taxon")
Comp_DA_results <- calc_fisher_q(Composite_DA_df= Comp_DA_results)
Comp_DA_results <- calc_random_effect_heterogeneity(Composite_DA_df= Comp_DA_results)

spearman_res <- cor.test(
  Comp_DA_results$logFC_b1,
  Comp_DA_results$logFC_b2,
  method = "spearman",
  exact = FALSE
)

LFCxLFC_plot <- ggplot(Comp_DA_results, aes(x = logFC_b1, y = logFC_b2)) +
    geom_point(aes(color = fisher_neglog10_q, shape = heterogeneity_class)) +
    geom_text_repel(
        data = Comp_DA_results |>
            filter(
                heterogeneity_class == "Low heterogeneity",
                fisher_q <= 0.05
            ),
        aes(label = taxon),
        size = 3.5,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.3,
        segment.size = 0.4,
        segment.color = "grey50",
        min.segment.length = 0
    ) +
    scale_shape_manual(
        values = c(
            "Low heterogeneity" = 16,   # ● filled circle
            "High heterogeneity" = 4    # × cross
        ),
        name = "Random Effect Heterogeneity"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    scale_color_gradient(
        low = "grey80",
        high = "firebrick",
        name = "-log10(Fisher q, BH adjusted)"
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
    filename = file.path("Exp_Output", args$out_dir, "batch_logFC_concordance.png"),
    plot = LFCxLFC_plot,
    width = 7,
    height = 6,
    dpi = 300
)

write.table(
    Comp_DA_results,
    file = file.path("Exp_Output", args$out_dir, "shared_taxa_meta.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)