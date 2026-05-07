library(tibble)
library(dplyr)
library(ggplot2)
library(argparse)
library(stringr)
library(metap)
library(ggrepel)
library(patchwork)
library(stringr)

parser <- ArgumentParser()

parser$add_argument("--batch1-Name",
                    type = "character",
                    help = "Name for batch 1 to locate differential abundance results within Exp_Output")
parser$add_argument("--batch2-Name",
                    type = "character",
                    help = "Name for batch 2 to locate differential abundance results within Exp_Output")
parser$add_argument("--DA-method",
                    type = "character",
                    help = "The tool used to generate differenital abundance results (Options: LIMMA_VOOM, ANCOMBC)")
parser$add_argument("--comparison",
                    type = "character",
                    help = "Comparison used for differential abundance (Options: CellLineControltoTumor, CellLineControltoNontumor, NegativeControl, PatientSample)")
parser$add_argument("--hetQ-p-cutoff",
                    type = "double",
                    help = "Threshold value for random effect model hetergenity Q p-value to differentiate between low- and high-heterogeneity taxa")
parser$add_argument("--fisher-q-cutoff",
                    type = "double",
                    help = "Threshold value for adjusted fisher composite p-values to differentiate significant differentially abundant taxa")
parser$add_argument("--out-trial",
                    type = "character",
                    help = "desired output directory within Exp_Ouput")

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
B1_DA_results <- readin_DA(file = paste0("Exp_Output/", args$batch1_Name, "/", 
                                         args$DA_method, "/",
                                         args$comparison, "/",
                                         B1_ID, "_", args$comparison, "_", args$DA_method, "Results.tsv"),
                           batch_label = B1_ExpID) |>
                    rename(
                        b1 = batch,
                        logFC_b1 = log2FoldChange,
                        p_b1 = p,
                        adj_p_b1 = padj,
                        se_b1 = se,
                        direction_b1 = direction
                    )
B2_DA_results <- readin_DA(file = paste0("Exp_Output/", args$batch2_Name, "/", 
                                         args$DA_method, "/",
                                         args$comparison, "/",
                                         B2_ID, "_", args$comparison, "_", args$DA_method, "Results.tsv"),
                           batch_label = B2_ExpID) |>
                    rename(
                        b2 = batch,
                        logFC_b2 = log2FoldChange,
                        p_b2 = p,
                        adj_p_b2 = padj,
                        se_b2 = se,
                        direction_b2 = direction
                    )

if (args$DA_method == "ANCOMBC") {
    B1_DA_results <- B1_DA_results |>
                        rename(struc0_b1 = struc0)
    B2_DA_results <- B2_DA_results |>
                        rename(struc0_b2 = struc0)
}

#------------------------------------------------
# Fisher composite q value calculator
#------------------------------------------------
calc_fisher_q <- function(Composite_DA_df,
                          fisher_q_cutoff) {
    required_cols <- c("p_b1", "p_b2")
    missing_cols <- setdiff(required_cols, names(Composite_DA_df))
    if (length(missing_cols) > 0) {
        stop("Missing required columns for heterogeneity calculation: ",
             paste(missing_cols, collapse = ", "))
    }

    Composite_DA_df |> 
    rowwise() |>
    mutate(
        fisher_p = sumlog(c(p_b1, p_b2))$p
    ) |> ungroup() |>
    mutate(
        fisher_q = p.adjust(fisher_p, method = "BH"),
        fisher_neglog10_p = -log10(pmax(fisher_p, .Machine$double.xmin)),
        fisher_neglog10_q = -log10(pmax(fisher_q, .Machine$double.xmin)),
        fisher_sig = fisher_q <= fisher_q_cutoff
    )
}

calc_random_effect_heterogeneity <- function(Composite_DA_df,
                                             hetQ_p_cutoff) {

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
            non_heterogeneous = het_Q_p > hetQ_p_cutoff,
            het_neglog10_q = -log10(pmax(het_Q_q, .Machine$double.xmin)),
            heterogeneity_class = ifelse(
                non_heterogeneous,
                "Low heterogeneity",
                "High heterogeneity")
        )
}

Comp_DA_results <- inner_join(B1_DA_results, B2_DA_results, by = "taxon")
# Extract Agreed Structural Zeros
agreed_struc0_df <- Comp_DA_results |>
    filter(!is.na(struc0_b1) & !is.na(struc0_b2))
# Filter out any taxa with any Structural Zero (disconcordant or concordant)
Comp_DA_results <- Comp_DA_results |> 
    filter(is.na(struc0_b1) & is.na(struc0_b2))
Comp_DA_results <- calc_fisher_q(Composite_DA_df = Comp_DA_results,
                                 fisher_q_cutoff = args$fisher_q_cutoff)
Comp_DA_results <- calc_random_effect_heterogeneity(Composite_DA_df= Comp_DA_results,
                                                    hetQ_p_cutoff = args$hetQ_p_cutoff)

# Create subset data frame with only taxa that meet our desired thresolds for concordance and signficance
Concord_DA_results <- subset(Comp_DA_results,
                                  heterogeneity_class == "Low heterogeneity" &
                                  fisher_sig)

Taxa_of_Interest <- Concord_DA_results$taxon

spearman_res <- cor.test(
  Comp_DA_results$logFC_b1,
  Comp_DA_results$logFC_b2,
  method = "spearman",
  exact = FALSE
)

ConcordancePlot <- ggplot(Comp_DA_results, aes(x = logFC_b1, y = logFC_b2)) +
    geom_point(
        aes(color = fisher_neglog10_q, 
            shape = heterogeneity_class)
    ) +
    geom_text_repel(
        data = Concord_DA_results,
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
        title = paste(args$comparison, "DA concordance"),
        subtitle = paste(B1_ExpID, "vs", B2_ExpID),
        caption = sprintf(
            "Shared taxa = %d; Spearman rho = %.2f (p = %.2g)",
            nrow(Comp_DA_results),
            unname(spearman_res$estimate),
            spearman_res$p.value
        ),
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    theme_classic(base_size = 12)

if (nrow(agreed_struc0_df) > 0) {
    Struc0Plot <- ggplot(agreed_struc0_df, aes(x = 0, y = factor(taxon))) +
    geom_text(aes(label = paste0(taxon, " (", direction_b1, ")")),
                hjust = 0, 
                size = 3) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_void() +
    labs(title = "Structural Zeros") +
    theme(plot.margin = margin(5.5, 20, 5.5, 5.5))

    ConcordancePlot <- (ConcordancePlot / Struc0Plot) + plot_layout(heights = c(3,0.5))
}



if (!dir.exists(paste0("Exp_Output/", args$out_trial, "/", args$DA_method))) {
    dir.create(paste0("Exp_Output/", args$out_trial, "/", args$DA_method),
    recursive = TRUE)
}

ggsave(
    filename = paste0("Exp_Output/", args$out_trial, "/",
                      args$DA_method, "/",
                      B1_ExpID, "v", B2_ExpID, "_", args$comparison, "_logFC_concordance.png"),
    plot = ConcordancePlot,
    width = 7,
    height = 6,
    dpi = 300
)

write.table(
    Comp_DA_results,
    file = paste0("Exp_Output/", args$out_trial, "/",
                  args$DA_method, "/",
                  B1_ExpID, "v", B2_ExpID, "_", args$comparison, "_shared_taxa_meta.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(Taxa_of_Interest, 
          file = paste0("Exp_Output/", args$out_trial, "/",
                        args$DA_method, "/",
                        B1_ExpID, "v", B2_ExpID, "_", args$comparison, "_taxa_of_interest.csv"),
          sep = ",",
          col.names = FALSE,
          row.names = FALSE)