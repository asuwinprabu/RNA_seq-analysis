#!/usr/bin/env Rscript

.libPaths("~/R/libs")

library(ggplot2)
library(dplyr)
library(DESeq2)
library(pheatmap)

cat("=== Differential Expression Analysis ===\n\n")

# Load filtered DESeq2 object from QC
cat("1. Loading data...\n")
dds <- readRDS("results/qc/dds_filtered.rds")
cat(sprintf("   Samples: %d, Genes: %d\n", ncol(dds), nrow(dds)))

# Create output directory
dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Run DESeq2
cat("\n2. Running DESeq2...\n")
dds <- DESeq(dds, fitType = "local")
cat("   Done!\n")

# Get results
cat("\n3. Extracting results...\n")
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)

# Convert to dataframe
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# Count significant genes
n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)

cat(sprintf("   Upregulated (FDR<0.05, log2FC>1): %d\n", n_up))
cat(sprintf("   Downregulated (FDR<0.05, log2FC<-1): %d\n", n_down))

# Save results
cat("\n4. Saving results...\n")
write.csv(res_df, "results/differential_expression/deseq2_results_all.csv", row.names = FALSE)

sig_res <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
write.csv(sig_res, "results/differential_expression/deseq2_results_significant.csv", row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "results/differential_expression/normalized_counts.csv")

saveRDS(dds, "results/differential_expression/dds_analyzed.rds")

# Volcano plot
cat("\n5. Creating volcano plot...\n")
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                             ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "NS")

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray50")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
    theme_minimal(base_size = 14) +
    labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)",
         title = "Volcano Plot: Tumor vs Normal",
         subtitle = sprintf("Up: %d, Down: %d", n_up, n_down))

ggsave("results/differential_expression/volcano_plot.pdf", p_volcano, width = 10, height = 8)
cat("   Saved: volcano_plot.pdf\n")

# Heatmap of top 50 genes
cat("\n6. Creating heatmap...\n")
vsd <- readRDS("results/qc/vsd.rds")
top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], 50)

mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

metadata <- as.data.frame(colData(dds))
annotation_col <- data.frame(Condition = metadata$condition, row.names = rownames(metadata))

pdf("results/differential_expression/heatmap_top50.pdf", width = 10, height = 12)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         fontsize_row = 8,
         main = "Top 50 Differentially Expressed Genes")
dev.off()
cat("   Saved: heatmap_top50.pdf\n")

# Summary
sink("results/differential_expression/de_summary.txt")
cat("=== Differential Expression Summary ===\n\n")
cat(sprintf("Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("Significant (FDR < 0.05, |log2FC| > 1): %d\n", n_up + n_down))
cat(sprintf("  - Upregulated: %d\n", n_up))
cat(sprintf("  - Downregulated: %d\n", n_down))
cat("\nTop 10 Upregulated:\n")
top_up <- head(res_df[res_df$log2FoldChange > 0 & !is.na(res_df$padj), ], 10)
print(top_up[, c("gene_id", "log2FoldChange", "padj")])
cat("\nTop 10 Downregulated:\n")
top_down <- head(res_df[res_df$log2FoldChange < 0 & !is.na(res_df$padj), ], 10)
print(top_down[, c("gene_id", "log2FoldChange", "padj")])
sink()

cat("\n=== Differential Expression Complete! ===\n")
cat("Results saved in: results/differential_expression/\n")
