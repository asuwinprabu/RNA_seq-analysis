#!/usr/bin/env Rscript

.libPaths("~/R/libs")

library(ggplot2)
library(dplyr)

cat("=== Pathway Analysis (Simplified) ===\n\n")

# Load DE results
cat("1. Loading DE results...\n")
de_results <- read.csv("results/differential_expression/deseq2_results_all.csv", stringsAsFactors = FALSE)

dir.create("results/pathway_analysis", recursive = TRUE, showWarnings = FALSE)

# Get significant genes
sig_up <- de_results %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
    arrange(padj)

sig_down <- de_results %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < -1) %>%
    arrange(padj)

cat(sprintf("   Upregulated: %d\n", nrow(sig_up)))
cat(sprintf("   Downregulated: %d\n", nrow(sig_down)))

# Save gene lists for external pathway analysis
cat("\n2. Saving gene lists...\n")
write.csv(sig_up, "results/pathway_analysis/upregulated_genes.csv", row.names = FALSE)
write.csv(sig_down, "results/pathway_analysis/downregulated_genes.csv", row.names = FALSE)

# Save just gene names for online tools (DAVID, Enrichr, etc.)
writeLines(sig_up$gene_id, "results/pathway_analysis/upregulated_gene_list.txt")
writeLines(sig_down$gene_id, "results/pathway_analysis/downregulated_gene_list.txt")
cat("   Gene lists saved for external analysis\n")

# Create bar plot of top genes
cat("\n3. Creating plots...\n")
top_genes <- rbind(head(sig_up, 15), head(sig_down, 15))
top_genes$direction <- ifelse(top_genes$log2FoldChange > 0, "Up", "Down")

p <- ggplot(top_genes, aes(x = reorder(gene_id, log2FoldChange), y = log2FoldChange, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    theme_minimal(base_size = 12) +
    labs(x = "Gene", y = "log2 Fold Change", 
         title = "Top 30 Differentially Expressed Genes",
         subtitle = "Tumor vs Normal")

ggsave("results/pathway_analysis/top_genes_barplot.pdf", p, width = 10, height = 10)
cat("   Saved: top_genes_barplot.pdf\n")

# Summary
sink("results/pathway_analysis/pathway_summary.txt")
cat("=== Pathway Analysis Summary ===\n\n")
cat(sprintf("Upregulated genes (log2FC > 1, FDR < 0.05): %d\n", nrow(sig_up)))
cat(sprintf("Downregulated genes (log2FC < -1, FDR < 0.05): %d\n", nrow(sig_down)))

cat("\n=== TOP 20 UPREGULATED ===\n")
print(head(sig_up[, c("gene_id", "log2FoldChange", "padj")], 20), row.names = FALSE)

cat("\n=== TOP 20 DOWNREGULATED ===\n")
print(head(sig_down[, c("gene_id", "log2FoldChange", "padj")], 20), row.names = FALSE)

cat("\n\n=== FOR FULL PATHWAY ANALYSIS ===\n")
cat("Upload gene lists to these free online tools:\n")
cat("  - Enrichr: https://maayanlab.cloud/Enrichr/\n")
cat("  - DAVID: https://david.ncifcrf.gov/\n")
cat("  - g:Profiler: https://biit.cs.ut.ee/gprofiler/\n")
cat("\nGene list files:\n")
cat("  - results/pathway_analysis/upregulated_gene_list.txt\n")
cat("  - results/pathway_analysis/downregulated_gene_list.txt\n")
sink()

cat("\n=== Pathway Analysis Complete! ===\n")
cat("Results saved in: results/pathway_analysis/\n")
cat("\nTip: Upload gene lists to Enrichr for full pathway analysis\n")
