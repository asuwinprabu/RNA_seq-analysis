#!/usr/bin/env Rscript

.libPaths("~/R/libs")

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

cat("=== Pathway Analysis ===\n\n")

# Load DE results
cat("1. Loading DE results...\n")
de_results <- read.csv("results/differential_expression/deseq2_results_all.csv", stringsAsFactors = FALSE)
cat(sprintf("   Total genes: %d\n", nrow(de_results)))

# Create output directory
dir.create("results/pathway_analysis", recursive = TRUE, showWarnings = FALSE)

# Get significant genes
sig_genes <- de_results %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(gene_id)

cat(sprintf("   Significant genes: %d\n", length(sig_genes)))

# Map gene IDs to Entrez
cat("\n2. Mapping gene IDs...\n")
gene_map <- bitr(sig_genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
cat(sprintf("   Mapped to Entrez: %d\n", nrow(gene_map)))

# GO Biological Process enrichment
cat("\n3. Running GO enrichment...\n")
ego <- enrichGO(gene = gene_map$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "results/pathway_analysis/go_bp_enrichment.csv", row.names = FALSE)
cat(sprintf("   Significant GO terms: %d\n", nrow(ego_df)))

# GO dotplot
if (nrow(ego_df) > 0) {
    pdf("results/pathway_analysis/go_bp_dotplot.pdf", width = 10, height = 10)
    print(dotplot(ego, showCategory = 20, title = "GO Biological Process"))
    dev.off()
    cat("   Saved: go_bp_dotplot.pdf\n")
}

# KEGG pathway enrichment
cat("\n4. Running KEGG enrichment...\n")
ekegg <- enrichKEGG(gene = gene_map$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

ekegg_df <- as.data.frame(ekegg)
write.csv(ekegg_df, "results/pathway_analysis/kegg_enrichment.csv", row.names = FALSE)
cat(sprintf("   Significant KEGG pathways: %d\n", nrow(ekegg_df)))

# KEGG dotplot
if (nrow(ekegg_df) > 0) {
    pdf("results/pathway_analysis/kegg_dotplot.pdf", width = 10, height = 8)
    print(dotplot(ekegg, showCategory = 20, title = "KEGG Pathways"))
    dev.off()
    cat("   Saved: kegg_dotplot.pdf\n")
}

# Create summary
sink("results/pathway_analysis/pathway_summary.txt")
cat("=== Pathway Analysis Summary ===\n\n")
cat(sprintf("Input genes: %d significant DE genes\n", length(sig_genes)))
cat(sprintf("Mapped to Entrez: %d\n", nrow(gene_map)))
cat(sprintf("\nGO Biological Process:\n"))
cat(sprintf("  Significant terms: %d\n", nrow(ego_df)))
if (nrow(ego_df) > 0) {
    cat("\n  Top 10 GO terms:\n")
    top_go <- head(ego_df, 10)
    for (i in 1:nrow(top_go)) {
        cat(sprintf("    - %s (p.adj=%.2e)\n", top_go$Description[i], top_go$p.adjust[i]))
    }
}
cat(sprintf("\nKEGG Pathways:\n"))
cat(sprintf("  Significant pathways: %d\n", nrow(ekegg_df)))
if (nrow(ekegg_df) > 0) {
    cat("\n  Top 10 KEGG pathways:\n")
    top_kegg <- head(ekegg_df, 10)
    for (i in 1:nrow(top_kegg)) {
        cat(sprintf("    - %s (p.adj=%.2e)\n", top_kegg$Description[i], top_kegg$p.adjust[i]))
    }
}
sink()

cat("\n=== Pathway Analysis Complete! ===\n")
cat("Results saved in: results/pathway_analysis/\n")
