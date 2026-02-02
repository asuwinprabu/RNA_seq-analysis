#!/usr/bin/env Rscript

# Set library path
.libPaths("~/R/libs")

library(ggplot2)
library(dplyr)
library(DESeq2)
library(pheatmap)

cat("=== Quality Control Analysis ===\n\n")

# Load data
cat("1. Loading data...\n")
counts <- read.csv("data/counts/raw_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("config/samplesheet_counts.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id

cat(sprintf("   Genes: %d, Samples: %d\n", nrow(counts), ncol(counts)))
cat(sprintf("   Tumor: %d, Normal: %d\n", 
            sum(metadata$condition == "Tumor"),
            sum(metadata$condition == "Normal")))

# Create output directory
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Basic stats
cat("\n2. Computing statistics...\n")
lib_sizes <- colSums(counts)
genes_detected <- colSums(counts > 0)

cat(sprintf("   Library size range: %.1f - %.1f million\n", 
            min(lib_sizes)/1e6, max(lib_sizes)/1e6))

# Create DESeq2 object
cat("\n3. Creating DESeq2 object...\n")
metadata$condition <- factor(metadata$condition, levels = c("Normal", "Tumor"))
counts <- counts[, rownames(metadata)]

dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~ condition
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat(sprintf("   Genes after filtering: %d\n", nrow(dds)))

# VST transformation
cat("\n4. Generating plots...\n")
vsd <- vst(dds, blind = TRUE, fitType = "local")

# PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
    theme_minimal(base_size = 14) +
    labs(x = paste0("PC1: ", percent_var[1], "%"),
         y = paste0("PC2: ", percent_var[2], "%"),
         title = "PCA: Tumor vs Normal")

ggsave("results/qc/pca_plot.pdf", p, width = 8, height = 6)
cat("   Saved: pca_plot.pdf\n")

# Sample correlation heatmap
sample_cor <- cor(assay(vsd))
annotation_col <- data.frame(Condition = metadata$condition, row.names = rownames(metadata))

pdf("results/qc/sample_correlation.pdf", width = 10, height = 10)
pheatmap(sample_cor,
         annotation_col = annotation_col,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         main = "Sample Correlation")
dev.off()
cat("   Saved: sample_correlation.pdf\n")

# Save for next step
saveRDS(dds, "results/qc/dds_filtered.rds")
saveRDS(vsd, "results/qc/vsd.rds")

# Write summary
sink("results/qc/qc_summary.txt")
cat("=== QC Summary ===\n\n")
cat(sprintf("Total samples: %d\n", ncol(dds)))
cat(sprintf("Tumor: %d, Normal: %d\n", 
    sum(metadata$condition == "Tumor"), 
    sum(metadata$condition == "Normal")))
cat(sprintf("Genes after filtering: %d\n", nrow(dds)))
cat(sprintf("Library size: %.1f - %.1f million\n", min(lib_sizes)/1e6, max(lib_sizes)/1e6))
cat(sprintf("\nPCA variance: PC1=%.1f%%, PC2=%.1f%%\n", percent_var[1], percent_var[2]))
sink()

cat("\n=== QC Complete! ===\n")
cat("Results saved in: results/qc/\n")
