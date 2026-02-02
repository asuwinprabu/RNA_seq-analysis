#!/usr/bin/env Rscript

.libPaths("~/R/libs")

library(TCGAbiolinks)
library(SummarizedExperiment)

cat("=== Downloading TCGA-BRCA Data ===\n\n")

# Query
cat("1. Querying GDC...\n")
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

# Download
cat("2. Downloading (30-60 min)...\n")
GDCdownload(query, method = "api", files.per.chunk = 20)

# Prepare
cat("3. Preparing data...\n")
data <- GDCprepare(query)

# Extract counts
counts <- assay(data, "unstranded")
sample_info <- as.data.frame(colData(data))

# Select 50 tumor + 50 normal
set.seed(42)
tumor_idx <- which(sample_info$sample_type == "Primary Tumor")
normal_idx <- which(sample_info$sample_type == "Solid Tissue Normal")

selected_tumor <- sample(tumor_idx, min(50, length(tumor_idx)))
selected_normal <- sample(normal_idx, min(50, length(normal_idx)))
selected <- c(selected_tumor, selected_normal)

counts_subset <- counts[, selected]
sample_info_subset <- sample_info[selected, ]

# Create metadata
metadata <- data.frame(
    sample_id = colnames(counts_subset),
    condition = ifelse(sample_info_subset$sample_type == "Primary Tumor", "Tumor", "Normal"),
    patient_id = sample_info_subset$patient,
    stringsAsFactors = FALSE
)

# Save
dir.create("data/counts", recursive = TRUE, showWarnings = FALSE)
dir.create("config", showWarnings = FALSE)

write.csv(counts_subset, "data/counts/raw_counts.csv")
write.csv(metadata, "config/samplesheet_counts.csv", row.names = FALSE)

cat("\n=== Done! ===\n")
cat(sprintf("Saved: %d tumor + %d normal samples\n", length(selected_tumor), length(selected_normal)))
