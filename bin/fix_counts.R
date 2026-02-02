#!/usr/bin/env Rscript

cat("=== Fixing Count Data ===\n\n")

# Load the raw downloaded data
cat("1. Loading expression data...\n")
expr_data <- read.table(gzfile("data/counts/expression.gz"), 
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

cat(sprintf("   Loaded: %d genes x %d samples\n", nrow(expr_data), ncol(expr_data)))
cat(sprintf("   Value range: %.2f to %.2f\n", min(expr_data), max(expr_data)))

# Data is log2(x+1) transformed, convert back to counts
cat("\n2. Converting to counts...\n")
# First, handle any negative values (set to 0)
expr_data[expr_data < 0] <- 0

# Convert: if log2(x+1), then x = 2^value - 1
expr_counts <- round((2^expr_data) - 1)

# Make sure no negative values
expr_counts[expr_counts < 0] <- 0

cat(sprintf("   Count range: %d to %d\n", min(expr_counts), max(expr_counts)))

# Select samples
cat("\n3. Selecting samples...\n")
sample_ids <- colnames(expr_counts)
sample_type <- substr(sample_ids, 14, 15)
tumor_samples <- sample_ids[sample_type == "01"]
normal_samples <- sample_ids[sample_type == "11"]

set.seed(42)
selected_tumor <- sample(tumor_samples, min(50, length(tumor_samples)))
selected_normal <- sample(normal_samples, min(50, length(normal_samples)))
selected_samples <- c(selected_tumor, selected_normal)

counts_subset <- expr_counts[, selected_samples]

cat(sprintf("   Selected: %d tumor + %d normal\n", 
            length(selected_tumor), length(selected_normal)))

# Save
cat("\n4. Saving corrected counts...\n")
write.csv(counts_subset, "data/counts/raw_counts.csv")

cat("\n=== Done! ===\n")
cat(sprintf("Final: %d genes x %d samples\n", nrow(counts_subset), ncol(counts_subset)))
cat(sprintf("Value range: %d to %d\n", min(counts_subset), max(counts_subset)))
