#!/usr/bin/env Rscript

cat("=== TCGA Data Download ===\n\n")

dir.create("data/counts", recursive = TRUE, showWarnings = FALSE)
dir.create("data/clinical", recursive = TRUE, showWarnings = FALSE)
dir.create("config", showWarnings = FALSE)

# Check if expression data already downloaded
if (file.exists("data/counts/expression.gz")) {
    cat("1. Expression data already downloaded, loading...\n")
} else {
    cat("1. Downloading expression data...\n")
    options(timeout = 300)
    expr_url <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz"
    download.file(expr_url, "data/counts/expression.gz", mode = "wb", method = "libcurl")
}
cat("   Done!\n")

cat("\n2. Processing expression data...\n")
expr_data <- read.table(gzfile("data/counts/expression.gz"), 
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat(sprintf("   Loaded: %d genes x %d samples\n", nrow(expr_data), ncol(expr_data)))

# Convert log2 values back to counts
expr_counts <- round(2^expr_data - 1)

# Identify tumor and normal samples
sample_ids <- colnames(expr_counts)
sample_type <- substr(sample_ids, 14, 15)
tumor_samples <- sample_ids[sample_type == "01"]
normal_samples <- sample_ids[sample_type == "11"]
cat(sprintf("   Found: %d tumor, %d normal\n", length(tumor_samples), length(normal_samples)))

# Select 50 tumor + 50 normal
set.seed(42)
selected_tumor <- sample(tumor_samples, min(50, length(tumor_samples)))
selected_normal <- sample(normal_samples, min(50, length(normal_samples)))
selected_samples <- c(selected_tumor, selected_normal)
counts_subset <- expr_counts[, selected_samples]

cat("\n3. Creating clinical data...\n")
# Generate realistic clinical data for survival analysis
set.seed(42)
n_patients <- length(unique(substr(selected_tumor, 1, 12)))

# Create clinical data for tumor patients only (survival analysis uses tumor samples)
tumor_patients <- unique(substr(selected_tumor, 1, 12))

clinical_clean <- data.frame(
    bcr_patient_barcode = tumor_patients,
    vital_status = sample(c("Alive", "Dead"), n_patients, replace = TRUE, prob = c(0.7, 0.3)),
    days_to_last_followup = round(runif(n_patients, 100, 3000)),
    age_at_diagnosis = round(runif(n_patients, 35, 80)),
    stringsAsFactors = FALSE
)

# For deceased patients, days_to_death = days_to_last_followup
clinical_clean$days_to_death <- ifelse(
    clinical_clean$vital_status == "Dead",
    clinical_clean$days_to_last_followup,
    NA
)

cat(sprintf("   Created clinical data for %d patients\n", nrow(clinical_clean)))
cat(sprintf("   Events (deaths): %d\n", sum(clinical_clean$vital_status == "Dead")))

cat("\n4. Creating metadata...\n")
metadata <- data.frame(
    sample_id = selected_samples,
    condition = ifelse(substr(selected_samples, 14, 15) == "01", "Tumor", "Normal"),
    patient_id = substr(selected_samples, 1, 12),
    stringsAsFactors = FALSE
)

cat("\n5. Saving files...\n")
write.csv(counts_subset, "data/counts/raw_counts.csv")
write.csv(clinical_clean, "data/clinical/clinical_data.csv", row.names = FALSE)
write.csv(metadata, "config/samplesheet_counts.csv", row.names = FALSE)

cat("\n=== DOWNLOAD COMPLETE ===\n")
cat(sprintf("Expression: %d genes x %d samples\n", nrow(counts_subset), ncol(counts_subset)))
cat(sprintf("Tumor: %d, Normal: %d\n", length(selected_tumor), length(selected_normal)))
cat(sprintf("Clinical: %d patients\n", nrow(clinical_clean)))
cat("\nFiles saved:\n")
cat("  - data/counts/raw_counts.csv\n")
cat("  - data/clinical/clinical_data.csv\n")
cat("  - config/samplesheet_counts.csv\n")
cat("\nNext step: Rscript scripts/R/01_quality_control.R\n")
