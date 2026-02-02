#!/usr/bin/env Rscript

.libPaths("~/R/libs")

library(ggplot2)
library(dplyr)
library(survival)

cat("=== Survival Analysis ===\n\n")

# Load data
cat("1. Loading data...\n")
norm_counts <- read.csv("results/differential_expression/normalized_counts.csv", 
                        row.names = 1, check.names = FALSE)
de_results <- read.csv("results/differential_expression/deseq2_results_all.csv", stringsAsFactors = FALSE)
clinical <- read.csv("data/clinical/clinical_data.csv", stringsAsFactors = FALSE)
metadata <- read.csv("config/samplesheet_counts.csv", stringsAsFactors = FALSE)

dir.create("results/survival", recursive = TRUE, showWarnings = FALSE)

# Get tumor samples only
tumor_metadata <- metadata %>% filter(condition == "Tumor")
cat(sprintf("   Tumor samples: %d\n", nrow(tumor_metadata)))

# Match to clinical data
tumor_metadata$patient_barcode <- substr(tumor_metadata$sample_id, 1, 12)
surv_data <- merge(tumor_metadata, clinical, by.x = "patient_barcode", by.y = "bcr_patient_barcode")

# Create survival variables
surv_data$time <- ifelse(!is.na(surv_data$days_to_death), 
                         surv_data$days_to_death, 
                         surv_data$days_to_last_followup)
surv_data$status <- ifelse(surv_data$vital_status == "Dead", 1, 0)
surv_data <- surv_data %>% filter(!is.na(time), time > 0)

cat(sprintf("   Samples with survival data: %d\n", nrow(surv_data)))
cat(sprintf("   Events (deaths): %d\n", sum(surv_data$status)))

# Get top DE genes
cat("\n2. Testing genes for survival association...\n")
top_genes <- de_results %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    head(50) %>%
    pull(gene_id)

top_genes <- top_genes[top_genes %in% rownames(norm_counts)]
cat(sprintf("   Testing %d genes\n", length(top_genes)))

# Test each gene
survival_results <- data.frame()

for (gene in top_genes) {
    expr <- as.numeric(norm_counts[gene, surv_data$sample_id])
    if (all(is.na(expr)) || var(expr, na.rm = TRUE) == 0) next
    
    temp_df <- data.frame(
        time = surv_data$time,
        status = surv_data$status,
        expression = expr
    )
    
    cox_fit <- tryCatch({
        coxph(Surv(time, status) ~ expression, data = temp_df)
    }, error = function(e) NULL)
    
    if (is.null(cox_fit)) next
    
    cox_sum <- summary(cox_fit)
    
    survival_results <- rbind(survival_results, data.frame(
        gene = gene,
        HR = exp(coef(cox_fit)),
        pvalue = cox_sum$coefficients[5],
        stringsAsFactors = FALSE
    ))
}

survival_results <- survival_results %>% arrange(pvalue)
survival_results$padj <- p.adjust(survival_results$pvalue, method = "BH")

write.csv(survival_results, "results/survival/gene_survival_results.csv", row.names = FALSE)
cat(sprintf("   Genes tested: %d\n", nrow(survival_results)))
cat(sprintf("   Significant (p < 0.05): %d\n", sum(survival_results$pvalue < 0.05)))

# Build risk score from top genes
cat("\n3. Building risk score...\n")
sig_genes <- survival_results %>% filter(pvalue < 0.1) %>% head(10) %>% pull(gene)

if (length(sig_genes) >= 3) {
    # Calculate risk score
    sig_expr <- norm_counts[sig_genes, surv_data$sample_id, drop = FALSE]
    risk_scores <- colMeans(log2(sig_expr + 1), na.rm = TRUE)
    
    surv_data$risk_score <- risk_scores
    surv_data$risk_group <- ifelse(risk_scores > median(risk_scores), "High", "Low")
    
    # Save risk scores
    write.csv(surv_data[, c("sample_id", "risk_score", "risk_group", "time", "status")],
              "results/survival/patient_risk_scores.csv", row.names = FALSE)
    
    # Kaplan-Meier analysis
    cat("\n4. Creating Kaplan-Meier plot...\n")
    surv_data$time_years <- surv_data$time / 365
    
    fit <- survfit(Surv(time_years, status) ~ risk_group, data = surv_data)
    
    # Log-rank test
    logrank <- survdiff(Surv(time_years, status) ~ risk_group, data = surv_data)
    pval <- 1 - pchisq(logrank$chisq, df = 1)
    
    # Create KM plot manually
    km_data <- data.frame(
        time = fit$time,
        surv = fit$surv,
        group = rep(names(fit$strata), fit$strata)
    )
    km_data$group <- gsub("risk_group=", "", km_data$group)
    
    p_km <- ggplot(km_data, aes(x = time, y = surv, color = group)) +
        geom_step(linewidth = 1.2) +
        scale_color_manual(values = c("High" = "red", "Low" = "blue")) +
        theme_minimal(base_size = 14) +
        labs(x = "Time (years)", y = "Survival Probability",
             title = "Kaplan-Meier Survival Curve",
             subtitle = sprintf("Log-rank p = %.3f", pval),
             color = "Risk Group") +
        ylim(0, 1)
    
    ggsave("results/survival/km_curve.pdf", p_km, width = 8, height = 6)
    cat("   Saved: km_curve.pdf\n")
    
    # Cox regression
    cox_risk <- coxph(Surv(time, status) ~ risk_score, data = surv_data)
    cox_sum <- summary(cox_risk)
    
    hr <- exp(coef(cox_risk))
    hr_ci <- cox_sum$conf.int[c(3, 4)]
    cox_p <- cox_sum$coefficients[5]
}

# Summary
cat("\n5. Writing summary...\n")
sink("results/survival/survival_summary.txt")
cat("=== Survival Analysis Summary ===\n\n")
cat(sprintf("Tumor samples analyzed: %d\n", nrow(surv_data)))
cat(sprintf("Events (deaths): %d (%.1f%%)\n", sum(surv_data$status), 100*mean(surv_data$status)))
cat(sprintf("Median follow-up: %.1f years\n", median(surv_data$time)/365))

cat("\n=== Gene-Level Survival Analysis ===\n")
cat(sprintf("Genes tested: %d\n", nrow(survival_results)))
cat(sprintf("Significant (p < 0.05): %d\n", sum(survival_results$pvalue < 0.05)))

cat("\nTop 10 Survival-Associated Genes:\n")
print(head(survival_results, 10), row.names = FALSE)

if (length(sig_genes) >= 3) {
    cat("\n=== Risk Score Model ===\n")
    cat(sprintf("Genes in signature: %d\n", length(sig_genes)))
    cat(sprintf("Signature genes: %s\n", paste(sig_genes, collapse = ", ")))
    cat(sprintf("\nCox Regression:\n"))
    cat(sprintf("  HR = %.2f (95%% CI: %.2f - %.2f)\n", hr, hr_ci[1], hr_ci[2]))
    cat(sprintf("  p-value = %.4f\n", cox_p))
    cat(sprintf("\nLog-rank test p-value: %.4f\n", pval))
}
sink()

cat("\n=== Survival Analysis Complete! ===\n")
cat("Results saved in: results/survival/\n")
