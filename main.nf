#!/usr/bin/env nextflow

/*
 * TCGA RNA-seq Biomarker Discovery Pipeline
 * 
 * Author: Asuwin Rajaganesh
 * Description: End-to-end RNA-seq analysis for cancer biomarker discovery
 */

nextflow.enable.dsl = 2

// Pipeline parameters
params.counts = null
params.metadata = null
params.clinical = null
params.outdir = "results"

// Print pipeline info
log.info """
    ╔══════════════════════════════════════════════════════╗
    ║   TCGA RNA-seq Biomarker Discovery Pipeline          ║
    ╠══════════════════════════════════════════════════════╣
    ║   Counts    : ${params.counts}
    ║   Metadata  : ${params.metadata}
    ║   Clinical  : ${params.clinical}
    ║   Output    : ${params.outdir}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()

// Validate inputs
if (!params.counts || !params.metadata) {
    error "Please provide --counts and --metadata"
}

/*
 * Quality Control
 */
process QC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "dds_filtered.rds", emit: dds
    path "vsd.rds", emit: vsd
    path "*.pdf"
    path "qc_summary.txt"

    script:
    """
    #!/usr/bin/env Rscript
    
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    
    counts <- read.csv("${counts}", row.names = 1, check.names = FALSE)
    metadata <- read.csv("${metadata}", stringsAsFactors = FALSE)
    rownames(metadata) <- metadata\$sample_id
    metadata\$condition <- factor(metadata\$condition, levels = c("Normal", "Tumor"))
    counts <- counts[, rownames(metadata)]
    
    dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = metadata, design = ~ condition)
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    vsd <- vst(dds, blind = TRUE, fitType = "local")
    
    # PCA Plot
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var <- round(100 * attr(pca_data, "percentVar"))
    
    p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
        geom_point(size = 4) +
        scale_color_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
        labs(x = paste0("PC1: ", pct_var[1], "%"), y = paste0("PC2: ", pct_var[2], "%")) +
        theme_minimal()
    ggsave("pca_plot.pdf", p, width = 8, height = 6)
    
    # Correlation heatmap
    pdf("sample_correlation.pdf", width = 10, height = 10)
    pheatmap(cor(assay(vsd)), annotation_col = data.frame(Condition = metadata\$condition, row.names = rownames(metadata)))
    dev.off()
    
    saveRDS(dds, "dds_filtered.rds")
    saveRDS(vsd, "vsd.rds")
    
    sink("qc_summary.txt")
    cat("Samples:", ncol(dds), "\\nGenes:", nrow(dds), "\\n")
    sink()
    """
}

/*
 * Differential Expression
 */
process DIFFERENTIAL_EXPRESSION {
    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
    path dds
    path vsd

    output:
    path "deseq2_results_all.csv", emit: results
    path "normalized_counts.csv", emit: norm_counts
    path "dds_analyzed.rds"
    path "*.pdf"
    path "de_summary.txt"

    script:
    """
    #!/usr/bin/env Rscript
    
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    
    dds <- readRDS("${dds}")
    vsd <- readRDS("${vsd}")
    
    dds <- DESeq(dds, fitType = "local")
    res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
    res_df <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)
    res_df <- res_df[order(res_df\$padj), ]
    
    write.csv(res_df, "deseq2_results_all.csv", row.names = FALSE)
    write.csv(counts(dds, normalized = TRUE), "normalized_counts.csv")
    saveRDS(dds, "dds_analyzed.rds")
    
    # Volcano plot
    res_df\$sig <- ifelse(!is.na(res_df\$padj) & res_df\$padj < 0.05 & abs(res_df\$log2FoldChange) > 1,
                         ifelse(res_df\$log2FoldChange > 0, "Up", "Down"), "NS")
    
    p <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
        theme_minimal() +
        labs(title = "Volcano Plot: Tumor vs Normal")
    ggsave("volcano_plot.pdf", p, width = 10, height = 8)
    
    # Heatmap
    top_genes <- head(res_df\$gene_id[!is.na(res_df\$padj)], 50)
    mat <- assay(vsd)[top_genes, ]
    mat_scaled <- t(scale(t(mat)))
    mat_scaled[mat_scaled > 2] <- 2
    mat_scaled[mat_scaled < -2] <- -2
    
    pdf("heatmap_top50.pdf", width = 10, height = 12)
    pheatmap(mat_scaled, show_colnames = FALSE)
    dev.off()
    
    n_up <- sum(res_df\$padj < 0.05 & res_df\$log2FoldChange > 1, na.rm = TRUE)
    n_down <- sum(res_df\$padj < 0.05 & res_df\$log2FoldChange < -1, na.rm = TRUE)
    
    sink("de_summary.txt")
    cat("Upregulated:", n_up, "\\nDownregulated:", n_down, "\\n")
    sink()
    """
}

/*
 * Survival Analysis
 */
process SURVIVAL_ANALYSIS {
    publishDir "${params.outdir}/survival", mode: 'copy'

    input:
    path norm_counts
    path de_results
    path clinical

    output:
    path "*.csv"
    path "*.pdf"
    path "survival_summary.txt"

    script:
    """
    #!/usr/bin/env Rscript
    
    library(survival)
    library(ggplot2)
    library(dplyr)
    
    norm_counts <- read.csv("${norm_counts}", row.names = 1, check.names = FALSE)
    de_results <- read.csv("${de_results}")
    clinical <- read.csv("${clinical}")
    
    # Prepare survival data
    clinical\$time <- ifelse(!is.na(clinical\$days_to_death), clinical\$days_to_death, clinical\$days_to_last_followup)
    clinical\$status <- ifelse(clinical\$vital_status == "Dead", 1, 0)
    clinical <- clinical[!is.na(clinical\$time) & clinical\$time > 0, ]
    
    # Get top genes
    top_genes <- de_results %>% filter(!is.na(padj), padj < 0.05) %>% head(20) %>% pull(gene_id)
    top_genes <- top_genes[top_genes %in% rownames(norm_counts)]
    
    # Build risk score
    if (length(top_genes) >= 3) {
        tumor_samples <- colnames(norm_counts)[grep("-01\$", colnames(norm_counts))]
        sig_expr <- norm_counts[top_genes[1:min(5, length(top_genes))], tumor_samples, drop = FALSE]
        risk_scores <- colMeans(log2(sig_expr + 1))
        
        surv_df <- data.frame(
            sample_id = names(risk_scores),
            patient = substr(names(risk_scores), 1, 12),
            risk_score = as.numeric(risk_scores)
        )
        surv_df <- merge(surv_df, clinical, by.x = "patient", by.y = "bcr_patient_barcode")
        surv_df\$risk_group <- ifelse(surv_df\$risk_score > median(surv_df\$risk_score), "High", "Low")
        surv_df\$time_years <- surv_df\$time / 365
        
        write.csv(surv_df, "patient_risk_scores.csv", row.names = FALSE)
        
        # KM plot
        fit <- survfit(Surv(time_years, status) ~ risk_group, data = surv_df)
        km_data <- data.frame(time = fit\$time, surv = fit\$surv, 
                              group = rep(names(fit\$strata), fit\$strata))
        km_data\$group <- gsub("risk_group=", "", km_data\$group)
        
        p <- ggplot(km_data, aes(time, surv, color = group)) +
            geom_step(linewidth = 1.2) +
            scale_color_manual(values = c("High" = "red", "Low" = "blue")) +
            labs(x = "Time (years)", y = "Survival Probability", title = "Kaplan-Meier Curve") +
            theme_minimal()
        ggsave("km_curve.pdf", p, width = 8, height = 6)
    }
    
    sink("survival_summary.txt")
    cat("Samples:", nrow(surv_df), "\\nEvents:", sum(surv_df\$status), "\\n")
    sink()
    """
}

/*
 * Main workflow
 */
workflow {
    counts_ch = Channel.fromPath(params.counts)
    metadata_ch = Channel.fromPath(params.metadata)
    clinical_ch = Channel.fromPath(params.clinical)

    QC(counts_ch, metadata_ch)
    DIFFERENTIAL_EXPRESSION(QC.out.dds, QC.out.vsd)
    SURVIVAL_ANALYSIS(DIFFERENTIAL_EXPRESSION.out.norm_counts, 
                      DIFFERENTIAL_EXPRESSION.out.results, 
                      clinical_ch)
}
