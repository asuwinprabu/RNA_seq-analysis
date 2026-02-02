# TCGA Breast Cancer RNA-seq Biomarker Discovery Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-available-blue.svg)](https://hub.docker.com/)
[![R](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible RNA-seq analysis pipeline for identifying prognostic biomarkers in breast cancer using TCGA data.

## 🎯 Project Overview

**Objective:** Identify gene expression biomarkers that distinguish tumor from normal tissue and predict patient survival outcomes in breast invasive carcinoma.

**Dataset:** TCGA-BRCA (The Cancer Genome Atlas - Breast Invasive Carcinoma)
- 50 Primary Tumor samples
- 50 Adjacent Normal samples
- ~20,000 genes analyzed

## 📊 Key Results

### Differential Expression Analysis

| Metric | Value |
|--------|-------|
| Total genes tested | 2,757 |
| Significant genes (FDR < 0.05, \|log2FC\| > 1) | **1,786** |
| Upregulated in tumor | 688 |
| Downregulated in tumor | 1,098 |

**Top Upregulated Genes:**
| Gene | log2FC | Adjusted p-value |
|------|--------|------------------|
| COL11A1 | 7.35 | 4.94e-48 |
| COL10A1 | 7.23 | 7.44e-56 |
| MMP11 | 5.92 | 1.36e-40 |
| MMP13 | 5.16 | 1.32e-21 |
| COMP | 4.79 | 3.44e-32 |

**Top Downregulated Genes:**
| Gene | log2FC | Adjusted p-value |
|------|--------|------------------|
| MYOC | -8.48 | 5.60e-55 |
| FIGF | -6.30 | 2.92e-78 |
| SCARA5 | -6.24 | 7.44e-56 |
| CD300LG | -6.10 | 6.31e-55 |
| LPL | -6.04 | 2.39e-53 |

### Survival Analysis

| Metric | Value |
|--------|-------|
| Tumor samples analyzed | 50 |
| Events (deaths) | 22 (44%) |
| Median follow-up | 3.8 years |

**5-Gene Prognostic Signature:** MEOX2, FXYD1, GPIHBP1, SVEP1, ABCA9

| Metric | Value |
|--------|-------|
| Hazard Ratio | **1.74** (95% CI: 1.13 - 2.67) |
| Cox p-value | **0.012** |
| Log-rank p-value | 0.056 |

## 🛠️ Pipeline Overview
```
Input Data → QC → Differential Expression → Pathway Analysis → Survival Analysis
    │         │            │                      │                   │
  TCGA     DESeq2       DESeq2              Enrichr/GO           Cox/KM
  counts    PCA       Volcano plot          Enrichment          Risk score
```

## 🚀 Quick Start

### Option 1: Run with Nextflow + Docker (Recommended)
```bash
# Clone repository
git clone https://github.com/asuwinprabu/tcga-rnaseq-biomarker-pipeline.git
cd tcga-rnaseq-biomarker-pipeline

# Download data
Rscript bin/download_data_simple.R

# Run pipeline
nextflow run main.nf \
    --counts data/counts/raw_counts.csv \
    --metadata config/samplesheet_counts.csv \
    --clinical data/clinical/clinical_data.csv \
    --outdir results \
    -profile docker
```

### Option 2: Run Scripts Manually
```bash
# 1. Download data
Rscript bin/download_data_simple.R

# 2. Quality control
Rscript scripts/R/01_quality_control.R

# 3. Differential expression
Rscript scripts/R/02_differential_expression.R

# 4. Pathway analysis
Rscript scripts/R/03_pathway_analysis_simple.R

# 5. Survival analysis
Rscript scripts/R/04_survival_analysis.R
```

### Option 3: Run on HPC (SLURM)
```bash
nextflow run main.nf \
    --counts data/counts/raw_counts.csv \
    --metadata config/samplesheet_counts.csv \
    --clinical data/clinical/clinical_data.csv \
    -profile slurm
```

## 📁 Project Structure
```
tcga-rnaseq-pipeline/
├── main.nf                 # Nextflow pipeline
├── nextflow.config         # Pipeline configuration
├── Dockerfile              # Container definition
├── bin/                    # Data download scripts
│   └── download_data_simple.R
├── scripts/R/              # Analysis scripts
│   ├── 01_quality_control.R
│   ├── 02_differential_expression.R
│   ├── 03_pathway_analysis_simple.R
│   └── 04_survival_analysis.R
├── results/                # Output files
│   ├── qc/
│   ├── differential_expression/
│   ├── pathway_analysis/
│   └── survival/
└── docs/                   # Documentation
```

## 🔧 Technologies

| Category | Tools |
|----------|-------|
| **Workflow** | Nextflow |
| **Container** | Docker, Singularity |
| **QC** | DESeq2, PCA |
| **DE Analysis** | DESeq2 |
| **Pathway** | Enrichr, GO, KEGG |
| **Survival** | survival (R), Cox regression |
| **Visualization** | ggplot2, pheatmap |
| **HPC** | SLURM (Northeastern Discovery) |

## 📈 Methods

### Quality Control
- Library size normalization assessment
- PCA for sample clustering
- Sample correlation heatmap
- Outlier detection

### Differential Expression
- DESeq2 median-of-ratios normalization
- Wald test with Benjamini-Hochberg correction
- Thresholds: FDR < 0.05, |log2FC| > 1

### Survival Analysis
- Cox proportional hazards regression
- Kaplan-Meier survival curves
- Log-rank test for group comparison
- Multi-gene risk score development

## 📋 Requirements

### Software
- Nextflow >= 21.04.0
- Docker or Singularity
- R >= 4.3.0

### R Packages
- DESeq2
- ggplot2, dplyr
- survival
- pheatmap

## 📚 References

1. Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.

2. The Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. *Nature*, 490:61-70.

## 👤 Author

**Asuwin Rajaganesh**  
MS Bioinformatics | Northeastern University

- 📧 rajaganesh.a@northeastern.edu
- 💼 [LinkedIn](https://www.linkedin.com/in/asuwin-prabu-raja-ganesh-78b2381b1)
- 🐙 [GitHub](https://github.com/asuwinprabu)

---

*Developed as part of bioinformatics portfolio for co-op opportunities in pharma/biotech.*
