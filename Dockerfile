# TCGA RNA-seq Biomarker Pipeline
# Dockerfile for reproducible analysis

FROM rocker/r-ver:4.3.1

LABEL maintainer="Asuwin Rajaganesh <rajaganesh.a@northeastern.edu>"
LABEL description="RNA-seq analysis pipeline for cancer biomarker discovery"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'dplyr', 'pheatmap', \
    'RColorBrewer', 'survival', 'ggrepel'), repos='https://cloud.r-project.org')"

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('DESeq2', 'clusterProfiler', 'org.Hs.eg.db', \
    'EnhancedVolcano', 'fgsea'), ask=FALSE, update=FALSE)"

# Set working directory
WORKDIR /pipeline

# Copy scripts
COPY scripts/ /pipeline/scripts/
COPY bin/ /pipeline/bin/

# Default command
CMD ["R"]
