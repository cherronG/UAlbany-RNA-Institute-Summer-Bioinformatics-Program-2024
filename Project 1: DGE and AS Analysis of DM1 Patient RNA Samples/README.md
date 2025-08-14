# Project 1: Differential Gene Expression and Alternative Splicing Analysis of Myotonic Dystrophy Patient RNA Samples

## Project Overview
This project focuses on conducting bioinformatic differential gene expression and alternative splicing analysis of myotonic dystrophy (DM1) patient RNA-Seq data. The primary objective was to identify significantly up-regulated genes in patients with myotonic dystrophy and practice fundamental bioinformatics programming skills.

**Presentation Date:** July 9, 2024

## Background
Myotonic Dystrophy (MD) is a multi-systemic autosomal dominant disorder affecting tissues throughout the body.
Key characteristics include:
- **DM1:** more severe form caused by CUG RNA repeat expansions
- **DM2:** less severe form caused by CCUG RNA repeat expansions
- **Pathogenesis:** varying lengths of simple repeat expansions (micro-satellites) lead to genomic shifts
- **Molecular Mechanism:** expansions result in muscle-blind (MBNL) proteins being sequestered by toxic RNA, leading to decreased MBNL protein concentration, splicing changes, and muscle defect symptoms

## Research Objectives
- Identify differentially expressed genes between DM1 patients and control samples
- Analyze alternative splicing events associated with myotonic dystrophy
- Examine MBNL1-related splicing patterns
- Practice essential bioinformatics analysis techniques

## Methods and Analytical Pipeline
### Data Processing (HPCC)
- **Input Data:** FASTQ files containing RNA sequences and Phred quality scores
- **Sequence Alignment:** STAR alignment tool to map FASTQ reads to human reference genome (resulting in SAM and BAM files)

### DGE Analysis (RStudio/DESeq)
- **PCA Plot:** sample variance visualization
- **MA Plot:** overall differential expression trends
- **Volcano Plot:** Log2FoldChange vs adjusted p-value
- **Heatmap:** top 20 differentially expressed genes
- **Normalized Counts Plot:** individual gene expression patterns

### AS Analysis (rMATS/HPCC)
- **Filtered PCA Plot:** genes with FDR <= 0.05 and delta PSI >= 0.1
- **Violin Plot:** MBNL1 splicing events analysis
- **Modified Sashimi Plot:** exon inclusion level visualization
  
## Key Findings
### DGE Analysis
- **Overall Pattern:** majority of genes significantly up-regulated in DM1 samples
- **Top Significance Gene:** MYH3 showed the highest statistical significance
- **Sample Clustering:** clear separation between control and DM1 samples in PCA analysis
- **Expression Patterns:** top 20 differentially expressed genes showed consistent up-regulation in DM1 samples

### AS Analysis
- **MBNL1 Splicing Events:** identified significant splicing differences between DM1 and control samples
- **Event 59859:** showed the most variation between sample groups with significantly different PSI values
- **Exon Inclusion:** higher PSI values in DM1 samples suggest splicing events promoted by the disease condition

## Repository Contents
This folder contains:
- complete presentation slides with results and interpretations
- R scripts for differential gene expression analysis and code generation for plots
- alternative splicing analysis results and visualization scripts

---
_This project served as fundamental training in RNA-seq data analysis and bioinformatics programming as part of the 2024 RNA Institute Summer Bioinformatics Program._
