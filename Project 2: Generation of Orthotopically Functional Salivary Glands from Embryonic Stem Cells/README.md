# Project 2: Generation of Orthotopically Functional Salivary Glands from Embryonic Stem Cells

## Project Overview
This independent research project investigates the development of mouse embryonic stem cells into functional salivary glands. The study focuses on comparing induced salivary glands (iSGs) transplanted with and without mesenchyme cells to determine if mesenchyme addition enhances maturation rates.

**Presentation Date:** August 2, 2024

**Laboratory:** Larson Lab - Developmental Biology

**Faculty Mentor:** Dr. Melinda Larson

## Background
### Scientific Context
- **Stem Cell Differentiation:** converting embryonic stem cells into specialized salivary gland tissue
- **Transcription Factor Induction:** iSGs were generated using Sox2 and Foxc1 transcription factors
### Sample Groups
- **iSG:** induced salivary glands without mesenchyme cells (control group)
- **iSG+mes**: induced salivary glands with added mesenchyme cells (experimental group)
### Hypothesis
Addition of mesenchyme cells would promote faster maturation of transplanted induced salivary glands, potentially through enhanced cellular signaling and structural support.

## Research Question
### Does the addition of mesenchyme cells to the transplanted induced Salivary Glands (iSGs) promote faster maturation?

## Methods and Analytical Pipeline
### Data Processing (HPCC)
- Gathered ReadsPerGene files from experimental RNA-seq data
### DGE Analysis (RStudio)
- Followed standard differential gene expression (DGE) analysis pipeline
- Generated comprehensive plots using DESeq2 in RStudio
  -  PCA plot for sample clustering
  -  MA and Volcano plots for differential expression patterns
  -  Heatmaps and normalized counts plot for detailed expression patterns
### GO Enrichment Analysis (Metascape)
- Performed gene ontology (GO) enrichment analysis using Metascape for functional characterization
  
## Key Findings
### DGE Analysis 
- Limited differential gene expression between iSG and iSG+mes groups
- Only one significantly up-regulated gene identified: **Mir668**
- Small variance between groups with no clear clustering by condition
### GO Enrichment Analysis
Immune System Enhancement:
- **H2-Q7** up-regulation in iSG+mes samples indicates enhanced T-cell tumor-targeted immune response
- Promoted immune health in iSG+mes samples suggests mesenchyme cells may contribute to immunological maturation
### Research Question Assessment
- Cannot confidently conclude that mesenchyme cells promote faster maturation
- Only minimal transcriptional differences observed

## Repository Contents
This folder contains:
- complete presentation slides with results and interpretations
- R scripts for differential gene expression analysis and code generation for plots

---
_This independent research project represents advanced application of bioinformatics techniques to address fundamental questions in developmental biology and regenerative medicine as part of the 2024 RNA Institute Summer Bioinformatics Program._
