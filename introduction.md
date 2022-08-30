---
layout: page
title: Introduction & QC
nav_order: 1
parent: Home
nav_order: 1
---

## What you will learn in this lesson?

1. How to assess the quality of my RNA-Seq experiment at the sample-level
   - Using FastQC to perform quality checks on each sample fastq sequencing file.
   - What are some count normalization methods? Why are RPKM, FPKM and TPM not adequate methods.
   - Creating a PCA plot to visualise the grouping of samples in relation to the experimental factors being investigated.
2. How to perform a differential expression analysis on RNA-Seq results using R
   - Raw counts are used for differential expression not scaled counts.
   - Creating a DESeq2 object requires 3 items: the raw gene counts, the sample to condition correspondence and a formula for testing.
   - How does the DESeq method works? What are the outputs obtained using DESeq2?
   - What are the typical outputs that one can obtain from a differential gene expression analysis?
   - A table of genes being differentially regulated between two conditions.
   - A volcano plot shows the relationship between log2 fold change and the adjusted p-value for each gene.
3. How to go beyond a list of differential genes and interpret its biological meaning
   - By performing an over-representation analysis (ORA), one can find pathways or categories where differential genes are significantly more abundant.
   - By performing a gene set enrichment analysis (GSEA), one first ranks differentially expressed genes before comparing enrichment scores for whole pathways.

## Dataset used 

We will make use of a published experimental dataset from a study made on the small model plant Arabidopsis thaliana by Vogel et al. (2016). This study compares the response of 4 weeks old plantlets to different bacteria that live on the leaves of different plant species:

- A known foliar pathogen called Pseudomonas syringae strain DC3000.
- A commensal (“neutral”) bacteria called Methylobacterium extorquens strain PA1.
- A commensal (“neutral”) bacteria called Sphingomonas melonis strain Fr1.
