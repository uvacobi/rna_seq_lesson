---
layout: page
title: Setup
description: Setup for completion of RNA-seq analysis
nav_order: 2
---

## Setup

### Computing on Rivanna

University of Virginia’s High-Performance Computing (HPC) system includes two large clusters named Rivanna and Afton. As a centralized resource the HPC has hundreds of pre-installed software packages available for computational research across many disciplines. All of us have logged on to Rivanna before this class. For today's lesson, each one of us will request 2 CPUs on Rivanna for the next 8 hours. We can do that by issuing the following shell command

```bash
ijob -c 2 --mem-per-cpu=9000 -A ${account} -p standard --time=08:00:00
```

We will replace the `${account}` in the command with the name of the group we will specifically create for the course.


### Working directory

We will create a folder named `my_rna_seq_analysis` when we log onto Rivanna. Remember that you can use the following command to create the folder

```bash
mkdir my_rna_seq_analysis
```

Now that we have the folder, we will `cd` into it, and do all our analyses there

```bash
cd my_rna_seq_analysis
```


### Software tools

We have installed all the tools that will be required for this analysis on the Rivanna cluster. Rivanna uses [Environment Modules](http://modules.sourceforge.net) to let users use tools without having to install them. 

### R and R packages

We have also installed the R packages we will need for this lesson in a directory "/standard/bims6000/R". 

### Data files

We have also downloaded the various files that you are going to need for this lesson on to Rivanna.

**What you need for the fastq QC, alignment, and counting**

The following files can be found at /standard/bims6000/data/morning/

- **Arabidopsis_sample1/2/3/4.fq.gz**: A `FASTQ` file containing a sample sequenced mRNA-seq reads in the FASTQ format.
- **AtChromosome1.fa**: the chromosome 1 sequence of the Arabidopsis thaliana genome in `FASTA` format.  
- **ath_annotation.gff3**: the genome annotation of Arabidopsis thaliana for chromosome 1 in the `GFF3` format. This indicates the positions of genes, their exons and 5' or 3' UTR on the chromosome and is used to generate the gene counts.   
- **adapters.fasta**: the Illumina adapter sequences used for read trimming using Trimmomatic. 
{: .prereq}

<br>

**What you need for the differential expression and enrichment analysis**

The following files can be found at /standard/bims6000/data/afternoon

- **Counts**: A `raw_counts.csv` dataframe of the sample raw counts. It is a tab separated file therefore data are in tabulated separated columns.
- **Samples to experimental conditions**: the `samples_to_conditions.csv` dataframe indicates the correspondence between samples and experimental conditions (e.g. control, treated).  
**Differentially expressed genes**: `differential_genes.csv` dataframe contains the result of the DESeq2 analysis.  

## Original study

This RNA-seq lesson will make use of a dataset from a study on the model plant Arabidopsis thaliana inoculated with commensal leaf bacteria (Methylobacterium extorquens or Sphingomonas melonis) and infected or not with a leaf bacterial pathogen called Pseudomonas syringae. Leaf samples were collected from Arabidopsis plantlets from plants inoculated or not with commensal bacteria and infected or not with the leaf pathogen either after two days (2 dpi, dpi: days post-inoculation) or seven days (6 dpi).

All details from the study are available in [Vogel et al. in 2016 and was published in New Phytologist](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14036).
