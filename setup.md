---
layout: page
title: Setup
description: Setup for completion of RNA-seq analysis
nav_order: 2
---

## Setup

### Computing on Rivanna

Rivanna is the University of Virginia’s High-Performance Computing (HPC) system. As a centralized resource it has hundreds of pre-installed software packages available for computational research across many disciplines. All of us have logged on to Rivanna before this class. For today's lesson, each one of us will request 2 CPUs on Rivanna for the next 8 hours. We can do that by issuing the following shell command

```bash
ijob -c 2 --mem-per-cpu=6000 -A ${account} -p standard --time=08:00:00
```

We will replace the `${account}` in the command with the name of the group we will specifically create for the course.


### Software tools

We have installed all the tools that will be required for this analysis on the Rivanna cluster. Rivanna uses [Environment Modules](http://modules.sourceforge.net) to let users use tools without having to install them. For example, let us say we want use `STAR` to align our sequences against a reference genome. If we type the following on the command line

```bash
STAR
```
 
the shell will return with a message

```text
-bash: STAR: command not found
```

The shell is informing us that it does not know if `STAR` is installed, or how to get to it. In order to check if `STAR` is available on Rivanna, we will request information about `STAR` by using the following command on the shell

```bash
module spider STAR
```

This informs us that Rivanna has versions 2.5.3a, 2.7.2b, and 2.7.9a for `STAR` available. If we want to use version 2.7.9a, we can next ask for the right way to load that particular version. We do that by using the following command on the shell

```bash
module spider star/2.7.9a
```

This command informs us of the modules that need to be loaded before `star/2.7.9a` can be loaded. Based on the response to the above command we can load `STAR` version 2.7.9a using the following on the command line

```bash
module load star/2.7.9a
```

Now if we type

```bash
STAR
```

we get a help message from the aligner. Success!


### R and R packages

We have also installed the R packages we will need for this lesson in a directory "/project/bims6000/R". In order to tell R where to look for those packages, we will do the following on the command line (shell) 

```bash
export R_LIBS=/project/bims6000/R
```

Now, we will use the modules to inform Rivanna that we want to use `R` version 4.1.1 by using the following command

```bash
module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1
```

Now, when we type 

```bash
R
```

on the command line, we should see a prompt similar to this

![prompt](../assets/images/prompt.png)


### Data files

We have also downloaded the various files that you are going to need for this lesson on to Rivanna.

**What you need for the fastq QC, alignment, and counting**

The following files can be found at /project/bims6000/data/morning/

- **Arabidopsis_sample1/2/3/4.fq.gz**: A `FASTQ` file containing a sample sequenced mRNA-seq reads in the FASTQ format.
- **AtChromosome1.fa**: the chromosome 1 sequence of the Arabidopsis thaliana genome in `FASTA` format.  
- **ath_annotation.gff3**: the genome annotation of Arabidopsis thaliana for chromosome 1 in the `GFF3` format. This indicates the positions of genes, their exons and 5' or 3' UTR on the chromosome and is used to generate the gene counts.   
- **adapters.fasta**: the Illumina adapter sequences used for read trimming using Trimmomatic. 
{: .prereq}

<br>

**What you need for the differential expression and enrichment analysis**

The following files can be found at /project/bims6000/data/afternoon

- **Counts**: A `raw_counts.csv` dataframe of the sample raw counts. It is a tab separated file therefore data are in tabulated separated columns.
- **Samples to experimental conditions**: the `samples_to_conditions.csv` dataframe indicates the correspondence between samples and experimental conditions (e.g. control, treated).  
**Differentially expressed genes**: `differential_genes.csv` dataframe contains the result of the DESeq2 analysis.  

## Original study

This RNA-seq lesson will make use of a dataset from a study on the model plant Arabidopsis thaliana inoculated with commensal leaf bacteria (Methylobacterium extorquens or Sphingomonas melonis) and infected or not with a leaf bacterial pathogen called Pseudomonas syringae. Leaf samples were collected from Arabidopsis plantlets from plants inoculated or not with commensal bacteria and infected or not with the leaf pathogen either after two days (2 dpi, dpi: days post-inoculation) or seven days (6 dpi).

All details from the study are available in [Vogel et al. in 2016 and was published in New Phytologist](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14036).




