---
layout: home
title: Home
nav_order: 1
has_children: true
has_toc: false
---

# {{ site.title }}: {{ site.tagline }}
{: .mb-2 }
{: .fs-6 .fw-300 }

## Welcome!

This lesson will introduce you to the basics of gene expression analysis using RNA-Seq (short for RNA sequencing). Due to the considerable progress and constant decreasing costs of RNA-Seq, this technique has became a standard technique in biology.

It is going to be fun and empowering! We will use the shell (covered on the first day of this submodule) and R (covered on the second day of the submodule) to perform our RNA-Seq analyses and visualisations. Before you begin, be sure you are all set up for the lesson. See and complete the [Setup]({{ site.baseurl }}{% link setup.md %})  section.

## Main learning objectives

After completing this lesson, you should be able to:

* Assess the quality of RNA-seq sequencing data (“reads”) using the command-line instructions.
* Align RNA-seq reads to a reference genome using a splice-aware aligner (e.g. STAR).
* Generate a count matrix from the RNA-seq data alignment.
* Perform a QC of your experiment through Principal Component Analysis (PCA) and sample clustering.
* Execute a differential gene expression analysis using R and the DESeq2 package.
* Be able to create key plots: volcano plot, heatmap and clustering of differentially expressed genes.
* Provide a biological interpretation to differentially expressed genes through ORA/GSEA analyses and data integration.

## Schedule

|     | Setup |     | 
| --- | ---   | --- |
| 09:00 - 10:00 | [Introduction & QC]({{ site.baseurl }}{% link introduction.md %}) | What can I learn by doing this RNA-Seq lesson?<br>What are the tools that I will be using?<br>How do I perform a quality check of my RNA-seq fastq files with FastQC?<br>How can I remove RNA-seq reads of low quality? |
| 10:00 - 11:00 | [Aligning]({{ site.baseurl }}{% link aligning.md %}) | How do I align my reads to a reference genome using STAR and hisat2? | 
| 11:00 - 11:45 | [Counting]({{ site.baseurl }}{% link counting.md %})  | What is a BAM file?<br>How do I determine the number of reads that map within genes? |
| 1:00 - 2:00  | [Differential expression]({{ site.baseurl }}{% link diffexp.md %}) | How do I know that my RNA-seq experiment has worked according to my experimental design?<br>What is a Principal Component Analysis (PCA) and how can I use it?<br>What are factor levels and why is it important for different expression analysis?<br>How can I call the genes differentially regulated in response to my experimental design?<br>What is a volcano plot and how can I create one?<br>What is a heatmap and how can it be informative for my comparison of interest? |
| 2:00 - 3:00 | [Over-representation analysis]({{ site.baseurl }}{% link ora.md %}) | Given a list of differentially expressed genes, how do I search for enriched functions? |
| 3:00 - 3:45 | [Gene set enrichment]({{ site.baseurl }}{% link gsea.md %}) | What is the difference between an over-representation analysis (ORA) and a gene set enrichment analysis (GSEA)? |

## Credits

This lesson is heavily based on teaching materials from the [Harvard Chan Bioinformatics Core (HBC) in-depth NGS data analysis course](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course) and [RNA-seq lesson from the ScienceParkStudyGroup](https://github.com/ScienceParkStudyGroup/rnaseq-lesson). 

