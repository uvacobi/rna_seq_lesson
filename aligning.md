---
layout: page
title: Aligning
nav_order: 2
parent: Home
mathjax: true
---

# 1. Alignment to a reference genome

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. In this tutorial we will be using STAR but also a tool like hisat2 does the job.

## STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

Seed searching Clustering, stitching, and scoring Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):

The different parts of the read that are mapped separately are called ‘seeds’. So the first MMP that is mapped to the genome is called seed1.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be seed2.

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

If STAR does not find an exact matching sequence for each part of the read due to mismatches or indels, the previous MMPs will be extended.

If extension does not give a good alignment, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

## Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).

## The alignment process consists of two steps:

Indexing the reference genome
Aligning the reads to the reference genome

## 1.1 Index the reference genome

Our first step is to index the reference genome for use by STAR. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment (index files are not exchangeable between tools).

Take note that depending on the genome size these index files produced by STAR can be pretty big. Make sure there’s enough disk space available.

~~~
# Make new directory genomeIndex and index Arabadopsis chromosome using STAR

mkdir genomeIndex

module load star

star --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles /project/bims6000/data/morning/AtChromosome1.fa --runThreadN 2

# The indexing should have produced 8 star index files. Use the following command to see if they’re really there.

ls -l genomeIndex/
~~~
{: .language-bash}

## 1.2 Align reads to reference genome

In some tools like hisat2 creating the sequence alignment files (bam-files) is done in two steps. first the aligning it self. After that the alignment file will be filtered for instance to only contain the reads that actualy map to the genome. This is done with sam flags in samtools view (with the ‘-F 4’ all the unmapped reads will be removed). STAR on the other hand has a build in filter and also a sort function. So the output is ready to use for downstream tools.

First of course we will need to create a directory to output the alignment files

~~~
mkdir mapped
~~~
{: .language-bash}

Running STAR to align ( or map ) the reads and optionaly filter and sort them.

In contrast to most tools, STAR does not have a help function. running STAR -h or STAR –help will result in an error. For information on what arguments to use you can use have a look at the STAR manual..

Here are some examples of common used arguments.

| argument |	meaning |
|----------|----------|
| `--runThreads` |	number of threads |
| `--genomeDir` | /path/to/genomeDir |
| `--readFilesIn` |	/path/to/read1 [/path/to/read2] |
| `--readFilesCommand zcat` |	when making use of gzipped fastq files |
| `--outFileNamePrefix` |	/path/to/output file name |
| `--outSAMtype` |	BAM/SAM or None [optional: SortedByCoordinate] |
| `--outReadsUnmapped` |	[default: None] Fastx ; output in separate fasta/fastq file |
| `--outFilterMultimapNmax` |	[default: 10] max number of alignments accepted |
| `--outFilterMismatchNmax` |	[default: 10] max number of mismatches accepted |
| `--outFilterMismatchNoverLmax` |	[default: 0.3] max fraction of mismatches mapped length |
| `--outFilterMismatchNoverReadLmax` |	[default: 1.0] max fraction of mismatches read length |
| `--alignEndsType` |	EndToEnd force end-to-end alignment, don’t soft-clip |

~~~
# For now we will be using STAR with the following arguments

STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/Arabidopsis_sample1_qc.fq --outFileNamePrefix mapped/Arabidopsis_sample1_qc --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All

# Next we want to make a loop to do all the files

# It’s good again to first start with a ‘dry’ run with the use of echo

for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)"_
   echo "STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All"
 done
 
 # If the commands look good, rerun but this time without the echo.
 
 for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)"_
   STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All
 done
 
 # The final.out file contains all the characteristics of the alignment, resulting in a table containing all the alignment values.
 
 less mapped/Arabidopsis_sample1_qc.final.out
 ~~~
 {: .language-bash}
 
 ## 1.3 Align reads to reference genome using hisat2
 
Alternatively it is possible to map the reads using hisat2. This tools works simular to star and gives a simular output. The commands are just a bit different. 
 
~~~
# Let's create a new genomeIndex and mapped directory and load hisat2

mkdir mapped_hisat2

mkdir genomeIndex_hisat2

module load hisat2

Just like with star the genome/chromosome needs to be indexed.

hisat2-build -p 2 /project/bims6000/data/morning/AtChromosome1.fa genomeIndex_hisat2/AtChromosome1
~~~
{: .language-bash}

