---
layout: page
title: Counting
nav_order: 3
parent: Home
mathjax: true
---

# 1. The SAM/BAM format
[SAM](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam) files, are tab-delimited text files that contain information for each individual read and its alignment to the genome. While we do not have time to go in detail of the features of the SAM format, the paper by [Heng Li et al.](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) provides a lot more detail on the specification.

The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

## 1.1 What’s in a SAM/BAM file

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.

<img src="../assets/images/sam_bam_1.png" width="600px" alt="alignment_star1">

Additionally tags (or attribute) can be added to each of of the lines. These tags give some additional information on the alignment. The number and type of tags varies between different alignment tools and the settings within these tools. Here a list of tags that are commonly used.

| Tag:Type |	Meaning |
|----------|----------|
| NM:i |	Edit distance |
| MD:i |	Mismatching positions/bases |
| AS:i |	Alignment score |
| BC:z |	Barcode sequence |
| X0:i |	Number of best hits |
| X1:i |	Number of suboptimal hits found by BWA |
| XN:i |	Number of ambiguous bases in the reference |
| XM:i |	Number of mismatches in the alignment |
| XO:i | Number of gap opens |
| XG:i |	Number of gap extensions |
| XT	 | Type: Unique/Repeat/N/Mate-sw |
| XA:z | Alternative hits; format: (chr,pos,CIGAR,NM;) |
| XS:i | Suboptimal alignment score |
| XF:i | Support from forward/reverse alignment |
|XE:i  | Number of supporting seeds |

To start of we’ll have a look at how to use samtools to have a peak at the contents of the bam files.

As these file are binary you can not simply use:

~~~
head Arabidopsis_sample1.bam 
~~~
{: .language-bash}

This will give an unreadable result. SAMtools can help us to make the content readable.

# 2. SAMtools

SAMtools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

Like many Unix commands, `samtools` commands follow a stream model, where data runs through each command as if carried on a conveyor belt. This allows combining multiple commands into a data processing pipeline. Although the final output can be very complex, only a limited number of simple commands are needed to produce it. If not specified, the standard streams (stdin, stdout, and stderr) are assumed. Data sent to stdout are printed to the screen by default but are easily redirected to another file using the normal Unix redirectors (`>` and `>>`), or to another command via a pipe (`|`).

## 2.1 SAMtools commands

SAMtools provides the following commands (in **bold**), each invoked as “samtools some_command”.

- **view** 
The view command filters SAM or BAM formatted data. Using options and arguments it understands what data to select (possibly all of it) and passes only that data through. Input is usually a sam or bam file specified as an argument, but could be sam or bam data piped from any other command. Possible uses include extracting a subset of data into a new file, converting between BAM and SAM formats, and just looking at the raw file contents. The order of extracted reads is preserved.

- **sort**
Sort alignments by leftmost coordinates, by read name when -n is used, by tag contents with -t, or a minimiser-based collation order with -M. The sorted output is dumped to a new file by default, although it can be directed to stdout (using the -o option). As sorting is memory intensive and BAM files can be large, this command supports a sectioning mode (with the -m options) to use at most a given amount of memory and generate multiple output file. These files can then be merged to produce a complete sorted BAM file.

- **index**
The index command creates a new index file that allows fast look-up of data in a (sorted) SAM or BAM. Like an index on a database, the generated \*.sam.sai or \*.bam.bai file allows programs that can read it to more efficiently work with the data in the associated files.

- **tview**
The tview command starts an interactive ASCII-based viewer that can be used to visualize how reads are aligned to specified small regions of the reference genome. Compared to a graphics based viewer like IGV,[3] it has few features. Within the view, it is possible to jumping to different positions along reference elements (using ‘g’) and display help information (‘?’).

- **mpileup**
The mpileup command produces a pileup format (or BCF) file giving, for each genomic coordinate, the overlapping read bases and indels at that position in the input BAM files(s). This can be used for SNP calling for example.

- **flagstat** 
Counts the number of alignments for each FLAG type.
Looking at the content of the file using samtools view:

~~~
# Let's cd into the mapped_hisat2 directory first

cd mapped_hisat2

samtools view Arabidopsis_sample1.bam | head
~~~
{: .language-bash}

SAMtools will make the data readable, this data is then piped through head to show the first 10 lines of the file.

## 2.2 Counting and sorting

SAMtools view can be used to filter the alignment based on characters like mapping quality, chromosome, orientation etc. When the -c option is added the filtered selection is counted.

~~~
# Count the total number of records.

samtools view -c Arabidopsis_sample1.bam 

# Count with flagstat for additional information.

samtools flagstat Arabidopsis_sample1.bam

# Count the records using the FLAG argument. Count the alignments that don’t align.

samtools view -f 4 -c Arabidopsis_sample1.bam

# The argument -f includes reads that fit samflag 4, read unmapped.

# Count the reads that do align.

samtools view -F 4 -c Arabidopsis_sample1.bam

# Here -F is used to exclude reads that fit samflag 4, read unmapped. Everything else is included.
~~~
{: .language-bash}

Question: Sometimes you will see that this number of alignments is higher then the number of sequences in your fastq file. How can this be?

Answer: When a read multimaps (aligned to multiple positions in the genome), each of these positions is included as a separate alignment.

~~~
# Write a .bam file with the above filters to a new file.

samtools view -Sb -F 4 -o Arabidopsis_sample1_mapped.bam Arabidopsis_sample1.bam

# In this command -Sb is needed to keep the file binairy(compressed), and -o specifies the output filename (a bam file again).

# Count the reads that align to the forward strand.

samtools view -F 20 -c Arabidopsis_sample1.bam

# Use -F 20 to exclude “read reverse strand” and “read unmapped”.

# Count the reads that align to the reverse strand.

samtools view -f 16 -c Arabidopsis_sample1.bam

# With -f 16 you select for “read reverse strand”.

# With SAMtools it is also posible to select for alignments with a minimal mapping quality.

# Alignments with a maximal score (60 for hisat2 output files and 255 for STAR output files) are truly unique.

samtools view -q 60 -c Arabidopsis_sample1.bam

# This number can never be bigger then the number of reads in the fastq file, as all reads in the output give a single alignment.
~~~
{: .language-bash}

Question: BAM files usually contain a tag or attribute that gives the number of mismatches between the read and the reference genome. With SAMtools it is unfortunately not possible to filter on these values. Could you think of an other way to select for alignments that align without any mismatches?
Hint: make use of `grep "XM:i:0"` among others.

~~~
# Answer:

samtools view Arabidopsis_sample1.bam | grep "XM:i:0" | wc -l
~~~
{: .language-bash}

# 3. Creating the counts file

For downstream applications e.g. differential expression analysis, the number of reads that maps within a gene has to be determined for each sample.

The `featureCounts` program from the [Subread](http://subread.sourceforge.net/) package can do this. The complete user guide is available here and `featureCounts` is in section 6.2.

`featureCounts` can…count (!) the number of reads that map within a feature. The Arabidopsis genome annotation in the GFF3 format contain three different features to choose from.

Depending on the downstream applications the choice is `gene`, `transcript` or `exon`. In this study we are just looking for differentially expressed genes so our feature of interest specified by the `-t` will be `gene`.

~~~
# Let's cd back into our working directory which is ~/my_rna_seq_analysis

cd ~/my_rna_seq_analysis

# Let's use featureCounts to count the reads in genes.  We downloaded it, and it's in /standard/bims6000/bin

/standard/bims6000/bin/featureCounts -O -t gene -g ID -a /standard/bims6000/data/morning/ath_annotation.gff3 -o counts.txt mapped_hisat2/*.bam
~~~
{: .language-bash}

Here is an explanation of the different arguments used:

- `-a <string>`: Name of an annotation file. GTF/GFF format by default.
- `-o <string>`: Name of the output file including read counts.
- `-O`: Assign reads to all their overlapping meta-features.
- `-t <string>`: Specify feature type in the GTF/GFF annotation to summarise the counts.
- `-g <string>`: Specify attribute type in GTF/GFF annotation. This GTF/GFF determines the name of the features.

The output file `counts.txt` produced by featureCounts is a tab-delimited file that can be opened in a spreadsheet program like Excel.

~~~
# Let's view the counts.txt file using less

less counts.txt

# Let's view the counts.txt.summary file using less

less counts.txt.summary
~~~
{: .language-bash}
