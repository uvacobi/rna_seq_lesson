---
layout: page
title: Counting
nav_order: 3
parent: Home
mathjax: true
---

# 1. The SAM/BAM format
SAM files, are tab-delimited text files that contain information for each individual read and its alignment to the genome. While we do not have time to go in detail of the features of the SAM format, the paper by Heng Li et al. provides a lot more detail on the specification.

The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to allow for indexing, which enables efficient random access of the data contained within the file.

## 1.1 What’s in a SAM/BAM file

The file begins with a header, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the alignment section. Each line that follows corresponds to alignment information for a single read. Each alignment line has 11 mandatory fields for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.

Additionally tags (or attribute) can be aded to each of of the lines. These tags give some aditional information on the alignment. The number and type of tags varies between different alinment tools and the settings within these tools. Here a list of tags that are commonly used.

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
| XG:i |	Number of gap extentions |
| XT	 | Type: Unique/Repeat/N/Mate-sw |
| XA:z | Alternative hits; format: (chr,pos,CIGAR,NM;) |
| XS:i | Suboptimal alignment score |
| XF:i | Support from forward/reverse alignment |
|XE:i  | Number of supporting seeds |
