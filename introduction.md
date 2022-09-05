---
layout: page
title: Introduction & QC
nav_order: 1
parent: Home
nav_order: 1
mathjax: true
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

<img src="../assets/images/04-workflow-overview.png" width="400px" alt="workflow">

# 1. The fastq format

The first step in the RNA-Seq workflow is to take the FASTQ files received from the sequencing facility and assess the quality of the sequencing reads.

The FASTQ file format is the defacto file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read) there are four lines, each of which are described below:

| Line |	Description |
|------|-------------|
| 1 |	Always begins with ‘@’ and then information about the read |
| 2 |	The actual DNA sequence |
| 3 |	Always begins with a ‘+’ and sometimes the same info in line 1 |
| 4 |	Has a string of characters which represent the quality scores; must have same number of characters as line 2 |

## 1.1 A first peek at our FASTQ files

Several sequencing files are available in the /datasets/ folder as it contains 4 fastq files. The files are generaly quite big (they usualy contain up to 40 milion reads), so it’s a smart thing to keep them zipped as they are.

~~~
# Let’s view the directory that contains the sequencing files (.fastq.gz) and other needed files e.g. genome reference sequence.

ls /project/bims6000/data/morning/

# zcat is a simular function as cat but works on zipped files. With the use of this function we can have a look at the files without having to unzip them.

zcat /project/bims6000/data/morning/Arabidopsis_sample2.fq.gz | head -n 20

# This will show the first 20 lines of the file, containing 5 reads.
~~~
{: .language-bash}

Let’s have a close look at the first read of this sample:

~~~
@ERR1406259.27450842
CATCGCTGAAGATCTGTGAACCAGCCTTGAACCAAACTGCCTCTCCAAACTTGACTCCGTTCCTGGCCAAAAGCTCAGGGAAGACGCAGCCTAGGGCTCCG
+
?ABEEEDCBFEDGHFJFJIHFEFCC=>BEC>FJ@GHCHBHCGFJHG;:F<AI;90F=E44:8FA>@8C;;33237-?84(>*$A#$#/B.5)->0%/8D=;
~~~
{: .output}

As mentioned previously, line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. **Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger**

~~~
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40
~~~    
{: .output}

Using the quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 30. The second base (A) has a quality of 32, etc.

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

$Q=−10 \times log_{10}(P)$

where $P$ is the probability that a base call is erroneous.

These probability values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

| Phred Quality Score |	Probability of incorrect base call | Base call accuracy |
|---------------------|------------------------------------|--------------------|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1000 | 99.9% |
| 40 | 1 in 10,000 |	99.99% |

Therefore, for the first nucleotide in the read (C), there is a 1 in 1000 chance that the base was called incorrectly. Also you can see that the second half of the read contains a lot of bases that have a more then 10% probabaility that the base is called incorrectly.

Question: How many reads do these samples contain?

~~~
# Answer: To get the number of reads, get the number of lines and divide by 4.

zcat /project/bims6000/data/morning/Arabidopsis_sample2.fq.gz | wc -l

# This gives 1,000,000 lines -> 250,000 reads.
~~~
{: .language-bash}

# 2. Quality control of FASTQ files

## 2.1. Running FastQC

We will create the quality reports of the reads that were downloaded.

First, we need to make an output directory for the fastqc results to be stored. This we want to do in the ‘home’ directory that contains all the needed files.

~~~
# create a new directory
mkdir fastqc

# let's load fastqc

module load fastqc

# Running fastqc uses the following command

fastqc -o fastqc /project/bims6000/data/morning/Arabidopsis_sample1.fq.gz

# Of course we don’t want to do this for all the samples seperately so we can loop through the list of samples and run them all sequentially. Using echo, you can start off with a “dry run”:

for filename in /project/bims6000/data/morning/*.fq.gz
do
  echo fastqc -o fastqc $filename
done

# What does this do?  Why do this?

#The echo command only prints the commands to the screen, and doesn’t really run it.

# If it looks good remove the echo and go for it.

for filename in /project/bims6000/data/morning/*.fq.gz
do
  fastqc -o fastqc $filename
done

# You will see an automatically updating output message telling you the progress of the analysis.

# In total, it should take about five minutes for FastQC to run on all four of our zipped FASTQ files.

# If the command doesn’t run or you want more information on fastqc, run the following to get the help page.

fastqc -h

# But if all went right, the FastQC program will have created several new files within our /home/fastqc directory.
~~~
{: .language-bash}

## 2.2. Viewing the FastQC results

For each of the samples there are two files. a .html and a .zip

If we were working on our local computer, we’d be able to display each of these HTML files as a webpage.

Instead, we have to transfer the html files to our local computer and visualize them there.

If you have Unix/Linux running on your local computer, you can execute the following command:

~~~
scp -r [uva_compute_id]@rivanna.hpc.virginia.edu:/home/[uva_compute_id]/fastqc ~

# Enter your Netbadge password and the fastqc directory and files contained in it will transfer to your local home directory.
~~~
{: .language-bash}

If you're running Windows, you can also use `scp` within MobaXterm or PuTTY to transfer the files.

Or, you can go to the following site: https://rivanna-portal.hpc.virginia.edu/pun/sys/dashboard

Click on Files -> Home Directory and drag and drop the fastqc folder to your local computer.

Open Arabidopsis_sample1_fastqc.html by clicking on it.  It should show up on your default browser.

## 2.3 Decoding the FastQC outputs

Upon opening the file Below we have provided a brief overview of interpretations for each of these plots. It’s important to keep in mind Now that we have run FASTQC and downloaded the report, we can take a look at the metrics and assess the quality of our sequencing data!

- Per tile sequence quality: the machines that perform sequencing are divided into tiles. This plot displays patterns in base quality along these tiles. Consistently low scores are often found around the edges, but hot spots can also occur in the middle if an air bubble was introduced at some point during the run.
- Per sequence quality scores: a density plot of quality for all reads at all positions. This plot shows what quality scores are most common.
- Per base sequence content: plots the proportion of each base position over all of the reads. Typically, we expect to see each base roughly 25% of the time at each position, but this often fails at the beginning or end of the read due to quality or adapter content.
- Per sequence GC content: a density plot of average GC content in each of the reads.
- Per base N content: the percent of times that ‘N’ occurs at a position in all reads. If there is an increase at a particular position, this might indicate that something went wrong during sequencing.
- Sequence Length Distribution: the distribution of sequence lengths of all reads in the file. If the data is raw, there is often on sharp peak, however if the reads have been trimmed, there may be a distribution of shorter lengths.
- Sequence Duplication Levels: A distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.
- Overrepresented sequences: A list of sequences that occur more frequently than would be expected by chance.
- Adapter Content: a graph indicating where adapater sequences occur in the reads.

FastQC has a really well documented manual page with detailed explanations about every plot in the report.

Within our report, a summary of all of the modules is given on the left-hand side of the report. Don’t take the yellow “WARNING”s and red “FAIL”s too seriously; they should be interpreted as flags for modules to check out.

<img src="../assets/images/fastqc_summary.png" width="400px" alt="fastqc_summary">

The first module gives the basic statistics for the sample. Generally it is a good idea to keep track of the total number of reads sequenced for each sample and to make sure the read length and %GC content is as expected.

One of the most important analysis modules is the “Per base sequence quality” plot. This plot provides the distribution of quality scores at each position in the read across all reads. This plot can alert us to whether there were any problems occuring during sequencing and whether we might need to contact the sequencing facility.

The “Per sequence quality scores” plot gives you the average quality score on the x-axis and the number of sequences with that average on the y-axis. We hope the majority of our reads have a high average quality score with no large bumps at the lower quality values.

The next plot gives the “Per base sequence content”, which always gives a FAIL for RNA-seq data. This is because the first 10-12 bases result from the ‘random’ hexamer priming that occurs during RNA-seq library preparation. This priming is not as random as we might hope giving an enrichment in particular bases for these intial nucleotides.

The “Per sequence GC content” plot gives the GC distribution over all sequences. Generally is a good idea to note whether the GC content of the central peak corresponds to the expected % GC for the organism. Also, the distribution should be normal unless over-represented sequences (sharp peaks on a normal distribution) or contamination with another organism (broad peak).

This plot would indicate some type of over-represented sequence with the sharp peaks, indicating either contamination or a highly over-expressed gene.

The next module explores numbers of duplicated sequences in the library. This plot can help identify a low complexity library, which could result from too many cycles of PCR amplification or too little starting material. For RNA-seq we don’t normally do anything to address this in the analysis, but if this were a pilot experiment, we might adjust the number of PCR cycles, amount of input, or amount of sequencing for future libraries. In this analysis we seem to have a large number of duplicated sequences, but this is can be expected due to the multiple copies of mRNA being duplicates.

The “Overrepresented sequences” table is another important module as it displays the sequences (at least 20 bp) that occur in more than 0.1% of the total number of sequences. This table aids in identifying contamination, such as vector or adapter sequences. If the %GC content was off in the above module, this table can help identify the source. If not listed as a known adapter or vector, it can help to BLAST the sequence to determine the identity.

## 2.4 Working with the FastQC text output

Now that we’ve looked at our HTML report to get a feel for the data, let’s look more closely at the other output files.

~~~
# Go back to your Unix terminal and cd into your fastqc folder

cd fastqc

ls
~~~
{: .language-bash}

Our .zip files are compressed files. They each contain multiple different types of output files for a single input FASTQ file. To view the contents of a .zip file, we can use the program `unzip` to decompress these files.  `unzip` expects to get only one zip file as input. We could go through and unzip each file one at a time, but this is very time consuming and error-prone. Someday you may have 500 files to unzip!

A more efficient way is to use a for loop like we learned in the Shell Genomics lesson to iterate through all of our .zip files. Let’s see what that looks like and then we’ll discuss what we’re doing with each line of our loop.

~~~
for filename in *.zip
do
 unzip $filename
done
~~~
{: .language-bash}

In this example, the input is four filenames (one filename for each of our `.zip` files). Each time the loop iterates, it will assign a file name to the variable `filename` and run the `unzip` command. The first time through the loop, `$filename` is `Arabidopsis_sample1_fastqc.zip`. The interpreter runs the command `unzip` on `Arabidopsis_sample1_fastqc.zip`. For the second iteration, `$filename` becomes `Arabidopsis_sample2_fastqc.zip`. This time, the shell runs `unzip` on `Arabidopsis_sample2_fastqc.zip`. It then repeats this process for the four other `.zip` files in our directory.

The `unzip` program is decompressing the .zip files and creating a new directory (with subdirectories) for each of our samples, to store all of the different output that is produced by FastQC. There are a lot of files here. The one we’re going to focus on is the `summary.txt` file.

If you list the files using `ls` in our directory now you will see a mix of files and directories.  The `.html` files and the uncompressed `.zip` files are still present, but now we also have a new directory for each of our samples. We can see for sure that it’s a directory if we use the `-F` flag for `ls`.

~~~ 
ls -F

# Let’s see what files are present within one of these output directories.

ls -F Arabidopsis_sample1_fastqc/

# Use less to preview the summary.txt file for this sample.

less Arabidopsis_sample1_fastqc/summary.txt

# The summary file gives us a list of tests that FastQC ran, and tells us whether this sample passed, failed, or is borderline (WARN). Remember to quit from less you enter q.

# We can make a record of the results we obtained for all our samples by concatenating all of our summary.txt files into a single file using the cat command. We’ll call this full_summaries.txt.

cat */summary.txt > fastqc_summaries.txt
~~~
{: .language-bash}

# 3. Trimming and filtering

Before we will do the alignment we need to remove sequences of low quality and sequences that are to short (below 25 bases). Also in this case we will trim down long sequences to 100 bases, quality of the Ion-torrent reads drops the further it gets. When making use of illumina reads this is not as much of a problem and 3’-trimming would then be a waste of data.

The trimming and quality filtering will be done with trimmomatic. In the programm the following arguments can be used.

| step              |	   meaning                          |
|-------------------|-------------------------------------|
| `SE` or `PE`      |	Reads are single end or paired end. |
| `ILLUMINACLIP`    |	Perform adapter removal.             |
| `SLIDINGWINDOW`   |	Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`	        | Cut bases off the start of a read, if below a threshold quality. |
| `TRAILING`	     | Cut bases off the end of a read, if below a threshold quality. |
| `CROP`	           | Cut the read to a specified length. |
| `HEADCROP`	     | Cut the specified number of bases from the start of the read. |
| `MINLEN`	        | Drop an entire read if it is below a specified length. |
| `TOPHRED33`	     | Convert quality scores to Phred-33. |
| `TOPHRED64`	     | Convert quality scores to Phred-64. |

~~~
# Let's make a new directory named trimmed and load trimmomatic

mkdir trimmed

module load trimmomatic

# To run trimmomatic on a single sample it looks something like this

trimmomatic SE -phred33 -threads 1 /project/bims6000/data/morning/Arabidopsis_sample1.fq.gz ~/trimmed/Arabidopsis_sample1_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
~~~
{: .language-bash}

Of course, we don’t want to do this for all the reads seperately so lets create a loop through all the fastq files.

When doing the fastqc only input files needed to be specified. In this case both the input and a matching output filenames need to be given. this can be done with the help of ‘basename’

~~~
for infile in /project/bims6000/data/morning/*.fq.gz
do
 echo inputfile $infile
 outfile="$(basename $infile .fq.gz)"_qc.fq
 echo outputfile $outfile
 echo
done

# Next we can start writing the trimmomatic loop. Again starting with a dry run with echo.

for infile in /project/bims6000/data/morning/*.fq.gz
do
  outfile="$(basename $infile .fq.gz)"_qc.fq
  echo "trimmomatic SE -phred33 -threads 2 $infile ~/trimmed/$outfile ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25"
done

# If it all looks ok, rerun with out echo

for infile in /project/bims6000/data/morning/*.fq.gz
do
  outfile="$(basename $infile .fq.gz)"_qc.fq
  trimmomatic SE -phred33 -threads 2 $infile ~/trimmed/$outfile ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done
~~~
{: .language-bash}

It’s possible to scroll up to check if the percentage of surviving & dropped is within the same range in all of the samples.
