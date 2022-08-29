---
layout: page
title: Over-representation analysis
nav_order: 5
parent: Home
mathjax: true
---

# 1. Introduction

## 1.1 From a list of genes to biological insights

You've finally managed to extract a list of differentially expressed genes from your comparison. Great job! But...now what??

Why did you do the experiment in the first place? Probably because you had an hypothesis or you were looking for new avenues to conduct experiments. 

A functional enrichment analysis will determine whether some functions are enriched in your set of differentially expressed genes.

In this tutorial, we are looking for Arabidopsis leaf genes that are induced or repressed upon inoculation by _Pseudomonas syringae_ DC3000 after 7 days.

One important goal is to gain a higher view and not only deal with individual genes but understand which pathways are involved in the response.

Once we obtain a list of genes, we have multiple analysis to perform to go beyond a simple list of genes:
- Annotating our list of genes with cross-databases identifiers and descriptions (Entrezid, Uniprot, KEGG, etc.)
- Performing Over-Representation Analysis (ORA)
- Performing Gene Set Enrichment Analysis (GSEA)
- Interpreting the results

These ORA and GSEA analysis require the use of external resources to assign functions to genes. Two resources are of particular importance and will be examined today.

## 1.2 Over Representation Analysis (ORA)

Over Representation Analysis is searching for biological functions or pathways that are enriched in a list obtained through experimental studies compared to the complete list of functions/pathways.

Usually, ORA makes use of so-called gene ontologies (abbreviated GO) where each gene receives one or multiple layers of information on their function, cellular localization, etc.

The ORA analysis rely on this mathematical equation to compute a p-value for a given gene set classified under a certain GO.

Source: [Wikipedia](https://en.wikipedia.org/wiki/Hypergeometric_distribution)

![formula](../assets/images/hypergeometric.png)

In this formula:
- **N** is the total number of genes in the background distribution. Also called the "universe" of our transcriptome.
- **M** is the number of genes within that distribution that are annotated (either directly or indirectly) to the gene set of interest.
- **n** is the size of the list of genes of interest (the size of your "drawing").
- **k** is the number of genes "drawn" within that list which are annotated to the gene set.

The background distribution by default is by default all genes that have annotation. You can change it to your specific background if you have a good reason for that (only genes with a detectable expression in your expression for instance). Also, p-values should be adjusted for multiple comparison.

See this [chapter](https://yulab-smu.github.io/clusterProfiler-book/chapter2.html) for more info.

## 1.3 The Gene Ontology (GO) resource

The __Gene Ontology (GO)__ produces a bird's-eye view of biological systems by building a tree of terms related to biological functions. This is particularly helpful when dealing with results from genome-wide experiments (e.g. transcriptomics) since classifying genes into groups of related functions can assist in the interpretation of results. Rather than focusing on each gene, one by one, the researcher gets access to metabolic pathways, functions related to development, etc.

> ## Note
> The GO resource is divided into 3 main subdomains:
1. __Biological Process (BP)__: a series of molecular events with a defined beginning and end relevant for the function of an organism, a cell, etc.
2. __Cellular Component (CC)__: the part of a cell.
3. __Molecular Function (MF)__: the enzymatic activites of a gene product.
[Source Wikipedia](https://en.wikipedia.org/wiki/Gene_ontology).
{: .callout}

Let's take an example. The At3g53260 gene codes for a phenylalanine ammonia-lyase (PAL) that catalyses the following reaction and is one of the first step of cell wall synthesis, flavonoid synthesis, etc. ; L-phenylalanine ⇌ trans-cinnamic acid + NH<sub>3</sub>. We can go to https://www.arabidopsis.org and search for this gene on the top right:

This gene has several GO terms associated:

<img src="../assets/images/At3g53260.png" height="600px" alt="AT3G53260 page on arabidopsis.org">

Here are an example term associated with each GO subdomain:
1. __BP:__ the "L-phenylalanine catabolic process" term with the GO:0006559 unique identifier.
2. __CC__: the "cytoplasm" term with the GO:0005737 unique identifier.
3. __MF__: the "ammonia-lyase activity" term with the GO:0016841 unique identifier.

## 1.4 The Kyoto Encyclopedia of Genes and Genomes (KEGG) database

KEGG stands for the "Kyoto Encyclopedia of Genes and Genomes". From the [KEGG website home page](https://www.genome.jp/kegg/):
> KEGG is a database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism and the ecosystem, from molecular-level information, especially large-scale molecular datasets generated by genome sequencing and other high-throughput experimental technologies.

Instead of using the Gene Ontology gene classification, one might be interested to use KEGG classification to view the
transcriptomic response of an organism. KEGG is not restricted to metabolic functions but has a great deal of metabolic maps that can help you.

> ## Important note
> While using a model organism such as _Arabidopsis thaliana_ makes ORA and GSEA analyses easier, it is noteworthy that the GO and KEGG resources are not restricted to model organisms but rather include a huge number of (plant) species.
{: .callout}

<br>

# 2. Gene Ontology ORA analysis using clusterProfiler (R code)

Gene information can be automatically queried directly from R to access the [Ensembl](https://www.ensembl.org) databases. Ensembl gathers a tremendous amount of genomic information which can be accessed through a web browser or programmatically.

> The Ensembl (https://www.ensembl.org) is a system for generating and distributing genome annotation such as genes, variation, regulation and comparative genomics across the vertebrate subphylum and key model organisms. The Ensembl annotation pipeline is capable of integrating experimental and reference data from multiple providers into a single integrated resource. Here, we present 94 newly annotated and re-annotated genomes, bringing the total number of genomes offered by Ensembl to 227.

<img src="../assets/images/ensembl.png" alt="Portal of Ensembl for Arabidopsis" width="1000px">

We are going to use two fantastic resources: the [Ensembl](https://www.ensembl.org) database and the [biomartr package](https://docs.ropensci.org/biomartr). Together, they will automate a lot of tedious and tiring steps when you want to retrieve gene annotations, sequences, etc.

This step is meant to retrieve the correspondence between organism-specific gene identifiers (e.g. At1g01020) and NCBI Entrez Gene ID (e.g. 839321) which are used by `clusterProfiler`.

We are going to load the required library first.
~~~
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
library("org.At.tair.db")
library("biomaRt")  # only use to remove cache bug
~~~
{: .language-r}

<br>

> ## Important note: troubleshooting
> If biomart refuses to query Ensembl again, run this command:
> ~~~
> biomaRt::biomartCacheClear() # to solve a known bug https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/335
> ~~~
> {: .language-r}
> This will clean the cache memory and allow to perform the Ensembl query again.
{: .callout}

## 2.1 Load the table of differential genes

If not done yet, load the table of differential genes.
~~~
diff_genes <- read_delim(file = "/project/bims6000/data/afternoon/differential_genes.csv", delim = ",")
~~~
{: .language-r}

All what we know about the differential genes are their locus identifier. Not much.... We are missing functional information which we will add.

## 2.2 Annotating your DE genes with Ensembl and biomartr

What purpose serves `biomartr`? From the documentation:
> The first step, however, of any genome based study is to retrieve genomes and their annotation from databases. To automate the retrieval process of this information on a meta-genomic scale, the biomartr package provides interface functions for genomic sequence retrieval and functional annotation retrieval. The major aim of biomartr is to facilitate computational reproducibility and large-scale handling of genomic data for (meta-)genomic analyses. In addition, biomartr aims to address the genome version crisis. With biomartr users can now control and be informed about the genome versions they retrieve automatically. Many large scale genomics studies lack this information and thus, reproducibility and data interpretation become nearly impossible when documentation of genome version information gets neglected.

What is available for _Arabidopsis thaliana_ in Ensembl?
~~~
# library("biomartr") (if not loaded already)
biomartr::organismBM(organism = "Arabidopsis thaliana")
~~~
{: .language-r}

~~~
organism_name description                                                  mart         dataset       version
  <chr>         <chr>                                                        <chr>        <chr>         <chr>
1 athaliana     Arabidopsis thaliana genes (TAIR10)                          plants_mart  athaliana_eg… TAIR10
2 athaliana     Arabidopsis thaliana Short Variants (SNPs and indels exclud… plants_vari… athaliana_eg… TAIR10
~~~
{: .output}

This indicates that we can get a dataset called `athaliana_eg_gene` of the genome annotation version `TAIR10` from the `plant_mart` mart.

Let's see how many different information fields we can retrieve from the `arabidopsis_eg_gene` dataset.
~~~
arabido_attributes <- 
  biomartr::organismAttributes("Arabidopsis thaliana") %>%
  filter(dataset == "athaliana_eg_gene")

arabido_attributes
~~~
{: .language-r}


~~~
# A tibble: 2,574 x 4
   name                  description              dataset           mart
   <chr>                 <chr>                    <chr>             <chr>
 1 ensembl_gene_id       Gene stable ID           athaliana_eg_gene plants_mart
 2 ensembl_transcript_id Transcript stable ID     athaliana_eg_gene plants_mart
 3 ensembl_peptide_id    Protein stable ID        athaliana_eg_gene plants_mart
 4 ensembl_exon_id       Exon stable ID           athaliana_eg_gene plants_mart
 5 description           Gene description         athaliana_eg_gene plants_mart
 6 chromosome_name       Chromosome/scaffold name athaliana_eg_gene plants_mart
 7 start_position        Gene start (bp)          athaliana_eg_gene plants_mart
 8 end_position          Gene end (bp)            athaliana_eg_gene plants_mart
 9 strand                Strand                   athaliana_eg_gene plants_mart
10 band                  Karyotype band           athaliana_eg_gene plants_mart
# … with 3,200 more rows
~~~
{: .output}

There is definitely a lot information in there! We should be able to get what we want = the correspondence between the Arabidopsis gene identifier and the NCBI Entrez Gene identifier.

~~~
attributes_to_retrieve = c("tair_symbol", "entrezgene_id")

result_BM <- biomartr::biomart( genes      = diff_genes$genes,                  # genes were retrieved using biomartr::getGenome()
                                mart       = "plants_mart",                     # marts were selected with biomartr::getMarts()
                                dataset    = "athaliana_eg_gene",               # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
                                filters =   "ensembl_gene_id" )# query key
head(result_BM)
~~~
{: .language-r}

We now have our original gene identifiers (column `ensembl_gene_id`) with the retrieved TAIR symbols (`tair_symbol`) and NCBI Entrez Gene Id (`entrezgene_id`).
~~~
  ensembl_gene_id  tair_symbol entrezgene_id
1       AT1G01030         NGA3        839321
2       AT1G01070                     839550
3       AT1G01090 PDH-E1 ALPHA        839429
4       AT1G01140        CIPK9        839349
5       AT1G01220         FKGP        839420
6       AT1G01225                     839358
~~~
{: .output}

> ## For other species
> If your species is not "Arabidopsis thaliana", simply change your R code here:
> ~~~
> # library("biomartr") (if not loaded already)
> biomartr::organismBM(organism = "[my favorite species]")
> ~~~
> {: .language-r}
> This will retrieve the information available for your species on Ensembl.
{: .callout}

## 2.3 ORA with clusterProfiler

To perform the ORA within R, we will use the [clusterProfiler Bioconductor package](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) that has an [extensive documentation available here](https://yulab-smu.github.io/clusterProfiler-book/index.html).

First, we need to annotate both genes that make up our "universe" and the genes that were identified as differentially expressed.
~~~
# building the universe!
all_arabidopsis_genes <- read.csv("/project/bims6000/data/afternoon/raw_counts.csv", header = TRUE, stringsAsFactors = FALSE)[,1] # directly selects the gene column

# we want the correspondence of TAIR/Ensembl symbols with NCBI Entrez gene ids
attributes_to_retrieve = c("tair_symbol", "uniprotswissprot","entrezgene_id")

# Query the Ensembl API
all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                     mart       = "plants_mart",
                                                     dataset    = "athaliana_eg_gene",
                                                     attributes = attributes_to_retrieve,
                                                     filters =  "ensembl_gene_id" )

# for compatibility with enrichGO universe
# genes in the universe need to be characters and not integers (Entrez gene id)
all_arabidopsis_genes_annotated$entrezgene_id = as.character(all_arabidopsis_genes_annotated$entrezgene_id)

~~~
{: .language-r}

We now have a correspondence for all our genes found in Arabidopsis.

~~~
# retrieving NCBI Entrez gene id for our genes called differential
diff_arabidopsis_genes_annotated <- biomartr::biomart(genes = diff_genes$genes,
                                                     mart       = "plants_mart",
                                                     dataset    = "athaliana_eg_gene",
                                                     attributes = attributes_to_retrieve,
                                                     filters =  "ensembl_gene_id" )
~~~
{: .language-r}

This gave us the second part which is the classification of genes "drawn" from the whole gene universe. The "drawing" is coming from the set of genes identified as differential (see [episode 06](../06-differential-analysis/index.html)).

~~~
# performing the ORA for Gene Ontology Biological Process class
ora_analysis_bp <- enrichGO(gene = diff_arabidopsis_genes_annotated$entrezgene_id,
                            universe = all_arabidopsis_genes_annotated$entrezgene_id,
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE,
                            pool = FALSE)
~~~
{: .language-r}

Since we have 3 classes for GO terms i.e. Molecular Function (MF), Cellular Component (CC) and Biological Processes (BP), we have to run this 3 times for each GO class.

> ## Exercise
> How many GO categories do you find overrepresented (padj < 0.05) for the Cellular Component and Molecular Function classes?
{: .challenge}

The Gene Ontology classification is very redundant meaning that parental terms overlap a lot with their related child terms. The `clusterProfiler` package comes with a dedicated function called `simplify` to solve this issue.

~~~
# clusterProfiler::simplify to disambiguate which simplify() function you want to use
ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp)
~~~
{: .language-r}

The `ora_analysis_bp_simplified` is a rich and complex R object. It contains various layers of information (R object from the S4 class). Layers can be accessed through the "@" notation.

You can extract a nice table of results for your next breakthrough publication like this.
~~~
write_delim(x = as.data.frame(ora_analysis_bp_simplified@result),
            file = "go_results.tsv",
            delim = "\t")

# have a look at a few columns and rows if you'd like.
ora_analysis_bp_simplified@result[1:5,1:8]
~~~
{: .language-r}

~~~
           ONTOLOGY         ID                         Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
GO:0009753       BP GO:0009753           response to jasmonic acid  257/3829 474/20450 1.099409e-68 1.772247e-65 1.160744e-65
GO:0009611       BP GO:0009611                response to wounding  187/3829 335/20450 1.252868e-52 1.009812e-49 6.613824e-50
GO:0006612       BP GO:0006612       protein targeting to membrane  199/3829 374/20450 1.683735e-51 4.690438e-49 3.072032e-49
GO:0010243       BP GO:0010243 response to organonitrogen compound  222/3829 443/20450 1.707866e-51 4.690438e-49 3.072032e-49
GO:0072657       BP GO:0072657    protein localization to membrane  200/3829 377/20450 1.745821e-51 4.690438e-49 3.072032e-49
~~~
{: .output}

## 2.4 Plots from the Gene Ontology ORA analysis

Nice to have all this textual information but an image is worth a thousand words so let's create some visual representations.

A dotplot can be created very easily.
~~~
tiff("ora_dotplot.tiff")
dotplot(ora_analysis_bp_simplified)
dev.off()
~~~
{: .language-r}

<img src="../assets/images/ora_dotplot.tiff" width="800px">

You can also create an enrichment map that connects GO terms with edges between overlapping gene sets.
This makes it easier to identify functional modules.

~~~
ora_analysis_bp <- pairwise_termsim(ora_analysis_bp, method = "JC")

tiff("enrichment_map.tiff", width=10, height=10, res=100, units="in")
emapplot(ora_analysis_bp, color = "qvalue")
dev.off()
~~~
{: .language-r}

<img src="../assets/images/enrichment_map.tiff" width="600px">

On this plot, we can see that one major module related to cell death, the immune response etc. is to be seen along with two minor modules
related to metabolism (upper left) and one related to jasmonic acid and wounding (bottom).


> ## Important note
> Remember to perform the analysis for all GO categories:
> - Biological Process (`ont = "BP"`),
> - Cellular Component (`ont = "CC"`),
> - Molecular Function (`ont = "MF"`).
{: .callout}

<br>

# 3. KEGG Over Representation Analysis using clusterProfiler 

## 3.1 Retrieving species-specific KEGG information

First things first, load the required libraries if not done yet.

~~~
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
suppressPackageStartupMessages(library("org.At.tair.db"))
library("biomaRt")  # only use to remove cache bug
~~~
{: .language-r}

To see if your organism is referenced in the KEGG database, you can search this page: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)
In our case, _Arabidopsis thaliana_ is referenced as "ath" in the KEGG database.

You can also do this programmatically using R and the `clusterProfiler` package.
~~~
search_kegg_organism('ath', by='kegg_code')
search_kegg_organism('Arabidopsis thaliana', by='scientific_name')
~~~
{: .language-r}

## 3.2 KEGG ORA analysis

Performing the ORA analysis is then quite similar to what we've done with the [`GO analysis with clusterProfiler`](/07-functional-enrichment/index.html#34-ora-with-clusterprofiler).

~~~
ora_analysis_kegg <- enrichKEGG(gene = diff_arabidopsis_genes_annotated$entrezgene_id,
                                universe = all_arabidopsis_genes_annotated$entrezgene_id,
                                organism = "ath",
                                keyType = "ncbi-geneid",
                                minGSSize = 10,
                                maxGSSize = 500,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05,
                                use_internal_data = FALSE) # force to query latest KEGG db

~~~
{: .language-r}

We can then create a dotplot to visualise the KEGG categories significantly enriched.
~~~
# create a simple dotplot graph
tiff("ora_kegg.tiff")
dotplot(ora_analysis_kegg,
    color = "qvalue",
    showCategory = 10,
    size = "Count")
dev.off()
~~~
{: .language-r}

<img src="../assets/images/ora_kegg.tiff" width="800px">

## 3.3 KEGG Modules ORA

[The KEGG MODULE datase](https://www.genome.jp/kegg/module.html) is a series of "manually defined functional units of gene sets". In particular, pathway modules are functional units of gene sets in metabolic pathways that can give a metabolic-centric view of differentially expressed genes.

The complete list of available modules is [available here](https://www.genome.jp/kegg-bin/show_brite?ko00002.keg).

~~~
ora_analysis_kegg_modules <- enrichMKEGG(gene = diff_arabidopsis_genes_annotated$entrezgene_id,
                                         universe = all_arabidopsis_genes_annotated$entrezgene_id,
                                         organism = "ath",
                                         keyType = "ncbi-geneid",
                                         minGSSize = 10,           # minimal size of genes annotated by Ontology term for testing.
                                         maxGSSize = 500,          # maximal size of genes annotated for testing
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 0.05)

~~~
{: .language-r}

Similarly, we can plot this ORA result as a dotplot.
~~~
# create a simple dotplot graph
tiff("module_dotplot.tiff")
dotplot(ora_analysis_kegg_modules,
    color = "qvalue",
    showCategory = 10,
    size = "Count")
dev.off()
~~~
{: .language-r}

<img src="../assets/images/module_dotplot.tiff" width="800px">


> ## Discussion
> Compare the two KEGG plots. Can you identify differences? Which metabolic functions have been grouped together?
{: .discussion}

# 8. Going further

## 8.1 Useful links
- [BiomartR](https://docs.ropensci.org/biomartr/)
- [Arabidopsis.org (TAIR) list of data mining tools](https://www.arabidopsis.org/portals/expression/microarray/microarrayExpressionV2.jsp)
- [ResearchGate related question](https://www.researchgate.net/post/How_can_I_analyze_a_set_of_DEGs_differentially_expressed_genes_to_obtain_information_from_them)

## 8.2. References
* [The Cluster Profiler companion book, a great place to start](https://yulab-smu.github.io/clusterProfiler-book/chapter2.html)
* Zhou et al. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 10, 1523 (2019). [link](https://doi.org/10.1038/s41467-019-09234-6)
* Yates et al. (2020) Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, [Link](https://doi.org/10.1093/nar/gkz966)
* Tian et al. (2017) agriGO v2.0: a GO analysis toolkit for the agricultural community. _Nucleic Acids Research_, Volume 45, Issue W1, Pages W122–W129.[Link](https://doi.org/10.1093/nar/gkx382)
* Drost et al. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ 33(8): 1216-1217. [doi:10.1093/bioinformatics/btw821](https://academic.oup.com/bioinformatics/article/33/8/1216/2931816).
