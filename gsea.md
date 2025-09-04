---
layout: page
title: Gene set enrichment analysis
nav_order: 6
parent: Home
---

# 1. Introduction

## ORA vs GSEA

The outcome of the over-representation test in the last section depends on the significance threshold used to declare the genes as differentially expressed. As you might have noticed, that analysis uses the number of genes in the statistical tests, but does not include the measured changes in determination of the enrichment. Functional categories in which many genes exhibit small changes may go undetected in such analyses. We also know that genes are not independent, which violates a key assumption of the Fisher's exact tests used in such over-representation analyses.

The Gene Set Enrichment Analysis (GSEA) is another way to investigate functional enrichment of genes and pathways using the Gene Ontology classification. It calculates a gene-level statistic (e.g, fold change, Wald statistic) and then aggregates the gene-level statistics for all genes in a pathway into a single pathway-level statistic. All genes can be used in GSEA; GSEA aggregates the per gene statistics across genes within a gene set, therefore making it possible to detect situations where all genes in a predefined set change in a small but coordinated way. This is important since it is likely that many relevant phenotypic differences are manifested by small but consistent changes in a set of genes. [This section](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#gsea-algorithm) provides a clear explanation of GSEA. 

Similar to GO enrichment, we can use clusterProfiler to run GSEA. First, let us load all the required libraries

~~~
library("biomartr")
library("clusterProfiler")
library("tidyverse")
library("enrichplot")
library("org.At.tair.db")
library("biomaRt")  # only use to remove cache bug
~~~
{: language-r}

Then, let's read in the DESeq2 result for all the genes. Remember we saved this file in the section on differential expression.

~~~
deseq_results <- read.csv("all_genes.csv", row.names=1)
~~~
{: language-r}

Similar to the analyses presented in ORA, we want the entrez-id for the genes in order to run the GSEA analysis. So let's fetch those using biomaRt

~~~
# Query the Ensembl API to get entrez id's for all the genes
attributes_to_retrieve = c("entrezgene_id")

all_arabidopsis_genes_annotated <- biomartr::biomart(genes = rownames(deseq_results),
                                                     mart       = "plants_mart",
                                                     dataset    = "athaliana_eg_gene",
                                                     attributes = attributes_to_retrieve,
                                                     filters =  "ensembl_gene_id" )

# for compatibility, genes need to be characters and not integers (Entrez gene id)
all_arabidopsis_genes_annotated$entrezgene_id = as.character(all_arabidopsis_genes_annotated$entrezgene_id)
~~~
{: language-r}

GSEA tests whether a pre-defined set of genes (ex: those belonging to a specific GO term or KEGG pathway) show up more frequently than expected by chance at the top or bottom of a sorted gene list from our experiment. Let's sort the genes in our experiment based on log-fold changes.

~~~
gene_list <- deseq_results %>% as_tibble(rownames="gene") %>% left_join(all_arabidopsis_genes_annotated, by=c("gene"="ensembl_gene_id")) %>% dplyr::select(entrezgene_id, log2FoldChange) %>% na.omit() %>% group_by(entrezgene_id) %>% summarise(log2FoldChange = mean(log2FoldChange)) %>% arrange(desc(log2FoldChange)) %>% deframe()
~~~
{: language-r}

Now, let's run GSEA using clusterProfiler

~~~
gsea_bp <- gseGO(gene_list,
                 ont = "BP",
                 keyType = "ENTREZID",
                 OrgDb = org.At.tair.db)
~~~

The Gene Ontology classification is very redundant meaning that parental terms overlap a lot with their related child terms. The `clusterProfiler` package comes with a dedicated function called `simplify` to solve this issue.

~~~
gsea_bp_simplified <- clusterProfiler::simplify(gsea_bp)

~~~
{: language-r}

We can plot the results of the analysis using a dot-plot similar to how we did for ORA

~~~
tiff("gsea_bp_dotplot.tiff", width=8, height=10, units="in", res=72)
dotplot(gsea_bp, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()
~~~

The dotplot below shows some of the same results as the ORA analysis (response to oxidative stress, response to wounding, ...), but we also see some new enrichments 

![dotplot](../assets/images/gsea_bp_dotplot.tiff)

We can also generate a plot of the Running Enrichment Score for a gene set as the analysis walks down the ranked gene list, including the location of the maximum enrichment score. Let's plot this for "response to oxidative stress", came out as significant in ORA and GSEA analysis.

~~~

which(gsea_bp$Description == "response to oxidative stress")

tiff("gsea_bp_gseaplot.tiff", width=8, height=10, units="in", res=72)
gseaplot(gsea_bp, by = "all", title = gsea_bp$Description[3], geneSetID = 3)
dev.off()
~~~
{: language-r}

![gseaplot](../assets/images/gsea_bp_gseaplot.tiff)

The green line shows the running enrichment score as more genes from the sorted list are included in the analysis. The red line shows the location of the maximum enrichment score. The black lines in the Running Enrichment Score show where the members of the gene set appear in the ranked list of genes, indicating the leading edge subset.

We can get more details, including the ENTREZ id's of the genes in the leading edge 
~~~
gsea_bp[3,]
~~~
{: language-r}

~~~
                   ID                  Description setSize enrichmentScore
GO:0006979 GO:0006979 response to oxidative stress     441       0.5512628
                NES pvalue     p.adjust      qvalues rank
GO:0006979 1.721731  1e-10 1.321667e-08 1.163743e-08 4129
                             leading_edge
GO:0006979 tags=29%, list=15%, signal=25%
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       core_enrichment
GO:0006979 828271/828273/839515/822488/842405/832546/842342/834244/28720712/826447/819331/838544/822462/834223/828272/821225/820264/827765/838483/841789/831241/830416/841082/826718/832145/815160/842120/836103/839666/831075/829623/835953/829304/840158/838016/830938/827463/836876/824072/823768/839969/818490/843584/837429/818289/825820/841747/836884/838212/823994/839622/818670/831076/818588/836437/837531/823752/826616/842368/838052/818280/832302/836531/831563/819032/835756/828481/821227/822119/2745879/838752/830562/834464/829795/833957/824807/817494/833755/828413/835603/832476/837252/844177/826730/818562/828053/837894/824073/823067/817499/843660/841062/832472/831650/827831/838017/839174/818420/821193/829940/836532/842912/841843/827622/824919/834678/830637/825818/829530/828414/825154/844073/836286/826884/836025/820564/828418/837662/833779/832672/837075/834035/816500/818949/830424/820728/824493
~~~
{: .language-r}

You can also run the GSEA analysis over any gene-set of your choice. Several other enrichment tests are possible using `clusterProfiler` and you can read about them in the [clusterProfiler book](https://yulab-smu.github.io/clusterProfiler-book/).

