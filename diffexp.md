---
layout: page
title: Differential expression
nav_order: 4
parent: Home
mathjax: true
---

# 1. Introduction

Differential expression analysis is the process of identifying the genes that are significantly affected by the experimental design. In the example study that we use, Arabidopsis plants were infected or not by a pathogenic bacteria called _Pseudomonas syringae_ DC3000. One comparison of interest could be to identify genes whose expression is affected by the infection with this pathogenic bacteria.

Here, we will see how to perform a simple one-condition experimental comparison using `DESeq2`. We will compare the transcriptome of Arabidopsis in response to infection by the leaf pathogenic bacteria _Pseudomonas syringae_ DC3000 after 7 days (7 dpi). This will yield a table containing genes $$log_{2}$$ fold change and their corrected p-values. We will also see how to create a few typical representations classically used to display RNA-seq results such as volcano plots and heatmaps. 

> ## Important note
> For differential expression analysis, you should use the __raw__ counts and __not__ the scaled counts. 
> As the DESeq2 model fit requires raw counts (integers), make sure that you use the `raw_counts.csv` file. 
{: .callout}

# 2. Differential expression analysis

We will do all our analysis in a folder, so we can organize our files. Let's first create a directory using the shell. Then let's make sure we tell `R` where the packages are located, and use environmental modules to load R

~~~
# create a new directory
mkdir rna_seq

# move into that directory
cd rna_seq

# let's tell R where our R packages are located. If you edited the ~/.Renviron
# file during day 2, then you should be fine. If not, then you can do
export R_LIBS=/standard/bims6000/R

# let's make sure R is available for us to use
module load goolf R/4.4.1

# let's load R
R
~~~
{: .language-bash}

Now let's make sure we are pointing to the directory that contains all the R libraries we want to use. We can check which locations are being used by R to look for libraries by typing the following in the R console

~~~
.libPaths()
~~~
{: .language-r}

You should see an output similar to the following with `/sfs/ceph/standard/bims6000/R` being the first option where R looks for libraries.

~~~
[1] "/sfs/ceph/standard/bims6000/R"                                                                              
[2] "/sfs/gpfs/tardis/home/nelle/R/goolf/4.4"                                                                    
[3] "/sfs/gpfs/tardis/applications/202506/software/standard/mpi/gcc/11.4.0/openmpi/4.1.4/R/4.4.1/lib64/R/library"
~~~
{: .output}


## 2.1 Creating the DESeqDataSet object

Since we do not want to work on all comparisons, we will filter out the samples and conditions that we do not need. Only the mock growth and the _P. syringae_ infected condition will remain.  

~~~
# Import libraries
library("DESeq2")
library("tidyverse")

# import the samples to conditions correspodence
xp_design <- read.csv("/standard/bims6000/data/afternoon/samples_to_conditions.csv",              
                      header = TRUE, 
                      stringsAsFactors = FALSE, 
                      colClasses = rep("character",4))

# filter design file to keep only "mock" and the "infected P. syringae at 7 dpi" conditions.
xp_design_mock_vs_infected <- xp_design %>% 
                              filter(growth == "MgCl2" & dpi == "7")
~~~
{: .language-r}

We then import the gene counting values and call it `raw_counts`. The gene names have to be changed to the names of the rows of the table for compatibility with `DESeq2`. This is done using the `column_to_rownames()` function from the `tibble` package (contained in `tidyverse` suite of packages).

~~~
# Import the gene raw counts
raw_counts <- read.csv("/standard/bims6000/data/afternoon/raw_counts.csv", 
                       header = TRUE, 
                       stringsAsFactors = FALSE) %>% 
              column_to_rownames("Geneid")


# reorder counts columns according to the complete list of samples 
raw_counts <- raw_counts[ , xp_design$sample]
~~~
{: .language-r}

We will now filter both the `raw_counts` and `xp_design` objects to keep a one-factor comparison and investigate the leaf transcriptome
of Arabidopsis plants whose seeds were MgCl2 treated and whose plants were infected or not with Pseudomonas syringae DC3000 at 7 dpi.

The corresponding code is available below.

~~~
# Filter count file accordingly (to keep only samples present in the filtered xp_design file)
raw_counts_filtered <- raw_counts[, colnames(raw_counts) %in% xp_design_mock_vs_infected$sample]

## Creation of the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)
~~~
{: .language-r}

You can have a glimpse at the DESeqDataSet `dds` object that you have created. It gives some useful information already. 

~~~
dds
~~~
{: .language-r}

~~~
class: DESeqDataSet 
dim: 33768 8 
metadata(1): version
assays(1): counts
rownames(33768): AT1G01010 AT1G01020 ... ATMG01400 ATMG01410
rowData names(0):
colnames(8): ERR1406305 ERR1406306 ... ERR1406265 ERR1406266
colData names(4): sample growth infected dpi
~~~
{: .output}

<br>

> ## Important note on factor levels
> It is important to make sure that levels are properly ordered so we are indeed using the _mock_ group as our reference level. A positive gene fold change means that the gene is upregulated in the _P. syringae_ condition relatively to the _mock_ condition.  
{: .callout}

Please consult [the dedicated section of the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#factorlevels) on factor levels. 

One way to see how levels are interpreted within the DESeqDataSet object is to display the factor levels. 
~~~
dds$infected
~~~
{: .language-r}

~~~
[1] mock  mock  mock  mock  Pseudomonas_syringae_DC3000
[6] Pseudomonas_syringae_DC3000 Pseudomonas_syringae_DC3000 Pseudomonas_syringae_DC3000
Levels: mock Pseudomonas_syringae_DC3000
~~~
{: .output}

This shows that the _mock_ level comes first before the _Pseudomonas_syringae_DC3000_ level. If this is not correct, you can change it following [the dedicated section of the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#factorlevels) on factor levels. 


## 2.2 Running the DE analysis

Differential gene expression analysis will consist of simply two lines of code:
1. The first will call the `DESeq` function on a `DESeqDataSet` object that you've just created under the name `dds`. It will be returned under the same `R` object name `dds`.
2. Then, results are extracted using the `results` function on the `dds` object and results will be extracted as a table under the name `res` (short for results). 

~~~
dds <- DESeq(dds)
~~~
{: .language-r}


~~~
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
~~~
{: .output}


~~~
res <- results(dds)

# have a peek at the DESeqResults object 
res
~~~
{: .language-r}

The theory beyond DESeq2 differential gene expression analysis is beyond this course but nicely explained [within the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory). 

> ## Beware of factor levels
> 
> If you do not supply any values to the contrast argument of the `DESeq` function, it will use the first value of the design variable from the design file.
> 
> In our case, we will perform a differential expression analysis between `mock` and `Pseudomonas_syringae_DC3000`. 
> 1. Which of these two is going to be used as the reference level?
> 2. How would you interpret a positive log2 fold change for a given gene?
>
> > ## Solution
> > 1. The `mock` condition is going to be used as the reference level since _m_ from `mock` comes before `P` from `Pseudomonas_syringae_DC3000`.
> > 2. A positive log2 fold change for a gene would mean that this gene is more abundant in `Pseudomonas_syringae_DC3000` than in the `mock` condition.
> {: .solution} 
{: .challenge}

The complete explanation comes from the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis):
> Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.

A possible preferred way is to specify the comparison of interest explicitly. We are going to name this new result object `all_genes_results` and compare it with the previous one called `res`.

~~~
all_genes_results <- results(dds, contrast = c("infected",                      # name of the factor
                                  "Pseudomonas_syringae_DC3000",    # name of the numerator level for fold change
                                  "mock"))                          # name of the denominator level    

~~~
{: .language-r}

If we now compare the `res` and `all_genes_results` DESeqResults objects, they should be exactly the same and return a `TRUE` value.
~~~
all_equal(res, as.data.frame(all_genes_results))
~~~
{: .language-r}

If not, that means that you should check your factor ordering. 

## 2.3 Extracting the table of differential genes 

We can now have a look at the result table that contains all information on _all_ genes (p-value, fold changes, etc).  

Let's take a peek at the first lines.
~~~
head(all_genes_results)                
~~~
{: .language-r}

~~~
log2 fold change (MLE): infected Pseudomonas_syringae_DC3000 vs mock 
Wald test p-value: infected Pseudomonas syringae DC3000 vs mock 
DataFrame with 6 rows and 6 columns
           baseMean log2FoldChange     lfcSE      stat      pvalue       padj
          <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
AT1G01010   87.4203      0.3672806  0.211702  1.734897 0.082759060 0.18722371
AT1G01020  477.1530      0.2663723  0.107898  2.468749 0.013558621 0.04572769
AT1G03987   14.6179      1.4707140  0.462673  3.178735 0.001479191 0.00740164
AT1G01030  194.0951      0.9166233  0.276959  3.309597 0.000934304 0.00506641
AT1G03993  175.9825     -0.1084689  0.142106 -0.763293 0.445288450 0.61400086
AT1G01040 1761.9499     -0.0519691  0.076330 -0.680848 0.495967753 0.65803411
~~~
{: .output}

<br>


> ## Question
> 1. What is the biological meaning of a **log2** fold change equal to 1 for gene X?
> 2. What is the biological meaning of a **log2** fold change equal to -1?
> 
> > ## Solution
> > 1. A **log2** equal to 1 means that gene X has a higher expression (two-fold) in the DC3000 infected condition compared to the mock condition. 
> > 2. A **log2** equal to -1 means that gene X has a smaller expression (0.5) in the DC3000 infected condition.   
> >  
> > {: .language-r}
> {: .solution}
{: .challenge}

<br>

Some explanations about this output:
> The results table when printed will provide the information about the comparison, e.g. "log2 fold change (MAP): condition treated vs untreated", meaning that the estimates are of log2(treated / untreated), as would be returned by contrast=c("condition","treated","untreated"). 

So in our case, since we specified `contrast = c("infected", "Pseudomonas_syringae_DC3000", "mock")`, the `log2FoldChange` will return the log2(Pseudomonas syringae DC3000 / mock).  

Additional information on the DESeqResult columns is available using the `mcols` function. 
~~~
mcols(all_genes_results)
~~~
{: .language-r}

This will indicate a few useful _metadata_ information about our results:

~~~
DataFrame with 6 rows and 2 columns
                       type                                                          description
                <character>                                                          <character>
baseMean       intermediate                            mean of normalized counts for all samples
log2FoldChange      results log2 fold change (MLE): infected Pseudomonas_syringae_DC3000 vs mock
lfcSE               results         standard error: infected Pseudomonas syringae DC3000 vs mock
stat                results         Wald statistic: infected Pseudomonas syringae DC3000 vs mock
pvalue              results      Wald test p-value: infected Pseudomonas syringae DC3000 vs mock
padj                results                                                 BH adjusted p-values

~~~
{: .output}


## 2.4 False discovery rates

When you perform thousands of statistical tests (one for each gene), you will by chance call genes differentially expressed while they are not (false positives). You can control for this by applying certain statistical procedures called _multiple hypothesis test correction_.  The selected &alpha; threshold controls for type I error rate: rejecting the _null_ hypothesis (H<sub>0</sub> no difference) and therefore affirming that there is a gene expression difference between conditions while there aren't any. This &alpha; value is often set at
at &alpha; = 0.01 (1%) or &alpha; = 0.001 (0.1%) in RNA-seq analyses.
 
We can count the number of genes that are differentially regulated at a certain  &alpha; level. 
~~~
library(dplyr)

# threshold of p = 0.01
all_genes_results %>% 
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  dim()

# threshold of p = 0.001
all_genes_results %>% 
  as.data.frame() %>% 
  filter(padj < 0.001) %>% 
  dim()
~~~
{: .language-r}

You should obtain __4979__ differentially expressed genes at 0.01 and __3249__ at 0.001 which are quite important numbers: indeed, it corresponds to respectively \~15% and \~10% of the whole number transcriptome (total number of mRNA is 33,768).    

> ## Histogram p-values
> This [blog post](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/) explains in detail what you can expect from each p-value distribution profile. 

## Extracting the table of differential genes
Ok, here's the moment you've been waiting for. How can I extract a nicely filtered final table of differential genes? Here it is!

~~~
diff_genes <- all_genes_results %>% 
              as.data.frame() %>% 
              rownames_to_column("genes") %>% 
              filter(padj < 0.01) %>% 
              arrange(desc(log2FoldChange), desc(padj))
head(diff_genes)
~~~
{: .language-r}


> ## Choosing thresholds
> Getting a list of differentially expressed genes means that you need to choose an __absolute__ threshold for the log2 fold change (column `log2FoldChange`) and the adjusted p-value (column `_padj_`). Therefore you can make different list of differential genes based on your selected thresholds. It is common to choose a log2 fold change threshold of |1| or |2| and an adjusted p-value of 0.01 for instance. 
{: .callout}

You could write this file on your disk with `write.csv()` for instance to save a comma-separated text file containing your results. Ideally, you should save the data-frame `all_genes_results` on your disk as well.

~~~
write.csv(all_genes_results, "all_genes.csv")
~~~
{: language-r}

# 3. Volcano plot
For each gene, this plot shows the gene fold change on the x-axis against the p-value plotted on the y-axis. 

Here, we make use of a library called _EnhancedVolcano_ which is available through [Bioconductor](http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) and described extensively on its [own GitHub page](https://github.com/kevinblighe/EnhancedVolcano).

First, we are going to "shrink" the $$\log2$$ fold changes to remove the noise associated with fold changes coming from genes with low count levels. Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. This helps to get more meaningful log2 fold changes for all genes independently of their expression level.

~~~
library("apeglm")

resLFC <- lfcShrink(dds = dds, 
                  res = all_genes_results,
                  type = "normal",
                  coef = "infected_Pseudomonas_syringae_DC3000_vs_mock") # name or number of the coefficient (LFC) to shrink
~~~
{: .language-r}

To see what coefficients can be extracted, type: 
~~~
resultsNames(dds)
~~~
{: .language-r}

~~~
[1] "Intercept"                                   
[2] "infected_Pseudomonas_syringae_DC3000_vs_mock"
~~~
{: .output}

We can build the Volcano plot rapidly without much customization. 
~~~
# load the library if not done yet
library("EnhancedVolcano")

# The main function is named after the package
tiff("volcano_plot.tiff", width=7, height=7, units="in", res=100)

EnhancedVolcano(toptable = resLFC,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(resLFC)
                )

dev.off()
~~~
{: .language-r}

<img src="../assets/images/volcano_plot.tiff" width="800px" alt="default volcano plot" >

Alternatively, the plot can be heavily customized to become a publication-grade figure.  
~~~
tiff("publication_ready_volcano_plot.tiff", width=7, height=7, units="in", res=200)

EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 1e-06,
                FCcutoff = 2, 
                title = "Pseudomonas syringae DC3000 versus mock \n (fold change cutoff = 2, p-value cutoff = 1e-06)",
                legendPosition = "bottom",
                legendLabSize = 10,
                legendLabels=c(
                  'Not sig.',
                  'Log2 fold-change',
                  'p-value',
                  'p-value & Log2 fold change')
                )

dev.off()
~~~
{: .language-r}
<img src="../assets/images/publication_ready_volcano_plot.tiff" width="800px" alt="customized volcano plot" >

# 4. Heatmap
Heatmap is a representation where values are represented on a color scale. It is usually one of the classic figures part of a transcriptomic study. 
One can also cluster samples and genes to identify groups of genes that show a coordinated behaviour. Let's build a nice looking heatmap to display our differential genes one step at a time.  

We are going to make use of a library called `pheatmap`. Here is a minimal example (`mtcars` is a dataset that comes included with R).
~~~
library(pheatmap)
df <- scale(mtcars)

tiff("mtcars_heatmap.tiff")
pheatmap(df)
dev.off()
~~~
{: .language-r}

<img src="../assets/images/mtcars_heatmap.tiff" alt="basic heatmap" height="400px">

> ## Troubleshooting
> If you have issues where your heatmap plot is not being shown, run `dev.off()` and try to plot again. It should solve your issue. 
{: .callout}


## 4.1 Function to scale the raw counts

Let's get the counts normalized by DESeq2 and look at the first few lines

~~~
normalized_counts <- counts(dds, normalized=TRUE)

head(normalized_counts)
~~~
{: .language-r}

~~~
          ERR1406305  ERR1406306  ERR1406307 ERR1406308 ERR1406263 ERR1406264 ERR1406265 ERR1406266
AT1G01010   85.83575   90.910197   69.891325   59.41828   74.15774  114.25728  106.48797   98.40356
AT1G01020  452.20786  456.549010  398.076675  426.22715  511.29812  588.94434  454.57769  529.34329
AT1G03987   13.95703    4.995066    5.064589    6.33795   16.91317   20.77405   21.75561   27.14581
AT1G01030  153.52736  168.833223  118.511376   97.44598  183.44284  404.05529  229.00639  197.93820
AT1G03993  174.46291  189.812499  190.428537  176.67036  148.31549  162.03760  158.01441  208.11788
AT1G01040 1811.62285 1800.221699 1874.910751 1689.85596 1592.43995 1698.27864 1745.02871 1883.24056
~~~
{: .output}


## 4.2 First version

~~~
normalised_counts_only_diff_genes <- normalized_counts %>%
                                     as_tibble(rownames="genes") %>%
                                     filter(genes %in% diff_genes$genes) %>%
                                     column_to_rownames("genes")
~~~
{: .language-r}

We indeed find that we have **4979 genes** (rows, p < 0.01) and **8 samples** (columns) which corresponds to the number of differential genes identified previously between Mock and DC3000 infected conditions at 7 dpi and with a MgCl2 seed coating. You can also use `head()` to show the first lines of this table. 

Let's plot our first version of the heatmap. 
~~~
png("normalized_count_heatmap.png")

pheatmap(normalised_counts_only_diff_genes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none",
         show_rownames = FALSE, 
         show_colnames = TRUE)

dev.off()
~~~
{: .language-r}

We have removed the genes names with `show_rownames = FALSE` since they are not readable anymore for such a high number of genes.


<img src="../assets/images/normalized_count_heatmap.png" alt="first heatmap version" height="400px">


Well....not very useful right?

> ## Question
> Do you have an idea of how to improve this heatmap?  
> {: .language-r}
> > ## Solution
> > The scale on which gene counts are represented is the (main) issue here.   
> > There are a lot of genes for which the number of counts are very low. 
> > > {: .solution}
{: .challenge}

## 4.2 Second version with scaling 

When creating a heatmap, it is vital to control how scaling is performed. We can perform a Z-score calculation for each gene so that $$Z = {x - \mu \over \sigma}$$   where $$x$$ is an individual gene count inside a given sample, $$\mu$$ the row mean of for that gene across all samples and $$\sigma$$ its standard deviation. We can specify `scale = "row"` to the `pheatmap()` function to perform row scaling since gene expression levels will become comparable. So let's see if this scaling improves our heatmap?
~~~
png("scaled_count_heatmap.png")

pheatmap(normalised_counts_only_diff_genes, 
         scale = "row",  
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = TRUE)

dev.off()
~~~
{: .language-r}


After applying the scaling procedure, the gene expression levels become more comparable. Still, this heatmap isn't really useful so far. 

<img src="../assets/images/scaled_count_heatmap.png" width="400px" alt="second heatmap (scaled)"  >

## 4.3 Third version with genes and samples grouped by profiles
One interesting feature of the heatmap visualisation is the ability to group genes and samples by their expression profile. Let's see how this heatmap looks with both gene and sample clustering.

~~~
png("clustered_count_heatmap.png")

pheatmap(normalised_counts_only_diff_genes,
         scale = "row",  
         cluster_rows = TRUE,                      
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         main = "Clustering on")

dev.off()
~~~
{: .language-r}

<img src="../assets/images/clustered_count_heatmap.png" alt="third heatmap version (clustered)" height="400px">

This is getting easier to read. Genes with similar profiles that distinguish different samples can be easily visualised. 

> ## Question
> Do you know how this gene and sample clustering was done? How can you find this out?
> > ## Solution
> > Check in the help page related to the `pheatmap` function (type `?pheatmap`) inside R. 
> > By default, the clustering distance is **euclidean** for both rows (genes) and columns (samples). The clustering_method is **complete**.
> {: .solution}
{: .challenge}

You can change this default behavior easily and try other clustering methods (see `?hclust` for supported methods).

# References
* [Kamil Slowikoski blog post about heatmap](https://slowkow.com/notes/pheatmap-tutorial/)
* Z-score calculations: [link 1](https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/z-score/) and [link 2](https://www.datatechnotes.com/2018/02/z-score-with-r.html).
* [Type I and type II error rates in gene expression studies](https://www.ncbi.nlm.nih.gov/pubmed/28637422)
* [p-value histograms explained](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/)

