---
title: "53_WGCNA_ElevCrabs"
author: "Aidan Coyle"
date: "Last compiled on 2021-05-06"
output: 
  html_document:
    keep_md: true
---




## Analysis of Elevated-Temperature Gene Expression using WGCNA

We will now run WGCNA on all individual libraries of elevated-temperature crabs over days 0, 2, and 17. 

We will include transcripts from both C. bairdi and Hematodinium by using kallisto alignments to an unfiltered transcriptome (cbai_transcriptomev2.0, AKA cbaihemat_transcriptomev2.0)

Table of crabs and libraries included in analysis:

| Crab ID | Treatment Group  | Day 0 Sample ID | Day 2 Sample ID |
|---------|------------------|-----------------|-----------------|
| G       | elevated         | 173             | 272             |
| H       | elevated         | 072             | 294             |
| I       | elevated         | 127             | 280             |

Again, script is based largely on [Yaamini's script](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/analyses/WGCNA/WGCNA.md), which is based largely on [the official WGCNA tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

We will first extract the TPM (transcripts per million) counts from the kallisto libraries created earlier in the pipeline. We will then change those to logTPM counts, and then begin the WGCNA analysis

All kallisto libraries should be available within the GitHub repo. However, one datafile - our blastx table, cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt,is not. This is a BLASTx annotation of cbai_transcriptome_v2.0.fasta (known by my notation as cbai_hemat_transcriptome_v2.0.fasta). Annotation description available [here](https://robertslab.github.io/sams-notebook/2020/05/02/Transcriptome-Assembly-C.bairdi-All-RNAseq-Data-Without-Taxonomic-Filters-with-Trinity-on-Mox.html), [direct file available here](https://gannet.fish.washington.edu/Atumefaciens/20200508_cbai_diamond_blastx_transcriptome-v2.0/20200507.C_bairdi.Trinity.blastx.outfmt6), and [original transcriptome available here](https://owl.fish.washington.edu/halfshell/genomic-databank/cbai_transcriptome_v2.0.fasta), 
md5sum: ace82a75cb947574ac807d868427253c



```r
library(tidyverse)
library(WGCNA)
library(DESeq2)
```

### Step 1: Extract TPM Counts from Kallisto Libraries

This portion of the script is based largely off 21_obtaining_TPM_for_DEGs.Rmd

First, set all variables.

Later in the script, you will need to specify a soft thresholding power to choose. But otherwise, minimal changes are necessary


```r
# Path to kallisto libraries
kallisto_path <- "../output/kallisto_libraries/cbaihemat_transcriptomev2.0/"

# Libraries we want to read in to our TPM matrix
libraries <- c("173", "072", "127", "272", "294", "280")

# For each row, crab and day should correspond to the order of libraries (ex: 4th row of crabTraits should match libraries[4])
crabTraits <- data.frame("crab" = rep(c("G", "H", "I"), times = 2),
                         "day" = factor(c(rep(0, times = 3),
                                   rep(2, times = 3))))

# Create clinical data trait matrix. Same rules as above, but both crab and day are numeric. Crab A will be noted as 1, B as 2, and C as 3
crabClinicalData <- data.frame("crab" = rep(c(1, 2, 3), times = 2),
                         "day" = c(rep(0, times = 3),
                                   rep(2, times = 3)))

# Variable being examined - should match column in two data frames above
variable <- "day"

# Start and ending we want for each file and graph saved
file_start <- "../output/WGCNA_output/cbai_transcriptome_v2.0/elev_crabs_no_filter/"

# Location of blastx table
blastx_table_site <- "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt"

# Number of samples we're examining
numsamples <- 6
```

Then, we begin creating our TPM matrix for all transcripts


```r
# Create character vector with all filenames for our libraries
kallisto_files <- paste0(kallisto_path, "id", libraries, "/abundance.tsv")

# Read first kallisto file in to start data frame
TPMcounts <- read.delim(file = kallisto_files[1],
                        header = TRUE,
                        sep = "\t")
# Eliminate all columns except transcript ID and TPM
TPMcounts <- TPMcounts %>%
  select(target_id, tpm)

# Rename columns for consistency and to ID TPM counts
colnames(TPMcounts)[1:2] <- c("Transcript_ID",
                              paste0("id", libraries[1], "_TPM"))

# Loop through remaining kallisto files, performing full joins to the kallisto file we read in
for (i in 2:length(kallisto_files)){
  idnum <- str_extract(kallisto_files[i], "id[0-9]+")
  kallisto_output <- read.delim(file = kallisto_files[i],
                                header = TRUE,
                                sep = "\t")
  # Select only transcript ID and TPM (transcripts per million) columns
  kallisto_output <- kallisto_output %>%
    select(target_id, tpm)
  # Rename kallisto column names to give ID to count column
  colnames(kallisto_output)[1:2] <- c("Transcript_ID", 
                                      paste0(idnum, "_TPM"))
  # Add TPM value to table of DEGs
  # Perform full join, keeping all transcript IDs
  TPMcounts <- full_join(TPMcounts, kallisto_output, by = "Transcript_ID")
}
```

WGCNA has several recommendations when it comes to RNAseq data, available in the FAQ [here](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html). First, they suggest removing all transcripts with counts below 10 in over 90% of samples. Since we have 6 samples total, we will only remove all transcripts with counts below 10 in all samples.

They also suggest a variance-stabilizing transformation or log-transforming the counts using log2(x+1). Unlike in our trial run of Crab A, we will be able to create a DESeq object and perform a variance-stabilizing transformation.

We will change our data to fit both of these recommendations. After, we will transpose the data frame so samples are rows and transcripts are columns


```r
# Create logical matrix for whole dataframe, comparing values to 10

# Move transcript ID to rownames
TPMcounts <- TPMcounts %>%
  column_to_rownames(var = "Transcript_ID")

# Get initial dimensions of data frame
dim(TPMcounts)
```

```
## [1] 1412254       6
```

```r
# Filter out all variables with no counts greater than 10
TPMcounts <- TPMcounts %>%
  filter_all(any_vars(. > 10))

# See how many transcripts we have left
dim(TPMcounts)
```

```
## [1] 26950     6
```

```r
# Round all counts to the nearest integer
TPMcounts <- round(TPMcounts, digits = 0)

# Normalize raw counts with DESeq()
crab.dds <- DESeqDataSetFromMatrix(countData = TPMcounts,
                                   colData = crabTraits,
                                   design = as.formula(paste0("~", variable)))
```

```
## converting counts to integer mode
```

```r
crab.dds <- DESeq(crab.dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
# Perform vst on DESeq object
vsd <- getVarianceStabilizedData(crab.dds)

# Transpose dataframe to format for WGCNA
CrabExpr0 <- as.data.frame(t(vsd))

# Check dataframe was transposed correctly
dim(CrabExpr0)
```

```
## [1]     6 26950
```

We will now begin analysis with WGCNA. Our script is based largely on [Yaamini's WGCNA script](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/analyses/WGCNA/WGCNA.md), which is based largely on [the WGCNA tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

### Check for Missing Values and Outliers


```r
# Check for genes and samples with too many missing values

gsg <- goodSamplesGenes(CrabExpr0, verbose = 3)
```

```
##  Flagging genes and samples with too many missing values...
##   ..step 1
```

```r
gsg$allOK      # should return TRUE if all genes pass test
```

```
## [1] TRUE
```

```r
sampleTree <- hclust(dist(CrabExpr0), method = "average")
```


```r
path <- paste0(file_start, "ClusterDendrogram.png")
png(path)
plot(sampleTree)
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plot(sampleTree)
```

![](53_WGCNA_ElevCrabs_files/figure-html/sample_tree-1.png)<!-- -->

### Create clinical trait data matrix


```r
# Print the crabTraits matrix we made earlier
head(crabTraits)
```

```
##   crab day
## 1    G   0
## 2    H   0
## 3    I   0
## 4    G   2
## 5    H   2
## 6    I   2
```

```r
# Use same rownames as expression data to create analogous  matrix
rownames(crabTraits) <- rownames(CrabExpr0)
# Make sure it looks good
head(crabTraits)
```

```
##           crab day
## id173_TPM    G   0
## id072_TPM    H   0
## id127_TPM    I   0
## id272_TPM    G   2
## id294_TPM    H   2
## id280_TPM    I   2
```

```r
# Create a dendrogram to look at sample and trait clustering
sampleTree2 <- hclust(dist(CrabExpr0), method = "average")
traitColors <- numbers2colors(crabClinicalData, signed = FALSE)
```


```r
path <- paste0(file_start, "ClusterDendrogram_W_Colors.png")
png(path)
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels = names(crabTraits))
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels = names(crabTraits))
```

![](53_WGCNA_ElevCrabs_files/figure-html/sample_dendrogram-1.png)<!-- -->

### Network Connection and Module Detection

[Relevant tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

#### Determine soft-thresholding power


```r
# Create set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Use network topology analysis function to eval soft-thresholding power vals
sft <- pickSoftThreshold(CrabExpr0, powerVector = powers, verbose = 5)
```

```
## pickSoftThreshold: will use block size 1660.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 1660 of 26950
```

```
## Warning: executing %dopar% sequentially: no parallel backend registered
```

```
##    ..working on genes 1661 through 3320 of 26950
##    ..working on genes 3321 through 4980 of 26950
##    ..working on genes 4981 through 6640 of 26950
##    ..working on genes 6641 through 8300 of 26950
##    ..working on genes 8301 through 9960 of 26950
##    ..working on genes 9961 through 11620 of 26950
##    ..working on genes 11621 through 13280 of 26950
##    ..working on genes 13281 through 14940 of 26950
##    ..working on genes 14941 through 16600 of 26950
##    ..working on genes 16601 through 18260 of 26950
##    ..working on genes 18261 through 19920 of 26950
##    ..working on genes 19921 through 21580 of 26950
##    ..working on genes 21581 through 23240 of 26950
##    ..working on genes 23241 through 24900 of 26950
##    ..working on genes 24901 through 26560 of 26950
##    ..working on genes 26561 through 26950 of 26950
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.5220  3.370         0.8100   12700     13800  16800
## 2      2   0.1830  1.030         0.4140    8090      9060  12200
## 3      3   0.1250  0.570         0.0779    5990      6670   9330
## 4      4   0.3680  0.566         0.1890    4850      5280   7770
## 5      5   0.0848  0.214        -0.1460    4160      4320   7220
## 6      6   0.0332 -0.195        -0.2370    3700      3620   6960
## 7      7   0.1170 -0.466        -0.1330    3370      3080   6780
## 8      8   0.1940 -0.684        -0.0343    3120      2650   6660
## 9      9   0.2400 -0.847         0.0235    2930      2300   6560
## 10    10   0.2720 -0.947         0.0640    2770      2000   6490
## 11    12   0.3210 -1.170         0.1300    2540      1550   6370
## 12    14   0.3460 -1.300         0.1620    2360      1240   6290
## 13    16   0.3660 -1.360         0.1870    2230      1020   6230
## 14    18   0.3800 -1.470         0.2080    2130       939   6180
## 15    20   0.3870 -1.530         0.2150    2040       881   6140
```

#### Plot scale-free topology fit as function of soft-thresholding power


```r
path <- paste0(file_start, "ScaleIndependence.png")
png(path)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               labels = powers,
                               cex = 1,
                               col = "red")
abline(h=0.80, col = "black")
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               labels = powers,
                               cex = 1,
                               col = "red")
abline(h=0.80, col = "black")
```

![](53_WGCNA_ElevCrabs_files/figure-html/scale_free_topology_vs_soft_thresholding_power-1.png)<!-- -->

#### Plot mean connectivity as function of soft-thresholding power


```r
path <- paste0(file_start, "MeanConnectivity.png")
png(path)
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
# Add sft values
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers,
     cex = 1,
     col = "red")
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
# Add sft values
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers,
     cex = 1,
     col = "red")
```

![](53_WGCNA_ElevCrabs_files/figure-html/mean_connectivity_vs_soft_thresholding_power-1.png)<!-- -->

Typically, we would choose the lowest power below 15 that reached an R2 value of 0.8 or higher. However, no eligible soft-thresholding power reached 0.8. According to the [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html), this indicates the following:

"If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) and the mean connectivity remains relatively high (in the hundreds or above) [note: which our data does], chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."

I chose to use an unsigned network, since the direction of correlation could be quite interesting (i.e. associating up-regulated C. bairdi with down-regulated Hemat.). Therefore, the WGCNA FAQ recommends choosing a soft-thresholding power of 9, since we have fewer than 20 samples total.


#### Co-expression similarity and adjacency, and Topological Overlap Matrix (TOM)


```r
softPower <- 9
adjacency <- adjacency(CrabExpr0, power = softPower)

# Minimize noise and spurious associations by transforming adjacency into TOM
TOM <- TOMsimilarity(adjacency)
```

```
## ..connectivity..
## ..matrix multiplication (system BLAS)..
## ..normalization..
## ..done.
```

```r
#Calculate dissimilarity matrix
dissTOM <- 1 - TOM

# Clustering using TOM

# Create hierarchical clustering object
geneTree <- hclust(as.dist(dissTOM), method = "average")
```

#### Plot initial dendrogram. Dissimilarity is based on topological overlap

```r
path <- paste0(file_start, "GeneDendrogram.png")
png(path)
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)
```

![](53_WGCNA_ElevCrabs_files/figure-html/gene_dendrogram-1.png)<!-- -->


```r
# Set minimum module size, AKA num of genes that need to be in a module. Here, using WGCNA default
minModuleSize <- 30
# Cut branches of dendrogram to ID WGCNA modules
dynamicMods <- cutreeDynamic(dendro =  geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
```

```
##  ..cutHeight not given, setting it to 0.973  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
# Look at table of modules
table(dynamicMods)
```

```
## dynamicMods
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 8215 3273 2591 2161 1910 1295  837  649  642  597  457  389  316  280  274  271 
##   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
##  262  220  203  177  171  168  167  166  155  143  134  123  123  111   76   58 
##   33   34   35   36   37   38   39 
##   55   55   51   49   48   46   32
```

```r
# Convert module numbers into colors
dynamicColors <- labels2colors(dynamicMods)
```


```r
path <- paste0(file_start, "GeneDendrogramWColors.png")
png(path)
# Plot dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
```

![](53_WGCNA_ElevCrabs_files/figure-html/gene_dendrogram_module_colors-1.png)<!-- -->

#### Merge modules with similar expression profiles

Merging lets us combine modules with genes that are highly co-expressed. To do this, we can create and cluster eigengenes


```r
# Calculate eigengenes
MElist <- moduleEigengenes(CrabExpr0, colors = dynamicColors)
# Save eigengenes as new object
MEs <- MElist$eigengenes
# Calculate dissimilarity of eigengenes
MEDiss <- 1-cor(MEs)
# Create cluster object
METree <- hclust(as.dist(MEDiss), method = "average")
```


```r
path <- paste0(file_start, "ClusteredEigengenes.png")
png(path)
# Plot dendrogram of clustered eigengenes
plot(METree, main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
# ID cut height based on sample number
dynamicMergeCut(numsamples)
```

```
## [1] 0.8850638
```

```r
MEDissThres <- dynamicMergeCut(numsamples)
abline(h = MEDissThres, col = "red")
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plot(METree, main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
# ID cut height based on sample number
dynamicMergeCut(numsamples)
```

```
## [1] 0.8850638
```

```r
MEDissThres <- dynamicMergeCut(numsamples)
abline(h = MEDissThres, col = "red")
```

![](53_WGCNA_ElevCrabs_files/figure-html/clustered_eigengenes-1.png)<!-- -->



```r
merge <- mergeCloseModules(CrabExpr0, dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.885063807768968
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 39 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 5 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 4 module eigengenes in given set.
##    Calculating new MEs...
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 4 module eigengenes in given set.
```

```r
# Extract merged colors and eigengenes
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
```


```r
path <- paste0(file_start, "ClusterDendrogramOrigAndMergedEigengenes.png")
png(path)
# Plot dendrogram with original and merged eigengenes
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
```

![](53_WGCNA_ElevCrabs_files/figure-html/dendrogram_w_orig_and_merged_eigengenes-1.png)<!-- -->


```r
# Rename and save variables for subsequent analysis

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50)) # Determine color order
moduleLabels <- match(moduleColors, colorOrder)-1 # Construct numerical labels based on colors
MEs <- mergedMEs # Replace unmerged MEs
```

### Relate modules to external traits


```r
# Count the number of genes and samples
nGenes <- ncol(CrabExpr0)
nSamples <- nrow(CrabExpr0)

# Recalculate MEs with color labels, order MEs based on MEs0
MEs0 <- moduleEigengenes(CrabExpr0, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Calculate trait correlations and obtain p-values
moduleTraitCor <- cor(MEs, crabClinicalData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue
```

```
##                       crab       day
## MEdarkolivegreen 0.4395886 0.3405617
## MEturquoise      0.2672544 0.3422808
## MEgreen          0.8837890 0.1061608
## MEorangered4     0.0205042 0.7926368
```

```r
# Create text matrix for correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "")
  
# Ensure matrix has same dimensions 
dim(textMatrix) == dim(moduleTraitCor)
```

```
## logical(0)
```


```r
# Create labeled heat map of correlation values from textMatrix. Red = positive correlation, blue = negative correlation
path <- paste0(file_start, "ModuleTreatmentHeatMap.png")
png(path)
par(mar = c(4, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(crabTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-Treatment relationships"))
dev.off()
```

```
## png 
##   2
```

```r
# Plot image again, so it shows up in knitted .Rmd
par(mar = c(4, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(crabTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-Treatment relationships"))
```

![](53_WGCNA_ElevCrabs_files/figure-html/heatmap-1.png)<!-- -->

### Gene Significance and Module Membership

Module membership info is needed for a possible downstream GO-MWU analysis. 


```r
# Define "day" using information from trait matrix
day <- as.data.frame(crabClinicalData$day)
# Modify names
names(day) <- "day"
# Save module names without "ME" at beginning of each entry
modNames <- substr(names(MEs), 3, nchar(names(MEs)))

# Determine gene significance

# Obtain gene significance statistics
geneTraitSignificance <- as.data.frame(cor(CrabExpr0, day, use = "p"))
# Add column names
names(geneTraitSignificance) <- paste("GS.", names(day), sep = "")
# Confirm formatting
head(geneTraitSignificance)
```

```
##                              GS.day
## TRINITY_DN88421_c0_g1_i1  0.7726992
## TRINITY_DN88418_c0_g1_i5  0.6168431
## TRINITY_DN88415_c0_g1_i2  0.6987604
## TRINITY_DN88420_c0_g1_i1 -0.5858301
## TRINITY_DN88474_c0_g1_i2  0.6216445
## TRINITY_DN88425_c0_g1_i2  0.7037341
```

```r
# Obtain p-values for each gene significance stat
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# Add column names
names(GSPvalue) <- paste("p.GS", names(day), sep = "")
# Confirm formatting
head(GSPvalue)
```

```
##                             p.GSday
## TRINITY_DN88421_c0_g1_i1 0.07162663
## TRINITY_DN88418_c0_g1_i5 0.19208835
## TRINITY_DN88415_c0_g1_i2 0.12244992
## TRINITY_DN88420_c0_g1_i1 0.22178237
## TRINITY_DN88474_c0_g1_i2 0.18764800
## TRINITY_DN88425_c0_g1_i2 0.11865806
```

```r
# Determine module membership

# Obtain gene module membership stats
geneModuleMembership <- as.data.frame(cor(CrabExpr0, MEs, use = "p"))
# Add column names
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
# Confirm formatting
head(geneModuleMembership)
```

```
##                          MMdarkolivegreen MMturquoise    MMgreen MMorangered4
## TRINITY_DN88421_c0_g1_i1     -0.700308755  -0.4381847  0.8865821    0.3026692
## TRINITY_DN88418_c0_g1_i5     -0.587534275  -0.2893182  0.8499210    0.1747400
## TRINITY_DN88415_c0_g1_i2     -0.620714218  -0.5314422  0.9339881    0.2145396
## TRINITY_DN88420_c0_g1_i1     -0.003895447   0.9626140 -0.4599639   -0.2642649
## TRINITY_DN88474_c0_g1_i2     -0.616634452  -0.7453272  0.3999210    0.7064188
## TRINITY_DN88425_c0_g1_i2     -0.485191802  -0.3280247  0.9857573   -0.2255965
```

```r
# Obtain p-values for each module membership statistic
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# Add column names
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
# Confirm formatting
head(MMPvalue)
```

```
##                          p.MMdarkolivegreen p.MMturquoise    p.MMgreen
## TRINITY_DN88421_c0_g1_i1          0.1212639    0.38478993 0.0185659317
## TRINITY_DN88418_c0_g1_i5          0.2201060    0.57813144 0.0320953802
## TRINITY_DN88415_c0_g1_i2          0.1885050    0.27788455 0.0063925279
## TRINITY_DN88420_c0_g1_i1          0.9941569    0.00207044 0.3587106387
## TRINITY_DN88474_c0_g1_i2          0.1922823    0.08902855 0.4320996026
## TRINITY_DN88425_c0_g1_i2          0.3293221    0.52561070 0.0003028374
##                          p.MMorangered4
## TRINITY_DN88421_c0_g1_i1      0.5598598
## TRINITY_DN88418_c0_g1_i5      0.7405577
## TRINITY_DN88415_c0_g1_i2      0.6831279
## TRINITY_DN88420_c0_g1_i1      0.6128303
## TRINITY_DN88474_c0_g1_i2      0.1166330
## TRINITY_DN88425_c0_g1_i2      0.6673460
```

### Obtain gene lists - for each module, and a master list with membership and gene significance info


```r
# Save gene names as probes
probes <- names(CrabExpr0)
# Write out the gene lists for each module of interest
for (module in modNames) {
  modGenes <- (moduleColors == module) # Select module probes
  modLLIDs <- probes[modGenes] # Get gene IDs
  fileName <- paste(file_start, "GeneList-", module, ".txt", sep = "") # Assign filename for each module
  write.table(as.data.frame(modLLIDs), file = fileName, sep = "\t", row.names = FALSE, col.names = FALSE) # Write out files
}
```

#### Master list with membership and gene significance information


```r
# Import gene annotation info
crabGeneAnnot <- read.delim(blastx_table_site, header = FALSE, sep = "\t")
# Remove unnecessary columns
crabGeneAnnot <- crabGeneAnnot[, -c(3:10, 12)]
# Name columns
colnames(crabGeneAnnot) <- c("seqIDs", "Uniprot", "e-value")
# Look at column formatting
head(crabGeneAnnot)
```

```
##                     seqIDs               Uniprot  e-value
## 1 TRINITY_DN88428_c0_g1_i1 sp|Q9VVH9|SO74D_DROME  1.6e-70
## 2 TRINITY_DN88425_c3_g1_i1  sp|P11309|PIM1_HUMAN 7.0e-183
## 3 TRINITY_DN88425_c3_g1_i3  sp|P11309|PIM1_HUMAN 4.1e-188
## 4 TRINITY_DN88425_c3_g1_i2  sp|P11309|PIM1_HUMAN  3.2e-64
## 5 TRINITY_DN88463_c0_g1_i3 sp|Q9EPL4|METL9_MOUSE  1.7e-49
## 6 TRINITY_DN88429_c1_g2_i1  sp|P06603|TBA1_DROME  5.2e-16
```

```r
# If pipes in Uniprot ID column, separate to specifically get Uniprot ID, and then remove those new columns with species info
gene_ids <- dplyr::pull(crabGeneAnnot, Uniprot)
if(any(grepl("|", gene_ids, fixed = TRUE))) {
  crabGeneAnnot <- separate(data = crabGeneAnnot, col = Uniprot, into = c("sp", "Uniprot", "Species"),
                   sep = "\\|")
  crabGeneAnnot <- crabGeneAnnot[,-c(2, 4)]
}

# Also get all GO terms
```
