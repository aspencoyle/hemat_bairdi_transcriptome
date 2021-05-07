---
title: "54_WGCNA_AmbCrabs"
author: "Aidan Coyle"
date: "Last compiled on 2021-05-06"
output: 
  html_document:
    keep_md: true
---




## Analysis of Ambient-Temperature Gene Expression using WGCNA

We will now run WGCNA on all individual libraries of ambient-temperature crabs over days 0, 2, and 17. 

We will include transcripts from only C. bairdi by using kallisto alignments to a transcriptome filtered to only include sequences BLASTed against a _Chionoecetes opilio_ genome.

Table of crabs and libraries included in analysis:

| Crab ID | Treatment Group | Day 0 Sample ID | Day 2 Sample ID | Day 17 Sample ID |
|---------|-----------------|-----------------|-----------------|------------------|
| A       | ambient         | 178             | 359             | 463              |
| B       | ambient         | 118             | 349             | 481              |
| C       | ambient         | 132             | 334             | 485              |

Again, script is based largely on [Yaamini's script](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/analyses/WGCNA/WGCNA.md), which is based largely on [the official WGCNA tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

We will first extract the TPM (transcripts per million) counts from the kallisto libraries created earlier in the pipeline. We will then change those to logTPM counts, and then begin the WGCNA analysis

All kallisto libraries should be available within the GitHub repo. However, one datafile - our blastx table, cbai_hemat_diamond_blastx_table_transcriptome_v4.0.txt,is not. This is a BLASTx annotation of cbai_transcriptome_v4.0.fasta.

Annotation description available [here](https://robertslab.github.io/resources/Genomic-Resources/), [direct file available here](https://gannet.fish.washington.edu/Atumefaciens/20210318_cbai_diamond_blastx_transcriptome-v4.0/cbai_transcriptome_v4.0.blastx.outfmt6), and [original transcriptome available here](https://gannet.fish.washington.edu/Atumefaciens/20210317_cbai_trinity_RNAseq_transcriptome-v4.0/cbai_transcriptome_v4.0.fasta_trinity_out_dir/cbai_transcriptome_v4.0.fasta), 
md5sum: 8fd2ab9c27e59653fcb4f32574be1f61



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
kallisto_path <- "../output/kallisto_libraries/cbai_transcriptomev4.0/"

# Libraries we want to read in to our TPM matrix
libraries <- c("178", "118", "132", "359", "349", "334", "463", "481", "485")

# For each row, crab and day should correspond to the order of libraries (ex: 4th row of crabTraits should match libraries[4])
crabTraits <- data.frame("crab" = rep(c("A", "B", "C"), times = 3),
                         "day" = factor(c(rep(0, times = 3),
                                   rep(2, times = 3),
                                   rep(17, times = 3))))

# Create clinical data trait matrix. Same rules as above, but both crab and day are numeric. Crab A will be noted as 1, B as 2, and C as 3
crabClinicalData <- data.frame("crab" = rep(c(1, 2, 3), times = 3),
                         "day" = c(rep(0, times = 3),
                                   rep(2, times = 3),
                                   rep(17, times = 3)))

# Variable being examined - should match column in two data frames above
variable <- "day"

# Start and ending we want for each file and graph saved
file_start <- "../output/WGCNA_output/cbai_transcriptome_v4.0/amb_crabs_no_filter/"

# Location of blastx table
blastx_table_site <- "../data/cbai_blastx_table_transcriptome_v4.0.txt"

# Number of samples we're examining
numsamples <- 9
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

WGCNA has several recommendations when it comes to RNAseq data, available in the FAQ [here](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html). First, they suggest removing all transcripts with counts below 10 in over 90% of samples. Since we have 9 samples total, we will only remove all transcripts with counts above 10 in 1 (or zero) sample

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
## [1] 88302     9
```

```r
# Get count of values >= 10 for each row
row_sums_tpm <- rowSums(TPMcounts >= 10)
# Filter rows without multiple rows with counts > 10
TPMcounts <- TPMcounts %>%
  filter(row_sums_tpm > 1)

# See how many transcripts we have left
dim(TPMcounts)
```

```
## [1] 7173    9
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
## [1]    9 7173
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/sample_tree-1.png)<!-- -->

### Create clinical trait data matrix


```r
# Print the crabTraits matrix we made earlier
head(crabTraits)
```

```
##   crab day
## 1    A   0
## 2    B   0
## 3    C   0
## 4    A   2
## 5    B   2
## 6    C   2
```

```r
# Use same rownames as expression data to create analogous  matrix
rownames(crabTraits) <- rownames(CrabExpr0)
# Make sure it looks good
head(crabTraits)
```

```
##           crab day
## id178_TPM    A   0
## id118_TPM    B   0
## id132_TPM    C   0
## id359_TPM    A   2
## id349_TPM    B   2
## id334_TPM    C   2
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/sample_dendrogram-1.png)<!-- -->

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
## pickSoftThreshold: will use block size 6237.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 6237 of 7173
```

```
## Warning: executing %dopar% sequentially: no parallel backend registered
```

```
##    ..working on genes 6238 through 7173 of 7173
##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.8410  1.5900          0.905    3990      4370   5210
## 2      2   0.5660  0.5190          0.817    2760      3040   4250
## 3      3   0.1100  0.1410          0.605    2090      2250   3640
## 4      4   0.0208 -0.0548          0.535    1680      1710   3220
## 5      5   0.1860 -0.1840          0.596    1390      1330   2890
## 6      6   0.3540 -0.2720          0.684    1180      1050   2630
## 7      7   0.4680 -0.3280          0.723    1030       849   2420
## 8      8   0.5310 -0.3770          0.741     902       694   2250
## 9      9   0.6070 -0.4210          0.783     803       582   2110
## 10    10   0.6470 -0.4570          0.774     722       497   1980
## 11    12   0.6910 -0.5270          0.765     597       369   1780
## 12    14   0.7170 -0.5760          0.745     505       283   1610
## 13    16   0.7360 -0.6140          0.735     435       220   1470
## 14    18   0.7480 -0.6530          0.714     380       173   1360
## 15    20   0.7610 -0.6800          0.715     335       139   1250
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/scale_free_topology_vs_soft_thresholding_power-1.png)<!-- -->

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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/mean_connectivity_vs_soft_thresholding_power-1.png)<!-- -->

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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/gene_dendrogram-1.png)<!-- -->


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
##  ..cutHeight not given, setting it to 0.988  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
# Look at table of modules
table(dynamicMods)
```

```
## dynamicMods
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 2225 1790  412  357  330  258  253  236  163  149  145   92   87   75   74   74 
##   17   18   19   20   21   22   23   24   25 
##   69   67   56   56   56   42   40   34   33
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/gene_dendrogram_module_colors-1.png)<!-- -->

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
## [1] 0.5278481
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
## [1] 0.5278481
```

```r
MEDissThres <- dynamicMergeCut(numsamples)
abline(h = MEDissThres, col = "red")
```

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/clustered_eigengenes-1.png)<!-- -->



```r
merge <- mergeCloseModules(CrabExpr0, dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.527848068913797
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 25 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 8 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 6 module eigengenes in given set.
##    Calculating new MEs...
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 6 module eigengenes in given set.
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/dendrogram_w_orig_and_merged_eigengenes-1.png)<!-- -->


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
##                    crab       day
## MEtan         0.3373140 0.1241068
## MEgreen       0.5180555 0.2855814
## MEgreenyellow 0.9935509 0.2094550
## MEpurple      0.5030387 0.8983927
## MEblack       0.5577774 0.1864723
## MEsalmon      0.1217332 0.3717224
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

![](54_WGCNA_cbaiv4.0_AmbCrabs_files/figure-html/heatmap-1.png)<!-- -->

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
## TRINITY_DN21478_c0_g1_i1  0.2998367
## TRINITY_DN21505_c1_g1_i1 -0.5521966
## TRINITY_DN21444_c0_g2_i1  0.5024371
## TRINITY_DN21517_c0_g2_i1 -0.3611415
## TRINITY_DN21517_c0_g1_i1 -0.3907039
## TRINITY_DN21490_c0_g1_i1  0.4220838
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
##                            p.GSday
## TRINITY_DN21478_c0_g1_i1 0.4331081
## TRINITY_DN21505_c1_g1_i1 0.1231665
## TRINITY_DN21444_c0_g2_i1 0.1680619
## TRINITY_DN21517_c0_g2_i1 0.3396310
## TRINITY_DN21517_c0_g1_i1 0.2984872
## TRINITY_DN21490_c0_g1_i1 0.2577817
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
##                               MMtan    MMgreen MMgreenyellow    MMpurple
## TRINITY_DN21478_c0_g1_i1  0.2196367 -0.2706075    -0.1362901  0.31352921
## TRINITY_DN21505_c1_g1_i1  0.4178235  0.9656494     0.4394392  0.02425170
## TRINITY_DN21444_c0_g2_i1  0.1948624 -0.5754775    -0.4060458 -0.11770122
## TRINITY_DN21517_c0_g2_i1  0.3296174  0.9912441     0.4606235 -0.15586610
## TRINITY_DN21517_c0_g1_i1  0.3763382  0.9925252     0.4246795 -0.02941401
## TRINITY_DN21490_c0_g1_i1 -0.3506487 -0.8544990    -0.4675833 -0.03283013
##                             MMblack    MMsalmon
## TRINITY_DN21478_c0_g1_i1  0.7041579  0.02626978
## TRINITY_DN21505_c1_g1_i1 -0.8182067 -0.38394597
## TRINITY_DN21444_c0_g2_i1  0.9610914  0.22674883
## TRINITY_DN21517_c0_g2_i1 -0.7846982 -0.56560367
## TRINITY_DN21517_c0_g1_i1 -0.7543732 -0.49235849
## TRINITY_DN21490_c0_g1_i1  0.8250409  0.55592027
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
##                            p.MMtan    p.MMgreen p.MMgreenyellow p.MMpurple
## TRINITY_DN21478_c0_g1_i1 0.5701588 4.812650e-01       0.7266139  0.4113182
## TRINITY_DN21505_c1_g1_i1 0.2631251 2.391730e-05       0.2366173  0.9506191
## TRINITY_DN21444_c0_g2_i1 0.6153771 1.049418e-01       0.2781976  0.7629720
## TRINITY_DN21517_c0_g2_i1 0.3863720 2.050646e-07       0.2121086  0.6888305
## TRINITY_DN21517_c0_g1_i1 0.3181482 1.180232e-07       0.2545547  0.9401215
## TRINITY_DN21490_c0_g1_i1 0.3548670 3.345634e-03       0.2043772  0.9331790
##                             p.MMblack p.MMsalmon
## TRINITY_DN21478_c0_g1_i1 3.421858e-02  0.9465144
## TRINITY_DN21505_c1_g1_i1 7.025836e-03  0.3076570
## TRINITY_DN21444_c0_g2_i1 3.682579e-05  0.5574077
## TRINITY_DN21517_c0_g2_i1 1.226363e-02  0.1124572
## TRINITY_DN21517_c0_g1_i1 1.883524e-02  0.1781506
## TRINITY_DN21490_c0_g1_i1 6.187779e-03  0.1201335
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
##                     seqIDs              Uniprot e-value
## 1 TRINITY_DN21515_c0_g2_i1 sp|Q9XYN1|INX2_SCHAM 1.3e-09
## 2 TRINITY_DN21515_c0_g1_i1 sp|Q9XYN1|INX2_SCHAM 6.5e-41
## 3 TRINITY_DN21478_c0_g1_i1 sp|Q6GNG3|TMX3_XENLA 5.2e-15
## 4 TRINITY_DN21478_c0_g1_i2 sp|Q6GNG3|TMX3_XENLA 2.0e-39
## 5 TRINITY_DN21438_c0_g1_i1 sp|P45594|CADF_DROME 1.0e-33
## 6 TRINITY_DN21438_c1_g1_i1 sp|P45594|CADF_DROME 5.9e-14
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
