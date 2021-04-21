---
title: "53_WGCNA_AmbCrabs"
author: "Aidan Coyle"
date: "Last compiled on 2021-04-20"
output: 
  html_document:
    keep_md: true
---




## Analysis of Ambient-Temperature Gene Expression using WGCNA

Now that we've completed a run through WGCNA (well, as much as possible with one crab) with a single crab (Crab A), we will scale up our run to include all individual libraries of ambient-temperature crabs over days 0, 2, and 17. 

We will include transcripts from both C. bairdi and Hematodinium by using kallisto alignments to an unfiltered transcriptome (cbai_transcriptomev2.0, AKA cbaihemat_transcriptomev2.0)

Table of crabs and libraries included in analysis:

| Crab ID | Treatment Group | Day 0 Sample ID | Day 2 Sample ID | Day 17 Sample ID |
|---------|-----------------|-----------------|-----------------|------------------|
| A       | ambient         | 178             | 359             | 463              |
| B       | ambient         | 118             | 349             | 481              |
| C       | ambient         | 132             | 334             | 485              |

Again, script is based largely on [Yaamini's script](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/analyses/WGCNA/WGCNA.md), which is based largely on [the official WGCNA tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)

We will first extract the TPM (transcripts per million) counts from the kallisto libraries created earlier in the pipeline. We will then change those to logTPM counts, and then begin the WGCNA analysis


```r
library(tidyverse)
library(WGCNA)
library(DESeq2)
```

### Step 1: Extract TPM Counts from Kallisto Libraries

This portion of the script is based largely off 21_obtaining_TPM_for_DEGs.Rmd

First, set all variables


```r
# Path to kallisto libraries
kallisto_path <- "../output/kallisto_libraries/cbaihemat_transcriptomev2.0/"

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

# Start and ending we want for each file and graph saved
file_start <- "../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/"
file_ending <- "OnlyCtsOver30"
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

WGCNA has several recommendations when it comes to RNAseq data, available in the FAQ [here](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html). First, they suggest removing all transcripts with counts below 10 in over 90% of samples. Since we have 9 samples total, we will only remove all transcripts with counts below 10 in all samples.

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
## [1] 1412254       9
```

```r
# Filter out all variables with no counts greater than 80. Should be 10, but testing if this works
TPMcounts <- TPMcounts %>%
  filter_all(any_vars(. > 30))

# See how many transcripts we have left
dim(TPMcounts)
```

```
## [1] 10719     9
```

```r
# Round all counts to the nearest integer
TPMcounts <- round(TPMcounts, digits = 0)

# Normalize raw counts with DESeq()
crab.dds <- DESeqDataSetFromMatrix(countData = TPMcounts,
                                   colData = crabTraits,
                                   design = ~day)
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
## [1]     9 10719
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
plot(sampleTree)
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30sample_tree-1.png)<!-- -->

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
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels = names(crabTraits))
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30sample_dendrogram-1.png)<!-- -->

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
## pickSoftThreshold: will use block size 4173.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 4173 of 10719
```

```
## Warning: executing %dopar% sequentially: no parallel backend registered
```

```
##    ..working on genes 4174 through 8346 of 10719
##    ..working on genes 8347 through 10719 of 10719
##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
## 1      1  0.87600  1.8300        0.84400    5510      5540   7420
## 2      2  0.64400  0.6450        0.54300    3620      3510   5740
## 3      3  0.21300  0.2530       -0.00883    2660      2410   4720
## 4      4  0.00657  0.0379       -0.24400    2090      1810   4030
## 5      5  0.07930 -0.0976       -0.18400    1720      1600   3530
## 6      6  0.27300 -0.1890        0.07030    1460      1480   3160
## 7      7  0.41500 -0.2610        0.25700    1270      1290   2870
## 8      8  0.41200 -0.3180        0.24400    1120      1100   2640
## 9      9  0.45100 -0.3750        0.29600    1010       944   2450
## 10    10  0.53800 -0.4250        0.41800     914       880   2280
## 11    12  0.51300 -0.4710        0.37700     773       783   2010
## 12    14  0.51600 -0.5360        0.38400     672       637   1810
## 13    16  0.50400 -0.5960        0.36700     595       520   1640
## 14    18  0.49700 -0.6260        0.35600     536       437   1490
## 15    20  0.48800 -0.6330        0.34200     489       373   1370
```

#### Plot scale-free topology fit as function of soft-thresholding power


```r
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               labels = powers,
                               cex = 1,
                               col = "red")
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30scale_free_topology_vs_soft_thresholding_power-1.png)<!-- -->

#### Plot mean connectivity as function of soft-thresholding power


```r
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

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30mean_connectivity_vs_soft_thresholding_power-1.png)<!-- -->

Typically, we would choose the lowest power that reached an R2 value of 0.8 or higher. However, no soft-thresholding power was even close to reaching 0.8. According to the [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html), this indicates the following:

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
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30gene_dendrogram-1.png)<!-- -->


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
##  ..cutHeight not given, setting it to 0.979  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
# Look at table of modules
table(dynamicMods)
```

```
## dynamicMods
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
## 2454 1889 1464 1197 1053  458  365  293  267  214  209  186  163  126  123  109 
##   17   18 
##  107   42
```

```r
# Convert module numbers into colors
dynamicColors <- labels2colors(dynamicMods)
```


```r
# Plot dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30gene_dendrogram_module_colors-1.png)<!-- -->

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
# Plot dendrogram of clustered eigengenes
plot(METree, main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")
# ID cut height based on sample number (3)
dynamicMergeCut(9)
```

```
## [1] 0.5278481
```

```r
MEDissThres <- dynamicMergeCut(3)
```

```
## Warning in function dynamicMergeCut: too few observations for the dynamic assignment of the merge threshold.
##     Will set the threshold to .35
```

```r
abline(h = MEDissThres, col = "red")
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30clustered_eigengenes-1.png)<!-- -->



```r
merge <- mergeCloseModules(CrabExpr0, dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
```

```
##  mergeCloseModules: Merging modules whose distance is less than 0.35
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 18 module eigengenes in given set.
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 8 module eigengenes in given set.
##    Calculating new MEs...
##    multiSetMEs: Calculating module MEs.
##      Working on set 1 ...
##      moduleEigengenes: Calculating 8 module eigengenes in given set.
```

```r
# Extract merged colors and eigengenes
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
```


```r
# Plot dendrogram with original and merged eigengenes
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
```

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30dendrogram_w_orig_and_merged_eigengenes-1.png)<!-- -->


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
##                     crab       day
## MEmagenta    0.001994249 0.7876861
## MEtan        0.911617610 0.2282640
## MEblue       0.794952306 0.3448483
## MElightcyan  0.040599761 0.3140314
## MEcyan       0.646653210 0.8642174
## MEyellow     0.757132492 0.6037052
## MEbrown      0.253235959 0.2113560
## MElightgreen 0.606465668 0.1904946
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

![](../output/WGCNA_output/AmbCrabs_cbai_transcriptome_v2.0_trial/OnlyCtsOver30heatmap-1.png)<!-- -->

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
##                               GS.day
## TRINITY_DN88421_c0_g1_i1  0.41194174
## TRINITY_DN88415_c0_g1_i2  0.50522995
## TRINITY_DN88474_c0_g1_i2 -0.28383361
## TRINITY_DN88416_c0_g2_i1 -0.20194517
## TRINITY_DN88416_c0_g2_i3 -0.20194517
## TRINITY_DN88416_c0_g2_i2  0.03149722
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
## TRINITY_DN88421_c0_g1_i1 0.2705973
## TRINITY_DN88415_c0_g1_i2 0.1653257
## TRINITY_DN88474_c0_g1_i2 0.4592039
## TRINITY_DN88416_c0_g2_i1 0.6023267
## TRINITY_DN88416_c0_g2_i3 0.6023267
## TRINITY_DN88416_c0_g2_i2 0.9358874
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
##                           MMmagenta      MMtan     MMblue MMlightcyan
## TRINITY_DN88421_c0_g1_i1 -0.5204650 -0.3384469 -0.6107172  -0.3097530
## TRINITY_DN88415_c0_g1_i2 -0.2486876 -0.4373686 -0.4277072  -0.3182513
## TRINITY_DN88474_c0_g1_i2  0.3774789  0.4798011  0.9875122   0.6331593
## TRINITY_DN88416_c0_g2_i1 -0.3235397 -0.2024963 -0.5581320  -0.2619714
## TRINITY_DN88416_c0_g2_i3 -0.3235397 -0.2024963 -0.5581320  -0.2619714
## TRINITY_DN88416_c0_g2_i2 -0.1722316 -0.3497127 -0.5808429  -0.4049623
##                               MMcyan   MMyellow    MMbrown MMlightgreen
## TRINITY_DN88421_c0_g1_i1  0.12337905  0.4499506  0.4694915    0.9651473
## TRINITY_DN88415_c0_g1_i2 -0.46195190  0.4822038  0.5032324    0.7914412
## TRINITY_DN88474_c0_g1_i2 -0.33850007 -0.5573664 -0.8523332   -0.6061232
## TRINITY_DN88416_c0_g2_i1  0.01984216  0.9894469  0.4328087    0.5143004
## TRINITY_DN88416_c0_g2_i3  0.01984216  0.9894469  0.4328087    0.5143004
## TRINITY_DN88416_c0_g2_i2  0.14957545  0.8938700  0.4366963    0.3968507
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
##                          p.MMmagenta   p.MMtan     p.MMblue p.MMlightcyan
## TRINITY_DN88421_c0_g1_i1   0.1508538 0.3729903 8.064760e-02    0.41727700
## TRINITY_DN88415_c0_g1_i2   0.5187668 0.2390912 2.508179e-01    0.40392177
## TRINITY_DN88474_c0_g1_i2   0.3165639 0.1911905 7.078760e-07    0.06718426
## TRINITY_DN88416_c0_g2_i1   0.3957114 0.6013152 1.183534e-01    0.49590328
## TRINITY_DN88416_c0_g2_i3   0.3957114 0.6013152 1.183534e-01    0.49590328
## TRINITY_DN88416_c0_g2_i2   0.6576907 0.3562419 1.009892e-01    0.27960634
##                           p.MMcyan   p.MMyellow   p.MMbrown p.MMlightgreen
## TRINITY_DN88421_c0_g1_i1 0.7518209 2.242733e-01 0.202285304   2.515108e-05
## TRINITY_DN88415_c0_g1_i2 0.2106206 1.886552e-01 0.167280120   1.104919e-02
## TRINITY_DN88474_c0_g1_i2 0.3729104 1.189678e-01 0.003515383   8.359387e-02
## TRINITY_DN88416_c0_g2_i1 0.9595911 3.934574e-07 0.244588262   1.566172e-01
## TRINITY_DN88416_c0_g2_i3 0.9595911 3.934574e-07 0.244588262   1.566172e-01
## TRINITY_DN88416_c0_g2_i2 0.7009114 1.154214e-03 0.239897510   2.902697e-01
```

### Obtain gene lists - for each module, and a master list with membership and gene significance info


```r
# Save gene names as probes
probes <- names(CrabExpr0)
# Write out the gene lists for each module of interest
for (module in modNames) {
  modGenes <- (moduleColors == module) # Select module probes
  modLLIDs <- probes[modGenes] # Get gene IDs
  fileName <- paste("GeneLi")
}
```
