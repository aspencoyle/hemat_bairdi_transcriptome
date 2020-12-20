####################################
# Aidan Coyle, afcoyle@uw.edu
# 12/14/20
# Running mixed bairdii/hemat transcriptome
# through Deseq2, partially based off of Steven's code in nb-2019

# Analysis: Day 0 vs Day 17 ambient temperature

# Github file found at:
# /sr320/nb-2019/C_bairdi/11-Deseq2.Rmd
##################################################


BiocManager::install("vsn")

library(DESeq2)
library(apeglm)
library(vsn)

# Read in the matrix created by Trinity from our Kallisto outputs
data <- read.table("kallisto.isoform.counts.matrix",
                   header = TRUE, sep = "\t")

# Rename first column of matrix
names(data)[1] <- "target_ID"

# Set row names equal to the first column
rownames(data) <- data$target_ID

# Remove the now-irrelevant first column
data <- data[,-1]

# Make sure everything looks okay
head(data)
str(data)

# Round counts to integers (needed for DESeqDataSetFromMatrix()
data <- round(data, digits = 0)

# Day and temperature data for libraries 
# 118, 132, 178, 463, 481, 485 (in order)
deseq2.colData <- data.frame(day = factor(c(0, 0, 0, 17, 17, 17)),
                             temp = factor(c("Amb", "Amb", "Amb", "Amb", "Amb", "Amb")))

# Rename rows to correspond to library numbers
rownames(deseq2.colData) <- colnames(data)

# Check that it looks alright
deseq2.colData

# Create DESeq object that looks at effect of day
deseq2.dds <- DESeqDataSetFromMatrix(countData = (data),
                                  colData = deseq2.colData,
                                  design = ~ day)

deseq2.dds <- DESeq(deseq2.dds)

#Look at results
deseq2.res <- results(deseq2.dds)
deseq2.res
summary(deseq2.res)

# Shrink LFC estimates - used in some but not all analyses
resultsNames(deseq2.dds)
resLFC <- lfcShrink(deseq2.dds, coef = "day_17_vs_0", 
                    type = "apeglm")
resLFC

# Dimensions of non-NA results with an adjusted p-value of <= 0.05
dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])

# Sum of non-NA results w/ adjusted p-values less than 0.1
sum(deseq2.res$padj < 0.1, na.rm = TRUE)

# Look specifically at results w/ adjusted p-value < 0.05
deseq_res05 <- results(deseq2.dds, alpha = 0.05)
summary(deseq_res05)
head(deseq_res05)
sum(deseq_res05$padj < 0.05, na.rm = TRUE)

# Create MA plots
# Plot of full results, not shrunken
plotMA(deseq2.res, ylim = c(-28, 28))
# Plot of full results, shrunken
plotMA(resLFC, ylim = c(-2, 2), main = "apeglm")
# Plot of res05 results, not shrunken
plotMA(deseq_res05, ylim = c(-20, 20))
       
# Create a plot of Log2 fold change vs. normalized counts
deseq2_tmp <- deseq2.res
plot(deseq2_tmp$baseMean, deseq2_tmp$log2FoldChange, pch = 20,
     cex = 0.45, ylim = c(-28, 28), log = "x", col = "darkgray",
     main = "Differences by Date (padj <= 0.005)",
     xlab = "mean of normalized counts",
     ylab = "Log2 Fold Change")
# Get significant points, plot again so they're a diff color
deseq2_tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) &
                               deseq2.res$padj <= 0.005, ]
points(deseq2_tmp.sig$baseMean, deseq2_tmp.sig$log2FoldChange,
       pch = 20, cex = 0.45, col = "red")
abline(h=c(-1,1), col = "blue")

# Plot count of minimum gene
plotCounts(deseq2.dds, gene = which.min(deseq2.res$padj),
           intgroup = "day")

# Plot PCA of samples
# Transform values
vsd <- vst(deseq2.dds, blind = FALSE)
head(assay(vsd), 3)
# Create plot
plotPCA(vsd, intgroup = "day")

# Plot dispersion estimates
plotDispEsts(deseq2.dds)


# Write significant day-differing genes to table
write.table(deseq2_tmp.sig, "0vs17_DEGlist.txt",
            row.names = TRUE, col.names = FALSE, quote = FALSE,
            sep = "\t")
write.table(deseq2_tmp.sig, "0vs17_DEGlist_wcols.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE,
            sep = "\t")
