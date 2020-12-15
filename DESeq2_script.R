####################################
# Aidan Coyle, afcoyle@uw.edu
# 12/14/20
# Attempt at running mixed bairdii/hemat transcriptome
# through Deseq2, using Steven's code
# Github file found at:
# /sr320/nb-2019/C_bairdi/11-Deseq2.Rmd
##################################################

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)

data <- read.table("output/kallisto_output_transcriptome_v3.0/kallisto.isoform.counts.matrix",
                              header = TRUE, sep = "\t")

names(data)[1] <- "target_ID"

rownames(data) <- data$target_ID

data <- data[,-1]

head(data)
str(data)

deseq2.colData <- data.frame(day = factor(c(2, 2, 0, 2, 17)),
                   temp = factor(c("Low", "High", "Amb", "All", "All")))

rownames(deseq2.colData) <- colnames(data)

deseq2.dds <- DESeqDataSetFromMatrix(countData = (data),
                                     colData = deseq2.colData,
                                     design = ~ day)
