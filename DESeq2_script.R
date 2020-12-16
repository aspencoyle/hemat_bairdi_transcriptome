####################################
# Aidan Coyle, afcoyle@uw.edu
# 12/14/20
# Attempt at running mixed bairdii/hemat transcriptome
# through Deseq2, using Steven's code
# Github file found at:
# /sr320/nb-2019/C_bairdi/11-Deseq2.Rmd
##################################################

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("rhdf5")

library(DESeq2)

# Read in the matrix created by Trinity from our Kallisto outputs
data <- read.table("output/kallisto_output_transcriptome_v3.0/kallisto.isoform.counts.matrix",
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


#### OPTION 1: SPLIT LIBRARIES INTO 2 GROUPS PRIOR TO MAKING DESEQ OBJECT ###################################
# Split data into libraries 2, 4, and 6 (~same day, diff temps)
# and libraries 8 and 10 (all temperatures, day 2 and 17)
temp_data <- data[1:3]
day_data <- data[4:5]

# Temperature info for libraries 2, 4, and 6 (in order)
temperature.colData <- data.frame(condition = factor(c("Low", "High", "Amb")),
                                  type = factor(rep("single-read", 3)))

# Rename rows to correspond to library numbers
rownames(temperature.colData) <- colnames(temp_data)

# Create DESeq object that looks at temperature
temp.dds <- DESeqDataSetFromMatrix(countData = (temp_data),
                                   colData = temperature.colData,
                                   design = ~ condition)
temp.dds <- DESeq(temp.dds)
# Error in checkForExperimentalReplocates(object, modelMatrix:
# The design matrix has the same number of samples and coefficients to fit,
# so estimation of dispersion is not possible. Treating samples as replicates
# was deprecated in v1.20 and no longer supported since v1.22

#### OPTION 2: CREATE 2 DESEQ OBJECTS, BOTH FROM ALL DATA ##############################
# ONE OBJECT EXAMINES EFFECT OF TEMPERATURE, THE OTHER OF DAY

# Day and temperature data for libraries 2, 4, 6, 8, 10 (in order)
deseq2.colData <- data.frame(day = factor(c(2, 2, 0, 2, 17)),
                   temp = factor(c("Low", "High", "Amb", "All", "All")))

# Rename rows to correspond to library numbers
rownames(deseq2.colData) <- colnames(data)

# Check that it looks alright
deseq2.colData

# Create 2 DESeq objects that look at day and temperature
day.dds <- DESeqDataSetFromMatrix(countData = (data),
                                   colData = deseq2.colData,
                                   design = ~ day)
temp.dds <- DESeqDataSetFromMatrix(countData = (data),
                                     colData = deseq2.colData,
                                     design = ~ temp)

day.dds <- DESeq(day.dds)
temp.dds <- DESeq(temp.dds)

# Look at results
day.res <- results(day.dds)
day.res
temp.res <- results(temp.dds)
temp.res

# Dimensions of non-NA results with an adjusted p-value of <= 0.05
dim(day.res[!is.na(day.res$padj) & day.res$padj <= 0.05, ])
dim(temp.res[!is.na(temp.res$padj) & temp.res$padj <= 0.05, ])

# Sum of non-NA results w/ adjusted p-values less than 0.1
sum(day.res$padj < 0.1, na.rm = TRUE)
sum(temp.res$padj < 0.1, na.rm = TRUE)

# Look specifically at results w/ adjusted p-value < 0.05
day_res05 <- results(day.dds, alpha = 0.05)
summary(day_res05)
head(day_res05)

temp_res05 <- results(temp.dds, alpha = 0.05)
summary(temp_res05)
head(temp_res05)
