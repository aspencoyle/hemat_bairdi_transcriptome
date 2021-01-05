#########################################
# Aidan Coyle, afcoyle@uw.edu
# Roberts Lab, UW-SAFS
# 2020/12/22
# Objective: Take file of GO terms, turn into GOslim terms
#########################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("ALL")
library(topGO)
library(ALL)

# Set filepath
filepath1 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/Amb_vsLow_GO_IDs.txt"

# Read file of GO terms into R
ALL_pro_GO <- read.table(file = filepath1, sep = "\t")

# Give our file column names
colnames(ALL_pro_GO) <- c("protein_ID", "GO_IDs")

geneID2GO_bkgd <- readMappings(file = filepath1)

topgo_sig <- ALL_pro_GO
# Switch semicolons between GO IDs to commas
topgo_sig$GO_IDs <- gsub(";",",",topgo_sig$GO_IDs)

genesOfInterest <- topgo_sig

genesOfInterest <- as.character(genesOfInterest$protein_ID)
