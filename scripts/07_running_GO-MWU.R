################################
# 07_Running_GO-MWU
# Aidan Coyle, afcoyle@uw.edu
# 2021-01-31
# Roberts lab, UW-SAFS
################################

# Running GO-MWU requires 2 tables
# 1. A 2-column CSV table of gene IDs and continuous measures assc.
#     with GO enrichment (ex: unadjusted p-val)
# 2. A 2-column tab-delimited table of gene IDs and GO terms 
#     separated by ;

# Install packages
require(dichromat)

# Import 2-col table of accession IDs and p-values

geneBackground <- read.csv("../output/input_for_GO-MWU/Amb_vsLow_day02.csv")
# Confirm import
head(geneBackground)

# Import 2-col table of gene IDs and GO terms
uniprotGOTerms <- read.delim("../output/input_for_GO-MWU/Amb_vsLow_All_GOIDs.txt",
                             header = FALSE)
# Rename columns
colnames(uniprotGOTerms) <- c("Uniprot", "GO")
#Confirm import
head(uniprotGOTerms)


# 