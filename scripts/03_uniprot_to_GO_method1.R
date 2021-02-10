###################################
# Aidan Coyle, afcoyle@uw.edu
# 2021-02-02
# Roberts Lab, UW-SAFS
###################################

# Inputs:
# - A newline-separated file of accession IDs
# - A downloaded uncompressed copy of the SwissProt database, available at https://www.uniprot.org/uniprot/ with all GO terms included

# Output:
# - A tab-separated two-column table of accession IDs and GO terms
#     - This is one of the two inputs needed for GO-MWU

# Explanation:
# This is one of two methods used to obtain GO terms from 
# accession IDs. The other - 03_uniprot_to_GO_altmethod.ipynb -
# runs a script which pulls GO terms from a newline-separated file
# of accession IDs. 

# Upside of this script (03_uniprot_to_GO_method1.R): 
#   GO terms are obtained in minutes rather than days (when ran on local machine)
# Upside of alternative method (03_uniprot_to_GO_altmethod.ipynb):
#   Does not require a database to be downloaded manually
#   GO terms are always the most up-to-date

# If you do use this script,  you should finish processing in 03_uniprot_to_GO_altmethod.ipynb
# This final processing eliminates duplicate lines, which is needed for GO-MWU
# The spot to input your files in that script is marked within the .ipynb file.

library(tidyverse)

#### First Comparison:
#### Elevated Day 2 vs. Ambient Day 0+2, indiv libraries only for all ------------

# Specify filepath of output
outputpath <- "../output/accession_n_GOids/allgenes_IDs/elev2_vs_amb02_indiv_only_All_GOIDs.txt"

# Read in file of accession IDs from elev2_vs_amb02_indiv
accessionIDs <- read.table(file = "../output/accession_n_GOids/allgenes_IDs/elev2_vs_amb02_indiv_only_All_GeneIDs.txt",
                           header = FALSE, col.names = "accessionID")

# Read in uniprot data table downloaded from https://www.uniprot.org/uniprot/
# When specifying data to download, select all GO terms
uniprot_info <- read.delim(file = "../data/all_uniprot_info_inc_GOterms.tab", 
                           header = TRUE,
                           fill = TRUE,
                           sep = '\t')

# Rename first column
colnames(uniprot_info)[1] <- "accessionID"

# Left join
all_terms <- left_join(accessionIDs, uniprot_info, by = "accessionID")

# See how many unmatched terms we have
sum(is.na(all_terms$Gene.ontology.IDs))

# Select those unmatched terms, assign them to a new table
unmatched_terms <- all_terms[is.na(all_terms$Gene.ontology.IDs),]

# Remove all unmatched terms from main table
all_terms <- all_terms[!is.na(all_terms$Gene.ontology.IDs),]

# Select only accession IDs and GO IDs in new table
GO_terms <- all_terms %>%
  select(accessionID, Gene.ontology.IDs)

write.table(x = GO_terms, file = outputpath, sep = "\t",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

#### Second comparison: 
#### ambient day 0+2+17 + elevated day 0 + lowered day 0 vs. elevated day 2 -------------------------

# Specify filepath of output
outputpath <- "../output/accession_n_GOids/allgenes_IDs/amb0217_elev0_low0_vs_elev2_All_GOIDs.txt"

# Read in file of accession IDs from elev2_vs_amb02_indiv
accessionIDs <- read.table(file = "../output/accession_n_GOids/allgenes_IDs/amb0217_elev0_low0_vs_elev2_All_GeneIDs.txt",
                           header = FALSE, col.names = "accessionID")

# Read in uniprot data table downloaded from https://www.uniprot.org/uniprot/
# When specifying data to download, select all GO terms
uniprot_info <- read.delim(file = "../data/all_uniprot_info_inc_GOterms.tab", 
                           header = TRUE,
                           sep = '\t')

# Rename first column
colnames(uniprot_info)[1] <- "accessionID"

# Left join
all_terms <- left_join(accessionIDs, uniprot_info, by = "accessionID")

# See how many unmatched terms we have
sum(is.na(all_terms$Gene.ontology.IDs))

# Select those unmatched terms, assign them to a new table
unmatched_terms <- all_terms[is.na(all_terms$Gene.ontology.IDs),]

# Remove all unmatched terms from main table
all_terms <- all_terms[!is.na(all_terms$Gene.ontology.IDs),]

# Select only accession IDs and GO IDs in new table
GO_terms <- all_terms %>%
  select(accessionID, Gene.ontology.IDs)

write.table(x = GO_terms, file = outputpath, sep = "\t",
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)


