###################################
# Aidan Coyle, afcoyle@uw.edu
# 2021-02-02
# Roberts Lab, UW-SAFS
###################################

# Explanation:
# This is one of two methods used to obtain GO terms from 
# accession IDs. The other - 03_uniprot_to_GO_altmethod.ipynb -
# runs a script which pulls GO terms from a newline-separated file
# of accession IDs. 

# You must have a recently-downloaded, uncompressed copy of the
# Swiss-Prot database with all GO terms included,
# available at https://www.uniprot.org/uniprot/ 

# Upside of this script (03_uniprot_to_GO_method1.R): 
#   GO terms are obtained in minutes rather than days (when ran on local machine)
# Upside of alternative method (03_uniprot_to_GO_altmethod.ipynb):
#   Does not require a database to be downloaded manually
#   GO terms are always the most up-to-date

# If you do use this script,  you should finish processing in 03_uniprot_to_GO_altmethod.ipynb
# This final processing eliminates duplicate lines, which is needed for GO-MWU
# The spot to input your files in that script is marked within the .ipynb file.

library(tidyverse)
source("hematodinium_analysis_functions.R")

# Elevated Day 2 vs. Ambient Day 0+2. Individual libraries 
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/elev0_vs_elev2_indiv_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab", 
              output_path = "../output/accession_n_GOids/allgenes_IDs/elev0_vs_elev2_indiv_All_GOIDs.txt")

# Ambient Day 0+2+17 + Elevated Day 0 + Lowered Day 2 vs. Elevated Day 2
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/amb0217_elev0_low0_vs_elev2_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab", 
              output_path = "../output/accession_n_GOids/allgenes_IDs/amb0217_elev0_low0_vs_elev2_All_GOIDs.txt")


# Elevated Day 0 vs. Elevated Day 2. Individual libraries 
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/elev0_vs_elev2_indiv_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab", 
              output_path = "../output/accession_n_GOids/allgenes_IDs/elev0_vs_elev2_indiv_All_GOIDs.txt")


# Ambient Day 0 vs. Ambient Day 2. Individual libraries
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb2_indiv_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab",
              output_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb2_indiv_All_GOIDs.txt")

# Ambient Day 0 vs. Ambient Day 17. Individual libraries
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb17_indiv_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab",
              output_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb17_indiv_All_GOIDs.txt")

# Ambient Day 2 vs. Ambient Day 17. Individual libraries
uniprot_to_GO(accession_path = "../output/accession_n_GOids/allgenes_IDs/amb2_vs_amb17_indiv_All_GeneIDs.txt",
              swissprot_path = "../data/all_uniprot_info_inc_GOterms.tab",
              output_path = "../output/accession_n_GOids/allgenes_IDs/amb2_vs_amb17_indiv_All_GOIDs.txt")
