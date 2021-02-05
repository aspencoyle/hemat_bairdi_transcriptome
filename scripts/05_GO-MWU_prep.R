############################################
# Aidan Coyle, afcoyle@uw.edu
# Roberts Lab, UW - School of Aquatic and Fishery Science
# 2021-01-29

# This script produces the CSV needed for GO-MWU, which is:
# A CSV 2-column table of gene IDs and unadjusted p-values
#     - Header line required (contents irrelevant)

# GO-MWU also requires another input file, which we already created
# 1. A tab-delim 2-column table of gene IDs and GO terms
#     - No header
#     - Only one line per gene

# We have three relevant tables:
# 1. A table of gene IDs and GO terms (already matches Table 1)
# 2. A table of transcript IDs and p-values
# 3. A table of transcript IDs and gene IDs.

# This utilizes a function created for this script

##########################################################
library(tidyverse)

# Functions are defined in hematodinium_analysis_functions.R
source("hematodinium_analysis_functions.R")

# Ambient and low-temp treatments
geneIDs_pvals(input_file = "../graphs/amb_v_low_day02/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
              output_file = "../scripts/07_running_GO-MWU/Amb_vsLow_day02pvals.csv")

# Day 0 and day 17 treatments
geneIDs_pvals(input_file = "../graphs/day0_day17_ambient/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
              output_file = "../scripts/07_running_GO-MWU/day0_day17_ambpvals.csv")

# Elevated and ambient treatments
geneIDs_pvals(input_file = "../graphs/elev_v_amb_day02/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
              output_file = "../scripts/07_running_GO-MWU/Elev_vsAmb_day02pvals.csv")

# Elevated and low-temp treatments
geneIDs_pvals(input_file = "../graphs/elev_v_low_day02/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
              output_file = "../scripts/07_running_GO-MWU/Elev_vsLow_day02pvals.csv")


