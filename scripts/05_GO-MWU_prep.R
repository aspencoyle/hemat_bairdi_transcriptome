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

# Elevated Day 2 vs. Ambient Day 0+2, individual libraries only
geneIDs_pvals(input_file = "../graphs/DESeq2_output/elev2_vs_amb02_indiv_only/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/elev2_vs_amb02_indiv_only_pvals.csv")

# Elevated Day 2 vs. Ambient Day 0+2+17 + Elevated Day 0 + Lowered Day 0
geneIDs_pvals(input_file = "../graphs/DESeq2_output/amb0217_elev0_low0_vs_elev2/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/amb0217_elev0_low0_vs_elev2_pvals.csv")

# Elevated Day 0 vs. Elevated Day 2, indiv. libraries only
geneIDs_pvals(input_file = "../graphs/DESeq2_output/elev0_vs_elev2_indiv/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/elev0_vs_elev2_indiv_pvals.csv")

# Ambient Day 0 vs. Ambient Day 2, indiv. libraries only
geneIDs_pvals(input_file = "../graphs/DESeq2_output/amb0_vs_amb2_indiv/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/amb0_vs_amb2_indiv_pvals.csv")

# Ambient Day 0 vs. Ambient Day 17, indiv. libraries only
geneIDs_pvals(input_file = "../graphs/DESeq2_output/amb0_vs_amb17_indiv/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/amb0_vs_amb17_indiv_pvals.csv")

# Ambient Day 2 vs. Ambient Day 17, indiv. libraries only
geneIDs_pvals(input_file = "../graphs/DESeq2_output/amb2_vs_amb17_indiv/AllGenes_wcols.txt",
              blast_file = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
              output_file = "../scripts/06_running_GO-MWU/amb2_vs_amb17_indiv_pvals.csv")


