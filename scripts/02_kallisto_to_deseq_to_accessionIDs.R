##############################
# Aidan Coyle, afcoyle@uw.edu
# 2021/01/18, Roberts Lab

# Script for 
#   - running Kallisto output files through DESeq2,
#   - turning DESeq2 output into newline-separated file of UniProt accessions
#   - creating Venn diagrams of transcript IDs and accession IDs  
##############################

#### Loading Libraries -------------------------------
library(apeglm)
library(DESeq2)
library(tidyverse)
library(VennDiagram)
library(vsn)

# Functions are defined in hematodinium_analysis_functions.R
source("hematodinium_analysis_functions.R")

#### Analyzing Trinity matrix from Kallisto output with DESeq2 -----------------------------
# Almost all analysis takes place in function deseq_analysis(). 
# Documentation in function file
# If attempting to redo a specific graph or table,
# the file DESeq_demonstration.R may be helpful

# ELEVATED DAY 2 VS. AMBIENT DAY 0+2, INDIV LIBRARIES ONLY

# Day and temperature data for libraries 
# 178, 118, 132, 359, 349, 334, 272, 294, 280 (in order)
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(temp = factor(c("Amb", "Amb", "Amb", "Amb", "Amb", "Amb",
                                         "Elev", "Elev", "Elev")),
                         day = factor(c(0, 0, 0, 2, 2, 2, 2, 2, 2)))
  
deseq_analysis(kallisto_path = "../output/kallisto_matrices/elev2_vs_amb02_indiv_only/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/elev2_vs_amb02_indiv_only",
               variable = "temp")

# ALL ELEVATED AND AMBIENT LIBRARIES, BOTH POOLED AND INDIVIDUAL

# Day and temperature data for libraries
# 178, 359, 463, 118, 349, 481, 132, 334, 485, 151, 173, 072, 
# 127, 380821, 272, 294, 280, 380825 (in order)
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(temp = factor(c("Amb", "Amb", "Amb", 
                                         "Amb", "Amb", "Amb",
                                         "Amb", "Amb", "Amb",
                                         "Amb", "Amb", "Amb",
                                         "Amb", "Amb", "Elev",
                                         "Elev", "Elev", "Elev")),
                             day = factor(c(0, 2, 17, 
                                            0, 2, 17, 
                                            0, 2, 17, 
                                            0, 0, 0, 
                                            0, 2, 2,
                                            2, 2, 2)))

deseq_analysis(kallisto_path = "../output/kallisto_matrices/amb0217_elev0_low0_vs_elev2/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/amb0217_elev0_low0_vs_elev2",
               variable = "temp")

# ELEVATED DAY 0 VS. ELEVATED DAY 2, INDIVIDUAL LIBRARIES ONLY
# 173, 072, 127, 272, 294, 280 (in order)
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(temp = factor(c("Amb", "Amb", "Amb",
                                         "Elev", "Elev", "Elev")),
                         day = factor(c(0, 0, 0, 2, 2, 2)))

deseq_analysis(kallisto_path = "../output/kallisto_matrices/elev0_vs_elev2_indiv/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/elev0_vs_elev2_indiv",
               variable = "temp")

# AMBIENT DAY 0 VS. AMBIENT DAY 2, INDIVIDUAL LIBRARIES ONLY
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(day = factor(c(0, 0, 0, 
                                        2, 2, 2)),
                         temp = factor(c("amb", "amb", "amb",
                                         "amb", "amb", "amb")))

deseq_analysis(kallisto_path = "../output/kallisto_matrices/amb0_vs_amb2_indiv/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/amb0_vs_amb2_indiv",
               variable = "day")

# AMBIENT DAY 0 VS. AMBIENT DAY 17, INDIVIDUAL LIBRARIES ONLY
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(day = factor(c(0, 0, 0, 
                                        17, 17, 17)),
                         temp = factor(c("amb", "amb", "amb",
                                         "amb", "amb", "amb")))

deseq_analysis(kallisto_path = "../output/kallisto_matrices/amb0_vs_amb17_indiv/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/amb0_vs_amb17_indiv",
               variable = "day")

# AMBIENT DAY 2 VS. AMBIENT DAY 17, INDIVIDUAL LIBRARIES ONLY
# Order should match columns of Kallisto output matrix created by Trinity (kallisto.isoform.counts.matrix)
exp_design <- data.frame(day = factor(c(2, 2, 2, 
                                        17, 17, 17)),
                         temp = factor(c("amb", "amb", "amb",
                                         "amb", "amb", "amb")))

deseq_analysis(kallisto_path = "../output/kallisto_matrices/amb2_vs_amb17_indiv/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/DESeq2_output/amb2_vs_amb17_indiv",
               variable = "day")



#### Turning transcripts into gene IDs -------------------
# Take DESeq2 output and turn it into a newline-separated 
# list of accession IDs

# Get all gene IDs for all DEGs

# Elevated Day 2 vs. Ambient Day 0+2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/elev2_vs_amb02_indiv_only/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/elev2_vs_amb02_indiv_only_DEG_IDs.txt")
# Ambient Day 0+2+17 + Elevated Day 0 + Lowered Day 0 vs. Elevated Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0217_elev0_low0_vs_elev2/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/amb0217_elev0_low0_vs_elev2_DEG_IDs.txt")
# Elevated Day 0 vs. Elevated Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/elev0_vs_elev2_indiv/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/elev0_vs_elev2_indiv_DEG_IDs.txt")
# Ambient Day 0 vs. Ambient Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0_vs_amb2_indiv/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/amb0_vs_amb2_indiv_DEG_IDs.txt")
# Ambient Day 0 vs. Ambient Day 17
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0_vs_amb17_indiv/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/amb0_vs_amb17_indiv_DEG_IDs.txt")
# Ambient Day 2 vs. Ambient Day 17
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb2_vs_amb17_indiv/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/DEG_IDs/amb2_vs_amb17_indiv_DEG_IDs.txt")



# Get all gene IDs for all genes, not just DEGs
# Ambient vs. Low Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/elev2_vs_amb02_indiv_only/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/elev2_vs_amb02_indiv_only_All_GeneIDs.txt")
# Ambient Day 0+2+17 + Elevated Day 0 + Lowered Day 0 vs. Elevated Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0217_elev0_low0_vs_elev2/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/amb0217_elev0_low0_vs_elev2_All_GeneIDs.txt")
# Elevated Day 0 vs. Elevated Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/elev0_vs_elev2_indiv/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/elev0_vs_elev2_indiv_All_GeneIDs.txt")
# Ambient Day 0 vs. Ambient Day 2
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0_vs_amb2_indiv/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb2_indiv_All_GeneIDs.txt")
# Ambient Day 0 vs. Ambient Day 17
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb0_vs_amb17_indiv/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/amb0_vs_amb17_indiv_All_GeneIDs.txt")
# Ambient Day 2 vs. Ambient Day 17
transcripts_to_geneIDs(deseq_filepath = "../graphs/DESeq2_output/amb2_vs_amb17_indiv/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt",
                       output_path = "../output/accession_n_GOids/allgenes_IDs/amb2_vs_amb17_indiv_All_GeneIDs.txt")



#### Creating Venn diagram showing overlap in DEG --------------
# 2 diagrams created in total
#       - Overlap in transcript ID
#       - Overlap in accession IDs (since not all transcript IDs matched to accession IDs)

# PART 1: OVERLAP IN TRANSCRIPT IDs

# Links to DESeq2 outputs containing transcript IDs and data with headers
# If no headers, change second line of import_DEGs() definition
filepath1 <- "../graphs/amb_v_low_day02/DEGlist_wcols.txt"
filepath2 <- "../graphs/day0_day17_ambient/DEGlist_wcols.txt"
filepath3 <- "../graphs/elev_v_amb_day02/DEGlist_wcols.txt"
filepath4 <- "../graphs/elev_v_low_day02/DEGlist_wcols.txt"

# Use our created function to import vector of DEG transcript IDs
# Function present in hematodinium_analysis_functions.R
trans_one <- import_DEGs(filepath1)
trans_two <- import_DEGs(filepath2)
trans_three <- import_DEGs(filepath3)
trans_four <- import_DEGs(filepath4)

# Check how much our time-comparison genes overlap with 
# temperature-comparison genes
sum(trans_two %in% trans_one)
sum(trans_two %in% trans_three)
sum(trans_two %in% trans_four)
# Results are 2, 0, 2 - minimal overlap.
# As a result, only including temp comparisons in Venn diagram

# Create list of all DEGs from temperature comparisons
temp_DEGs <- list(trans_one, trans_three, trans_four)
# Create Venn diagram
venn.diagram(x = temp_DEGs,
             
             # Title and labels
             main = "Day 0/2 Individual DEGs - Transcript ID Overlap",
             main.cex = 1.25,
             cat.pos = c(0, 0, 180),
             category.names = c("Ambient vs. Low",
                                "Elevated vs. Ambient",
                                "Elevated vs. Low"),
             
             # Circles
             lwd = 2, 
             lty = "blank",
             fill = c("#FF949B", "#358AB8", "#FFFDAD"),
             euler.d = TRUE,
             scaled = TRUE,
             
             # Numbers
             cex = 1.5,
             fontface = "bold", 
             fontfamily = "sans",
             
             # Output features
             filename = "TranscriptID_DEGs.png",
             height = 3000, 
             width = 3000,
             resolution = 500,
             output = TRUE)

# How many unique transcript IDs do we have for all our DEGs?
all_IDs <- c(trans_one, trans_two, trans_three, trans_four)
length(unique(all_IDs))
# How many total transcript IDs do we have for our DEGs?
length(all_IDs)
# Percentage that are unique
length(unique(all_IDs)) / length(all_IDs)


# PART 2: OVERLAP IN DEG ACCESSION IDS

# Links to newline-separated file of accession IDs
filepath1 <- "../output/accession_n_GOids/DEG_IDs/Amb_vsLow_DEG_IDs.txt"
filepath2 <- "../output/accession_n_GOids/DEG_IDs/day0_day17_amb_DEG_IDs.txt"
filepath3 <- "../output/accession_n_GOids/DEG_IDs/Elev_vsAmb_DEG_IDs.txt"
filepath4 <- "../output/accession_n_GOids/DEG_IDs/Elev_vsLow_DEG_IDs.txt"

# Use our created function to import vector of DEG Accession IDs
# Using unlist to make it a vector, since we're reading in directly
trans_one <- unlist(read.table(file = filepath1, header = FALSE, sep = "\t"))
trans_two <- unlist(read.table(file = filepath2, header = FALSE, sep = "\t"))
trans_three <- unlist(read.table(file = filepath3, header = FALSE, sep = "\t"))
trans_four <- unlist(read.table(file = filepath4, header = FALSE, sep = "\t"))

# Check how much our time-comparison genes overlap with 
# temperature-comparison genes
sum(trans_two %in% trans_one)
sum(trans_two %in% trans_three)
sum(trans_two %in% trans_four)
# Results are 2, 1, 3. Still minimal, 
# so only including temp comparisons in Venn diagram

# Create list of all DEGs from temperature comparisons
temp_DEGs <- list(trans_one, trans_three, trans_four)
# Create Venn diagram
venn.diagram(x = temp_DEGs,
             
             # Title and labels
             main = "Day 0/2 Individual DEGs - Accession ID Overlap",
             main.cex = 1.25,
             cat.pos = c(-30, 30, 180),
             cat.fontface = "bold",
             category.names = c("Ambient vs. Low",
                                "Elevated vs. Ambient",
                                "Elevated vs. Low"),
             
             # Circles
             lwd = 2, 
             lty = "blank",
             fill = c("#FF949B", "#358AB8", "#FFFDAD"),
             euler.d = TRUE,
             scaled = TRUE,
             
             # Numbers
             cex = 1.5,
             fontface = "bold", 
             fontfamily = "sans",
             
             # Output features
             filename = "AccessionIDs_DEGs.png",
             height = 3000, 
             width = 3000,
             resolution = 500,
             output = TRUE)

# How many unique Accession IDs do we have for our DEGs?
all_IDs <- c(trans_one, trans_two, trans_three, trans_four)
length(unique(all_IDs))
# How many total Accession IDs do we have for our DEGs?
length(all_IDs)
# Percentage that are unique
length(unique(all_IDs)) / length(all_IDs)
