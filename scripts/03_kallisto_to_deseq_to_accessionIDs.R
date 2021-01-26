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

# DAY 0 VS DAY 17, AMBIENT TEMPERATURE, INDIVIDUAL LIBRARIES

# Day and temperature data for libraries 
# 118, 132, 178, 463, 481, 485 (in order)
# Order should match columns of Kallisto output matrix created by Trinity
exp_design <- data.frame(day = factor(c(0, 0, 0, 17, 17, 17)),
           temp = factor(c("Amb", "Amb", "Amb", "Amb", "Amb", "Amb")))

deseq_analysis(kallisto_path = "../output/kallisto_indivlibs_transcriptome_v3.0/day0_day17_matrix/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/day0_day17_ambient",
               variable = as.formula("day"))

# AMBIENT VS. LOW-TEMPERATURE, DAY 0-2, INDIVIDUAL LIBRARIES

# Day and temperature data for libraries
# 118, 132, 178, 334, 349, 359, 151, 254
exp_design <- data.frame(temp = factor(c("Amb", "Amb", "Amb",
                                             "Amb", "Amb", "Amb",
                                             "Low", "Low")),
                             day = factor(c(0, 0, 0, 2, 2, 2, 
                                            0, 2)))

deseq_analysis(kallisto_path = "../output/kallisto_indivlibs_transcriptome_v3.0/day02_amb_vs_low_matrix/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/amb_v_low_day02",
               variable = as.formula("temp"))

# ELEVATED VS. LOW-TEMPERATURE, DAY 0-2, INDIVIDUAL LIBRARIES

# Day and temperature data for libraries
# 127, 173, 72, 272, 280, 294, 151, 254
exp_design <- data.frame(temp = factor(c("Elev", "Elev", "Elev",
                                             "Elev", "Elev", "Elev",
                                             "Low", "Low")),
                             day = factor(c(0, 0, 0, 2, 2, 2, 0, 2)))

deseq_analysis(kallisto_path = "../output/kallisto_indivlibs_transcriptome_v3.0/day02_elev_vs_low_matrix/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/elev_v_low_day02",
               variable = as.formula("temp"))

# ELEVATED VS. AMBIENT TEMPERATURE, DAY 0-2, INDIVIDUAL LIBRARIES

# Day and temperature data for libraries
# 127, 173, 72, 272, 280, 294, 118, 132, 178, 334, 349, 359
exp_design <- data.frame(temp = factor(c("Elev", "Elev", "Elev",
                                             "Elev", "Elev", "Elev",
                                             "Amb", "Amb", "Amb", "Amb", 
                                             "Amb", "Amb")),
                             day = factor(c(0, 0, 0, 2, 2, 2, 
                                            0, 0, 0, 2, 2, 2)))

deseq_analysis(kallisto_path = "../output/kallisto_indivlibs_transcriptome_v3.0/day02_elev_vs_amb_matrix/kallisto.isoform.counts.matrix",
               experiment_table = exp_design,
               output_path = "../graphs/elev_v_amb_day02",
               variable = as.formula("temp"))

#### Turning transcripts into gene IDs -------------------
# Take DESeq2 output and turn it into a newline-separated 
# list of accession IDs

# Get all gene IDs for all DEGs

# Ambient vs. Low Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/amb_v_low_day02/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Amb_vsLow_DEG_IDs.txt")
# Day 0 vs. Day 17, Ambient Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/day0_day17_ambient/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/day0_day17_amb_DEG_IDs.txt")
# Elevated vs. Ambient Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/elev_v_amb_day02/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Elev_vsAmb_DEG_IDs.txt")
# Elevated vs. Low Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/elev_v_low_day02/DEGlist_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Elev_vsLow_DEG_IDs.txt")

# Get all gene IDs for all genes, not just DEGs
# Ambient vs. Low Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/amb_v_low_day02/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Amb_vsLow_All_GeneIDs.txt")
# Day 0 vs. Day 17, Ambient Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/day0_day17_ambient/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/day0_day17_All_GeneIDs.txt")
# Elevated vs. Ambient Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/elev_v_amb_day02/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Elev_vsAmb_All_GeneIDs.txt")
# Elevated vs. Low Temperature
transcripts_to_geneIDs(deseq_filepath = "../graphs/elev_v_low_day02/AllGenes_wcols.txt", 
                       blast_filepath = "../data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt",
                       output_path = "../output/signif_accession_ids/Elev_vsLow_All_GeneIDs.txt")

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
filepath1 <- "../output/signif_accession_ids/Amb_vsLow_DEG_IDs.txt"
filepath2 <- "../output/signif_accession_ids/day0_day17_amb_DEG_IDs.txt"
filepath3 <- "../output/signif_accession_ids/Elev_vsAmb_DEG_IDs.txt"
filepath4 <- "../output/signif_accession_ids/Elev_vsLow_DEG_IDs.txt"

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
