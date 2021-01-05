#####################################################
# Aidan Coyle, afcoyle@uw.edu
# Lab: Steven Roberts, UW-SAFS
# 2020/12/21
# Project: Differentially-expressed genes in hematodinium/bairdii
# at different temperatures/times
# Goal: Create two Venn diagrams showing DEG overlap
  # 1: Overlap in transcript IDs
  # 2: Overlap in accession IDs (since not all transcript IDs matched to accession IDs)
#####################################################


library(VennDiagram)
library(tidyverse)

BiocManager::install("VennDiagram")


#### FUNCTION DEFINITIONS ###########################
import_DEGs <- function(filepath) {
  full_data <- read.table(file = filepath, header = TRUE, sep = "\t")
  full_data <- tibble::rownames_to_column(full_data, "Transcript_ID")
  transcripts <- full_data$Transcript_ID
  return(transcripts)
} 

#### PART 1: OVERLAP IN TRANSCRIPT IDs #############################
# Links to DESeq2 outputs containing transcript IDs and data with headers
# If no headers, change second line of import_DEGs() definition
filepath1 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/graphs/amb_v_low_day02/Amb_vsLow_DEGlist_wcols.txt"
filepath2 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/graphs/day0_day17_ambient/0vs17_DEGlist_wcols.txt"
filepath3 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/graphs/elev_v_amb_day02/Elev_vsAmb_DEGlist_wcols.txt"
filepath4 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/graphs/elev_v_low_day02/Elev_vsLow_DEGlist_wcols.txt"

# Use our created function to import vector of DEG transcript IDs
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


#### PART 2: OVERLAP IN ACCESSION IDS ###############################
# Links to newline-separated file of accession IDs
filepath1 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/Amb_vsLow_DEG_IDs.txt"
filepath2 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/day0_day17_amb_DEG_IDs.txt"
filepath3 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/Elev_vsAmb_DEG_IDs.txt"
filepath4 <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/Elev_vsLow_DEG_IDs.txt"

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
