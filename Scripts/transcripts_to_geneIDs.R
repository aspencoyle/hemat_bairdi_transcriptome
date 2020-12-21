##############################################
# 2020/12/19
# Aidan Coyle, afcoyle@uw.edu
# Lab: Steven Roberts

# Turns DESeq2 output file into newline-separated file
# of UniProt accessions
# For each output file, change variables:
  # deseq_output_filepath: Path that leads to DESeq2 output file
  # output_path: Path to your new newline-separated file

##############################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(tidyverse)
library(GSEABase)
library(topGO)

# Absolute path that leads to DESeq2 output file
deseq_output_filepath <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/graphs/elev_v_low_day02/Elev_vsLow_DEGlist_wcols.txt"
# Absolute path that leads to transcript ID/gene ID table
blast_filepath <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt"

# Read output file into R
transcript_data <- read.table(deseq_output_filepath, 
                              header = TRUE, sep = "\t")
# Transcript IDs are rownames - move them into first column
transcript_data <- tibble::rownames_to_column(transcript_data, 
                                              "Transcript_ID")
# Read BLAST data into R
blast_data <- read.table(blast_filepath, header = FALSE,
                         sep = "\t")
# Columns have no names - add names for first two columns
colnames(blast_data)[1:2] <- c("Transcript_ID", "Gene_ID")

# Turn the first two columns of BLAST data into a Transcript ID/Gene ID key
blastkey <- blast_data %>%
  select(Transcript_ID, Gene_ID)

# Add Gene ID column to transcript data, using Transcript ID column to match
transcript_data <- left_join(transcript_data, blastkey, by = "Transcript_ID")

# Select only the Transcript ID and Gene ID columns
transcript_key <- transcript_data[,c("Transcript_ID", "Gene_ID")]
length(transcript_key$Transcript_ID)
sum(is.na(transcript_key$Gene_ID))

# Separate Gene ID to specifically get Uniprot accession ID
transcript_key <- separate(data = transcript_key, col = Gene_ID, into = c("sp", "Accession_ID", "species"), 
           sep = "\\|")

# Create vector of non-NA accession IDs
accession_IDs <- na.omit(transcript_key$Accession_ID)

# Path where we want to put our file of accession IDs
output_path <- "C:/Users/acoyl/Documents/GitHub/hemat_bairdii_transcriptome/output/signif_accession_ids/Elev_vsLow_DEG_IDs.txt"
# Create our file of accession IDs separated by a newline
write_lines(x = accession_IDs, file = output_path, sep = "\n")

test <- read.table(text = accession_IDs)
