##############################
# Aidan Coyle, afcoyle@uw.edu
# Roberts Lab, UW-SAFS
# 2021-02-12

# This script takes a file of differentially-expressed
# genes, and extracts the TPM (transcripts per million)
# counts from the kallisto indices created earlier
#
# Counts are only extracted from individual libraries
##############################


library(tidyverse)


# Import DEG list for Amb 2 vs. Elev 2
DEG_list <- read.table("../graphs/DESeq2_output/amb2_vs_elev2_indiv/DEGlist_wcols.txt",
                       header = TRUE,
                       sep = "\t")

# Specify file you want to write results to
outpath <- "../output/TPM_counts/amb2_vs_elev2_DEG_TPMs.txt"

# Transcript IDs are rownames - move them into first column
DEG_list <- rownames_to_column(DEG_list,
                               "Transcript_ID")

# Remove all columns that aren't transcript ID
DEG_list <- DEG_list[,1, drop = FALSE]

# Read BLAST data into R
blast_data <- read.table("../data/cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt", 
                         header = FALSE,
                         sep = "\t")

# Columns have no names - add names for first two columns
colnames(blast_data)[1:2] <- c("Transcript_ID", "Gene_ID")

# Turn the first two columns of BLAST data into a Transcript ID/Gene ID key
blast_data <- blast_data %>%
  select(Transcript_ID, Gene_ID)

# Add Gene ID column to transcript data, using Transcript ID column to match
DEG_list <- left_join(DEG_list, blast_data, by = "Transcript_ID")

# Select only DEGs with transcript IDs that match to genes
DEG_list <- DEG_list[!is.na(DEG_list$Gene_ID),]

# List all kallisto indices for individual libraries only. 
# To get indices from all libraries, change id??? to id*
kallisto_files <- Sys.glob("../output/kallisto_libraries_bairdihemat_transcriptomev2.0/id???/abundance.tsv")



for (i in 1:length(kallisto_files)) {
  # Extract the ID number from the kallisto file
  idnum <- str_extract(kallisto_files[i], "id...")
  # Read in the kallisto file
  kallisto_output <- read.delim(file = kallisto_files[i], 
                                header = TRUE,
                                sep = "\t")
  # Select only transcript ID and TPM (transcripts per million) columns
  kallisto_output <- kallisto_output %>%
    select(target_id, tpm)
  
  # Rename kallisto column names
  colnames(kallisto_output)[1:2] <- c("Transcript_ID", 
                                      paste0(idnum, "_TPM"))
  # Add TPM value to table of DEGs
  DEG_list <- left_join(DEG_list, kallisto_output, by = "Transcript_ID")
}

# Write results to table
write.table(DEG_list,
            file = outpath,
            quote = FALSE,
            row.names = FALSE)
