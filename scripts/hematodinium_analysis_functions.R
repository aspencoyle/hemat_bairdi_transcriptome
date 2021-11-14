###################################
# Aidan Coyle, afcoyle@uw.edu
# Roberts Lab, UW-SAFS
# 2021-01-19
# File containing functions created during analysis of Hematodinium transcriptome
###################################


#### deseq_analysis(): Analyzing Kallisto output with DESeq2 -------------------------------------
# # A function to analyze Kallisto output with DESeq2. 
# Will likely take quite long to run - there's some 
# computationally-heavy steps involved

# Takes 4 inputs
#       kallisto_path: path to matrix created by Trinity from Kallisto outputs
#       experiment_table: table with treatments in experiment. 
#               Order should match columns of matrix from Trinity
#               Example: Day and temperature data for libraries 
#                        118, 132, 178, 463, 481, 485 (in order)
#                        experiment_table <- data.frame(day = factor(c(0, 0, 0, 17, 17, 17)),
#                             temp = factor(c("Amb", "Amb", "Amb", "Amb", "Amb", "Amb")))
#       output_path: path to the output folder
#       variable: variable being examined. Used in creation of deseq object.
#           Should be 'day' or 'temp' (without quotes)
deseq_analysis <- function(kallisto_path, 
                           experiment_table, 
                           output_path,
                           variable) {
  # Read in matrix created by Trinity from Kallisto outputs
  data <- read.table(kallisto_path, header = TRUE,
                     sep = "\t")
  # Rename first column
  names(data)[1] <- "target_ID"
  
  # Set row names equal to the first column
  rownames(data) <- data$target_ID
  
  # Remove the now-irrelevant first column
  data <- data[,-1]
  
  # Make sure everything looks okay
  print("HEAD")
  print(head(data))
  
  print("STRUCTURE")
  print(str(data))
  
  # Round counts to integers - needed for DESeqDataSetFromMatrix()
  data <- round(data, digits = 0)
  
  # Rename rows to correspond to library numbers
  rownames(experiment_table) <- colnames(data)
  
  # Check that experiment_table appears to match columns with matrix from Trinity
  print("EXPERIMENTAL DESIGN")
  print(experiment_table)
  
  # Create DESeq object that looks at effect of variable
  deseq2.dds <- DESeqDataSetFromMatrix(countData = (data),
                                       colData = experiment_table,
                                       design = as.formula(paste0("~", variable)))
 
  deseq2.dds <- DESeq(deseq2.dds)
  
  #Look at results
  deseq2.res <- results(deseq2.dds)
  print("SUMMARY:")
  print(summary(deseq2.res))
  
  # Shrink LFC estimates - used in shrunken MA plot
  lfcnames <- resultsNames(deseq2.dds)
  print(lfcnames[2])
  resLFC <- lfcShrink(deseq2.dds, coef = lfcnames[2], 
                      type = "apeglm")
  
  # Look specifically at results w/ adjusted p-value < 0.05
  deseq_res05 <- results(deseq2.dds, alpha = 0.05)
  print("Number of DEGs (padj <= 0.05)")
  print(sum(deseq_res05$padj < 0.05, na.rm = TRUE))
  
  # Plot of full results, not shrunken
  plotMA(deseq2.res, ylim = c(-28, 28))
  dev.copy(png, file.path(output_path, "allres_MAplot.png"))
  dev.off()
  
  # Plot of full results, shrunken
  plotMA(resLFC, ylim = c(-2, 2), main = "apeglm")
  dev.copy(png, file.path(output_path, "allres_shrunken_MAplot.png"))
  dev.off()
  
  # Plot of res05 results, not shrunken
  plotMA(deseq_res05, ylim = c(-20, 20))
  dev.copy(png, file.path(output_path, "res05_MAplot.png"))
  dev.off()
  
  # Create a plot of Log2 fold change vs. normalized counts
  deseq2_tmp <- deseq2.res
  plot(deseq2_tmp$baseMean, deseq2_tmp$log2FoldChange, pch = 20,
       cex = 0.45, ylim = c(-28, 28), log = "x", col = "darkgray",
       main = paste("Differences by", variable, "(padj <= 0.005)"),
       xlab = "mean of normalized counts",
       ylab = "Log2 Fold Change")
  # Get significant points, plot again so they're a diff color
  deseq2_tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) &
                                 deseq2.res$padj <= 0.005, ]
  points(deseq2_tmp.sig$baseMean, deseq2_tmp.sig$log2FoldChange,
         pch = 20, cex = 0.45, col = "red")
  abline(h=c(-1,1), col = "blue")
  dev.copy(png, file.path(output_path, "normalizedcts_v_log2foldchange.png"))
  dev.off()
  
  # Plot PCA of samples
  # Transform values
  vsd <- vst(deseq2.dds, blind = FALSE)
  head(assay(vsd), 3)

  pcaData <- DESeq2::plotPCA(vsd, intgroup=variable, returnData=TRUE)
  # Round the percent variance for the labels
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  # Create the actual plot 
  PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color = pcaData[, variable])) +
    geom_point(size=3) + 
    scale_color_discrete(variable) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  
  # Print the PCA plot and save to output path
  print(PCA_plot)
  dev.copy(png, file.path(output_path, "PCA_plot.png"))
  dev.off()
  
  # Plot dispersion estimates
  plotDispEsts(deseq2.dds)
  dev.copy(png, file.path(output_path, "dispersion_estimates.png"))
  dev.off()
  
  # write all genes to table
  write.table(deseq2_tmp, file.path(output_path, "AllGenes.txt"), 
              row.names = TRUE, col.names = FALSE, 
              quote = FALSE, sep = "\t")
  write.table(deseq2_tmp, file.path(output_path, "AllGenes_wcols.txt"),
              row.names = TRUE, col.names = TRUE,
              quote = FALSE, sep = "\t")
  
  # Write significant day-differing genes to table
  write.table(deseq2_tmp.sig, file.path(output_path, "DEGlist.txt"),
              row.names = TRUE, col.names = FALSE, quote = FALSE,
              sep = "\t")
  write.table(deseq2_tmp.sig, file.path(output_path, "DEGlist_wcols.txt"),
              row.names = TRUE, col.names = TRUE, quote = FALSE,
              sep = "\t")
  
  
}

#### transcripts_to_geneIDs(): turns DESeq2 output to UniProt accessions --------------------------------
# deseq_filepath: leads to file containing gene list from DESeq2
# blast_filepath: path that leads to transcript ID/gene ID table
# output_path: path to a new newline-separated file

transcripts_to_geneIDs <- function(deseq_filepath, 
                                   blast_filepath,
                                   output_path) {
  # Import gene list
  transcript_data <- read.table(deseq_filepath,
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
  
  # If pipes in accession ID column, separate to get accession ID
  gene_ids <- dplyr::pull(transcript_key, Gene_ID)
  if(any(grepl("|", gene_ids, fixed = TRUE))){
    transcript_key <- separate(data = transcript_key, col = Gene_ID,
                               into = c("sp", "Accession_ID", "species"), 
                               sep = "\\|")
  } else {
    transcript_key <- transcript_key %>%
      dplyr::rename(Accession_ID = Gene_ID)
        
      }
  

  
  # Create vector of non-NA accession IDs
  accession_IDs <- na.omit(transcript_key$Accession_ID)
  
  # Create file of accession IDs separated by newline
  write_lines(x = accession_IDs, file = output_path, sep = "\n")
}


#### import_DEGs(): Import file of DEGs with column headers -------------
# Used for creating Venn diagrams
import_DEGs <- function(filepath) {
  full_data <- read.table(file = filepath, header = TRUE, sep = "\t")
  full_data <- tibble::rownames_to_column(full_data, "Transcript_ID")
  transcripts <- full_data$Transcript_ID
  return(transcripts)
} 

#### geneIDs_foldchange(): Turn DESeq2 output into CSV of gene IDs and log2 fold change --------------------
# P-values are unadjusted
# Creates one of the input files for GO-MWU

# input_file: DESeq2 output file containing transcript IDs and log2 fold change
# blast_file: path that leads to transcript ID/gene ID table
# output_file: path to the output file, ending in .csv

geneIDs_foldchange <- function(input_file, blast_file, output_file) {
  # Import gene list
  transcript_data <- read.table(input_file,
                                header = TRUE, sep = "\t")
  # Transcript IDs are rownames - move them into first column
  transcript_data <- tibble::rownames_to_column(transcript_data, 
                                                "Transcript_ID")
  # Read BLAST data into R
  blast_data <- read.table(blast_file, header = FALSE,
                           sep = "\t")
  # Columns have no names - add names for first two columns
  colnames(blast_data)[1:2] <- c("Transcript_ID", "Gene_ID")
  
  # Turn the first two columns of BLAST data into a Transcript ID/Gene ID key
  blastkey <- blast_data %>%
    select(Transcript_ID, Gene_ID)
  
  # Add Gene ID column to transcript data, using Transcript ID column to match
  transcript_data <- left_join(transcript_data, blastkey, by = "Transcript_ID")
  
  # Select only the Transcript ID, p-value, and Gene ID columns
  transcript_key <- transcript_data[,c("Transcript_ID", "log2FoldChange", "Gene_ID")]
  
  
  # If pipes in accession ID column, separate to get accession ID,
  # remove all columns except transcript ID, accession ID, and log2-fold change
  # and remove all rows with an NA accession ID
  gene_ids <- dplyr::pull(transcript_key, Gene_ID)
  if(any(grepl("|", gene_ids, fixed = TRUE))){
    transcript_key <- separate(data = transcript_key, col = Gene_ID,
                               into = c("sp", "Accession_ID", "species"), 
                               sep = "\\|")
    transcript_key <- transcript_key[!is.na(transcript_key$Accession_ID), c(4, 2)]
  } else {
    # If no pipes in accession ID column, remove all rows with an NA accession ID,
    # remove the first column (transcript ID), and reorder the 2nd and 3rd (should be accession ID, then log2-fold change)
    transcript_key <- transcript_key %>%
      rename("Gene_ID" = "Accession_ID")
    transcript_key <- transcript_key[!is.na(transcript_key$Accession_ID), c(3, 2)]
    
  }
  
  rownames(transcript_key) <- NULL
  write.csv(transcript_key, output_file, row.names = FALSE,
            quote = FALSE)
}  


#### uniprot_to_GO() ------------------------------------------
# Inputs:
# - A newline-separated file of accession IDs
# - A downloaded uncompressed copy of the SwissProt database, available at https://www.uniprot.org/uniprot/ with all GO terms included

# Output:
# - A tab-separated two-column table of accession IDs and GO terms
#     - This is one of the two inputs needed for GO-MWU

# Function arguments:
#   - accession_path: path to newline-separated file of accession IDs
#   - swissprot_db: path to a recently-downloaded uncompressed copy of 
#     the SwissProt database with all GO terms included, 
#     available here: https://www.uniprot.org/uniprot/ 

uniprot_to_GO <- function(accession_path, swissprot_path, output_path) {
  # Read in file of accession IDs
  accessionIDs <- read.table(file = accession_path,
                             header = FALSE,
                             col.names = "accessionID")
  
  # Read in uniprot data table containing all GO terms
  uniprot_info <- read.delim(file = swissprot_path, 
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
  
  write.table(x = GO_terms, file = output_path, sep = "\t",
              row.names = FALSE, 
              col.names = FALSE,
              quote = FALSE)
  
}

#### WGCNA_modules_accessions_kMEs ------------------------------------------
# Inputs:
# - The transcripts within a module produced by WGCNA
# - All kMEs for all modules in your WGCNA run (should be in same directory as significant module)
# - A BLAST file for your relevant transcriptome


# Output:
# - A tab-separated two-column table of accession IDs and module kME values for all genes (kME values for genes not in module should be set to 0)
#     - This goes in the output_GOMWU directory
# - A single column of accession IDs (this will be used to create the second input for GO-MWU)
#     - This will be put in the output_accessions directory

# Function arguments:
# - module_name: color of the module you want to examine
# - transcriptome: name of the relevant transcriptome
# - compar: libraries used in WGCNA. Should be all_crabs_no_filter or amb_vs_elev_no_filter
# - blast_path: path to the BLAST file for your relevant transcriptome
# - output_GOMWU: path to the folder you want to put your files in (should be the one you run GO-MWU in)
# - output_accessions: path to the folder you want to put your single column of accession IDs in

WGCNA_modules_accessions_kMEs <- function(module_name, transcriptome, 
                                          compar, blast_path, output_GOMWU, 
                                          output_accessions) {
  # Load in module data, which is just one column of transcript IDs in module
  module.dat <- read.delim(file = paste0("../output/WGCNA_output/", transcriptome, "/", 
                                         compar, "/hemat_level_as_var/GeneList-", 
                                         module_name, ".txt"),
                           sep = "\t",
                           col.names = "Transcript_ID")
  
  # Add column of 1s to module.dat. This'll be used to signify the transcript is in the relevant module later.
  module.dat <- module.dat %>%
    add_column(member = rep(1, times = nrow(module.dat)))
  
  # Alright, module data is now ready for joining. Let's read in the kME file and prep it
  
  # Load kME data, which is a column of all transcript IDs in WGCNA, plus the kME scores for each module
  # Each module is a separate column, but we just care about the one being examined here.
  kME.dat <- read.delim(file = paste0("../output/WGCNA_output/", transcriptome, "/", 
                                      compar, "/hemat_level_as_var/ kME_table.txt"),
                        sep = "\t")
  
  # Move rownames (the transcript IDs) to columns
  kME.dat <- kME.dat %>%
    rownames_to_column(var = "Transcript_ID")
  
  # Remove all columns in our kME data table that don't pertain to the relevant module
  kME.dat <- kME.dat %>%
    select(Transcript_ID, paste0("kME", module_name))
  
  # Rename 2nd column (the kME values) to "kME_vals" so we don't need to keep messing with the module color
  colnames(kME.dat)[2] <- "kME_vals"
  
  # kME data is now ready for joining! Let's do that
  
  full.dat <- left_join(kME.dat, module.dat, by = "Transcript_ID")
  
  # To analyze a WGCNA module using GO-MWU, all genes that aren't members of the relevant module should have their kME set to 0
  full.dat[is.na(full.dat$member), ]$kME_vals <- 0
  
  # We can now remove the membership column, since it's served its purpose of the kMEs of non-members to 0
  full.dat <- full.dat %>%
    select(-member)
  
  # We now have a two-column file of transcript IDs and kMEs for the relevant module
  # However, GO-MWU needs accession IDs, not transcript IDs.
  # Therefore, we'll read in the BLAST table and use that as a key to match accession IDs and transcript IDs
  
  # Read in BLAST table
  blast_data <- read.table(blast_path, header = FALSE,
                           sep = "\t")
  
  # Add names for first two columns of BLAST data
  colnames(blast_data)[1:2] <- c("Transcript_ID", "Gene_ID")
  
  # Select only the first two columns
  blast_data <- blast_data %>%
    select(Transcript_ID, Gene_ID)
  
  # Add Gene ID column to transcript data, using Transcript ID column to match
  full.dat <- left_join(full.dat, blast_data, by = "Transcript_ID")
  
  # If pipes in accession ID column, separate to get accession ID,
  # remove all columns except transcript ID, accession ID, and kME
  # and remove all rows with an NA accession ID
  gene_ids <- dplyr::pull(full.dat, Gene_ID)
  if(any(grepl("|", gene_ids, fixed = TRUE))){
    transcript_key <- separate(data = full.dat, col = Gene_ID,
                               into = c("sp", "Accession_ID", "species"), 
                               sep = "\\|")
    transcript_key <- transcript_key[!is.na(transcript_key$Accession_ID), c(4, 2)]
  } else {
    # If no pipes in accession ID column, remove all rows with an NA accession ID,
    # remove the first column (transcript ID), and reorder the 2nd and 3rd (should be accession ID, then kME)
    transcript_key <- transcript_key %>%
      rename("Gene_ID" = "Accession_ID")
    transcript_key <- transcript_key[!is.na(transcript_key$Accession_ID), c(3, 2)]
  }
  
  end_loc <- str_split(transcriptome, pattern = "_transcriptome")
  
  output_file <- paste0(output_GOMWU,
                        paste0(end_loc[[1]][1]),
                        paste0(end_loc[[1]][2]),
                        "_", compar, "_",
                        module_name, 
                        "_module_kMEs.csv")
  
  # Write to GOMWUoutput path
  rownames(transcript_key) <- NULL
  write.csv(transcript_key, output_file, row.names = FALSE, quote = FALSE)
  
  # At this point, we have a two-column file of accession IDs and kME values
  # We can remove the kME values column, thus producing the starting point for our 
  # two-column table of accession IDs and GO terms
  
  accession_IDs <- transcript_key %>%
    select(-kME_vals)
  
  accessions_out <- paste0(output_accessions,
                              transcriptome, "/",
                              compar, "_", 
                              module_name, "_GeneIDs.txt")
  
  # Create vector of non-NA accession IDs
  accession_IDs <- na.omit(transcript_key$Accession_ID)
  
  # Create file of accession IDs separated by newline
  write_lines(x = accession_IDs, file = accessions_out, sep = "\n")
  
  
  
}
