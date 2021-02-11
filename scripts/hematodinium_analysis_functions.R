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
  
  # Create DESeq object that looks at effect of day
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
  # Create plot
  PCA_plot <- plotPCA(vsd, intgroup = variable)
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
  
  # Separate Gene ID to specifically get Uniprot accession ID
  transcript_key <- separate(data = transcript_key, col = Gene_ID, into = c("sp", "Accession_ID", "species"), 
                             sep = "\\|")
  
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

#### geneIDs_pvals(): Turn DESeq2 output into CSV of gene IDs and p-values --------------------
# P-values are unadjusted
# Creates one of the input files for GO-MWU

# input_file: DESeq2 output file containing transcript IDs and unadjusted p-values
# blast_file: path that leads to transcript ID/gene ID table
# output_file: path to the output file, ending in .csv

geneIDs_pvals <- function(input_file, blast_file, output_file) {
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
  transcript_key <- transcript_data[,c("Transcript_ID", "pvalue", "Gene_ID")]
  
  # Separate Gene ID to specifically get Uniprot accession ID
  transcript_key <- separate(data = transcript_key, col = Gene_ID, into = c("sp", "Accession_ID", "species"), 
                             sep = "\\|")
  
  # Remove all columns except transcript ID, accession ID, and p-value
  # and remove all rows with an NA accession ID
  transcript_key <- transcript_key[!is.na(transcript_key$Accession_ID), c(4, 2)]
  
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
