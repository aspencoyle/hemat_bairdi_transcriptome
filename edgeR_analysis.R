#####################################
# Aidan Coyle, acoyle907@gmail.com
# 2020/12/10
# Differential analysis of hematodinium/C. bairdii libraries
# Looking at Libraries 2, 4, 6, 8, 10
# Inputting data from Kallisto, running through edgeR
# ###################################

BiocManager::install("DESeq2")
#### Loading Packages ---------------------------------
library(BiocManager)
library(edgeR)
library(rhdf5)
library(tximport)
library(tidyverse)
library(csaw)

#### Defining Functions -------------------------------

# kallisto_to_DGE(): read in Kallisto abundance.h5 output file,
# return DGEList of counts
# Three inputs
  # Name of folder containing abundance.h5 output file
  # Absolute path to folder containing kallisto output folders
  # Absolute path to blast results with transcript ID and gene ID
kallisto_to_txi <- function(folder, output_path, blast_path) {
  # Import BLAST results
  gene_matches_raw <- read_tsv(file = blast_path,
                               col_names = FALSE)
  # Select only first two columns of BLAST results - transcript and gene IDs
  gene_matches <- gene_matches_raw[,1:2]
  
  # Define full absolute path to abundance.h5 output file
  kallistopath <- paste(output_path, folder, "/abundance.h5", sep = "")
  # Import abundance.h5 file using tximport, creating a list of abundance/counts/length/countsfromabundance
  # txOut = TRUE gives transcript-level summarization.
  # Get gene-level summarization by setting txOut = FALSE
  imported <- tximport(files = kallistopath, type = "kallisto",
                       txOut = TRUE,
                       tx2gene = gene_matches)
  return(imported)
}

#### Get data properly formatted ---------------------

# Define absolute path to the folder containing your Kallisto outputs
output_location <- "C:/Users/acoyl/Documents/GradSchool/RobertsLab/crab_vs_hemat_diff_expression/output/kallisto_output_transcriptome_3.0/"
# Define absolute path to BLAST results
blast_location <- "C:/Users/acoyl/Documents/GradSchool/RobertsLab/crab_vs_hemat_diff_expression/data/cbai_hemat_diamond_blastx_table_transcriptome_v3.0.txt"

# Repeatedly run previously-defined extract_kallisto()
# Argument 1: folder name of kallisto results to be read in
# Argument 2: Absolute path to folder with Kallisto outputs (defined above)
# Argument 3: Absolute path to BLAST results (defined above)

library_02.txi <- kallisto_to_txi(folder = "library02", output_path = output_location, blast_path = blast_location)
library_04.txi <- kallisto_to_txi(folder = "library04", output_path = output_location, blast_path = blast_location)
library_06.txi <- kallisto_to_txi(folder = "library06", output_path = output_location, blast_path = blast_location)
library_08.txi <- kallisto_to_txi(folder = "library08", output_path = output_location, blast_path = blast_location)
library_10.txi <- kallisto_to_txi(folder = "library10", output_path = output_location, blast_path = blast_location)

class(library_02.txi$counts)


folder <- "library_02"

cts <- library_02.txi$counts
normMat <- library_02.txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts / normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
keep <- filterByExpr(y)
y <- y[keep, ]
y
class(y)



# Kallisto vignette
output_folder <- paste(output_location, 
                       "library02/abundance.h5", sep = "")
output_folder
txi <- tximport(files = output_folder, type = "kallisto",
                tx2gene = gene_matches)
names(txi)
head(txi$counts)
