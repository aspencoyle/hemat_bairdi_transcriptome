Aidan Coyle, afcoyle@uw.edu
Updated 2021-02-10

**Outline of each script and contents:**

**01_download_libraries_run_kallisto.ipynb** 

Downloads individual and pooled libraries, downloads two transcriptomes of mixed _C. bairdi_ and _Hematodinium_, creates kallisto index files from transcriptomes, runs kallisto on libraries, and then combines output into matrix. Kallisto is quite computationally intensive, and as a result, this was done on a remote machine. Therefore, commands are copied and pasted, not run directly in the notebook.

**02_kallisto_to_deseq_to_accessionIDs.R**

This script utilizes a series of custom-built functions for the following: 
- running a Trinity matrix of kallisto transcript counts through DESeq2
- turning DESeq2 output into newline-separated file of UniProt accessions
- creating Venn diagrams of transcript IDs and accession IDs (not needed for future steps - can skip if not wanted)

**03_uniprot_to_GO_method1.R**

This is one of two methods used to obtain one of the inputs for GO-MWU (a tab-delimited 2-column table of gene IDs and GO terms). The other is 03_uniprot_to_GO_altmethod.ipynb. This method is much faster (takes minutes), and thus preferred if you have a large number of accession IDs. However, it also requires you to download the SwissProt database. If you choose this method, you will need to use the latter portion of 03_uniprot_to_GO_altmethod.ipynb to remove duplicate rows prior to inputting to GO-MWU.

**03_uniprot_to_GO_altmethod.ipynb**

Alternative method of obtaining GO terms. Much slower, but completely automated - preferred if examining small numbers of DEGs. Uses shell script created by Sam White, which is saved as 05_uniprot2go.sh. Takes newline-separated Uniprot accession IDs, created in the previous script, and obtains GO terms. Also removes duplicate rows from your file of gene IDs and GO terms to prepare for input to GO-MWU.

**04_uniprot2go.sh**

Shell script created by Sam White to take newline-separated Uniprot accession IDs and obtain GO terms. Called in 03_uniprot_to_GO_altmethod.ipynb.

**05_GO-MWU_prep.R**

GO-MWU needs two input files - a tab-delimited 2-column table of gene IDs and GO terms with no header and one line per gene, and a 2-column CSV of gene IDs and GO terms with a header. We already created the former when we ran 04_uniprot2go.sh in our 03_uniprot_to_go.ipynb script (although we do have some genes on multiple lines, which we'll clean up later). This script creates the latter file

**06_running_GO-MWU/**

This folder contains: 

    1. an R script I created that runs GO-MWU (titled 06_running_GO-MWU)

    2. all files needed for running GO-MWU, available [here](https://github.com/z0on/GO_MWU)

    3. The specific input files for GO-MWU that we created (suboptimal structurally, but required for GO-MWU to run)

**11_obtaining_TPM_for_DEGs.Rmd**

Takes a file of genes that are differentially-expressed (output from DESeq2) and cross-references it with the transcriptome used to create the kallisto index. It then extracts the TPM (transcripts per million) counts from the kallisto libraries created earlier in the pipeline. This produces a single table containing the following:

- transcript and accession IDs for all inputted DEGs

- all TPM counts from individual kallisto libraries for those DEGs

**12_DEG_blast.ipynb**

This script takes a newline-separated file of accession IDs for all DEGs and BLASTs it (using BLASTn)against all Alveolata nucleotide sequences, and then again for all Arthropoda nucleotide sequences. The goal is to determine whether a particular DEG is more likely to be _C. bairdi_ or _Hematodinium sp._


**13_blastn_analysis.Rmd**

Takes the results from our BLASTs in 12_DEG_blast.ipynb and examines the e-values to determine whether sequences are more likely to have originated from _C. bairdi_ or _Hematodinium sp._

**21_ncbi_genome_blasts.ipynb**

Takes the closest available genome to _C. bairdi_ (_C. opilio_), and the closest available genome to _Hematodinium sp_. (_Amoebophrya sp._) and BLASTs cbai_transcriptome_v2.0 (AKA cbai_hemat_transcriptomev2.0) in a variety of ways

**Scripts 31-35**

These scripts perform the same processes as 01-05. However, while 01-05 analyze cbai_transcriptome_v2.0 - which is not filtered by taxa - these scripts analyze hemat_transcriptome_v1.6, which takes all reads and filters to only include those from Alveolata (the taxa that contains _Hematodinium_). 

**hematodinium_analysis_functions.R**

Contains a variety of functions created for the purpose of this project. Functions are called inside other R and Rmd scripts.