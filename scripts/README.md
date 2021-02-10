Aidan Coyle, afcoyle@uw.edu
Updated 2021-02-10

**Outline of each script and contents:**

**01_download_libraries_run_kallisto.ipynb** 

Downloads individual and pooled libraries, downloads two transcriptomes, creates kallisto index files from transcriptomes, runs kallisto on libraries, and then combines output into matrix. Kallisto is quite computationally intensive, and as a result, this was done on a remote machine. Therefore, commands are copied and pasted, not run directly in the notebook.

**02_kallisto_to_deseq_to_accessionIDs.R**

As the name implies, takes our matrices of Kallisto counts (created in 01_download_libraries...) and analyzes using DESeq2, producing two files - one containing all transcript IDs, and another of differentially-expressed transcripts. Those files are then taken and compared to a previously-existing table of transcripts and genes, and genes are matched to transcripts, producing accession IDs. For visualization, it then creates Venn diagrams of the overlap in differentially-expressed transcript IDs and accession IDs between our treatment conditions.

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