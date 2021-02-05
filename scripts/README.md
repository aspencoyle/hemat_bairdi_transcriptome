Outline of each script and contents:

01_download_libraries_run_kallisto.ipynb: 
Downloads individual and pooled libraries, downloads two transcriptomes, creates kallisto index files from transcriptomes, runs kallisto on libraries, and then combines output into matrix. Kallisto is quite computationally intensive, and as a result, this was done on a remote machine. Therefore, commands are copied and pasted, not run directly in the notebook.

02_kallisto_to_deseq_to_accessionIDs.R:
As the name implies, takes our matrices of Kallisto counts (created in 01_download_libraries...) and analyzes using DESeq2, producing two files - one containing all transcript IDs, and another of differentially-expressed transcripts. Those files are then taken and compared to a previously-existing table of transcripts and genes, and genes are matched to transcripts, producing accession IDs. For visualization, it then creates Venn diagrams of the overlap in differentially-expressed transcript IDs and accession IDs between our treatment conditions.

03_uniprot_to_GO.ipynb
Uses shell script created by Sam White, which is saved as 05_uniprot2go.sh. Takes newline-separated Uniprot accession IDs, created in the previous script, and obtains GO terms

04_uniprot2go.sh
Shell script created by Sam White to take newline-separated Uniprot accession IDs and obtain GO terms.

05_GO-MWU_prep
GO-MWU needs two input files - a tab-delimited 2-column table of gene IDs and GO terms with no header and one line per gene, and a 2-column CSV of gene IDs and GO terms with a header. We already created the former when we ran 04_uniprot2go.sh in our 03_uniprot_to_go.ipynb script (although we do have some genes on multiple lines, which we'll clean up later). This script creates the latter file