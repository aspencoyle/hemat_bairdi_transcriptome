Outline of each script and contents:

01_cbai_hemat_pooled_libraries.ipynb: 
Jupyter notebook containing my attempts to analyze pooled libraries. Takes pooled libraries, runs them through kallisto, and combines output into a matrix. This workflow was discarded, as pooled libraries couldn't sample widely enough to provide a meaningful DESeq2 analysis

02_cbai_hemat_indiv_daycomparisons:
Takes individual libraries, runs them through kallisto, and combines output into a matrix. Due to issues with kallisto being able to write to folders, this was largely done via copying commands from bash, but the workflow is the same.

03_kallisto_to_deseq_to_accessionIDs:
As the name implies, takes our matrices of Kallisto counts (created in 02_cbai-hemat...) and analyzes using DESeq2, producing two files - one containing all transcript IDs, and another of differentially-expressed transcripts. Those files are then taken and compared to a previously-existing table of transcripts and genes, and genes are matched to transcripts, producing accession IDs. For visualization, it then creates Venn diagrams of the overlap in differentially-expressed transcript IDs and accession IDs between our three treatment conditions.

04_uniprot_to_GO
Uses shell script created by Sam White, which is saved as 05_uniprot2go.sh. Takes newline-separated Uniprot accession IDs, created in the previous script, and obtains GO terms

05_uniprot2go.sh
Shell script created by Sam White to take newline-separated Uniprot accession IDs and obtain GO terms.