# data

Contact: Aidan Coyle, afcoyle@uw.edu

Roberts Lab, UW-SAFS

Last edited README: 2021-12-06

## Folders

- **blast_db**: BLAST databases from all NCBI Alveolata nucleotide sequences and all Arthropoda nucleotide sequences. Sequences downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files. _Alveolata_ sequences downloaded 2021-02-17, _Arthropoda_ sequences downloaded 2021-02-24

- **blast_tables**: Tables for our three transcriptomes annotated using BLASTx. Available [here](https://robertslab.github.io/resources/Genomic-Resources/). All downloaded last on 2021-11-22

- **indiv_libraries**: Trimmed individual libraries of infected crab. Downloaded from Gannet, available [here](https://gannet.fish.washington.edu/Atumefaciens/20200318_cbai_RNAseq_fastp_trimming/). Due to large file sizes, libraries were deleted from local machine after kallisto indices were created

- **jensen_archived_samples**: Sample data for collections made by Pam Jensen. Not, strictly speaking, part of this specific transcriptome analysis. Instead, performed exploratory analysis for future project. Leaving a description in the README for future analysis, but will move the folder to a dedicated repo for the new project. That repo is available at https://github.com/afcoyle/historical_hemat

- **ncbi_genomes**: Genomes of related species downloaded from NCBI database. Sequences downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files. [Link to _Chionoecetes opilio_ protein sequences](https://www.ncbi.nlm.nih.gov/genome/?term=txid41210%5BOrganism:exp%5D), and [link to _Amoebophrya sp. genome sequences](https://www.ncbi.nlm.nih.gov/genome/?term=txid1775427%5BOrganism:exp%5D). Genomes downloaded on 2021-03-02

- **pooled_libraries**: Trimmed pooled libraries of infected crab. Downloaded from Gannet, some from [here](https://gannet.fish.washington.edu/Atumefaciens/20200414_cbai_RNAseq_fastp_trimming/). Due to large file sizes, libraries were deleted from local machine after kallisto indices were created.

- **transcriptomes**: Transcriptomes created from libraries by Sam White and Grace Crandall prior to my involvement in this project. Available [here](https://robertslab.github.io/resources/Genomic-Resources/)

- **uniprot_taxa_seqs**: Downloads of all publicly-available nucleotide sequences within particular taxa. Downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files.

## Files:

**all_uniprot_info_inc_GOterms.tab**: The Swiss-Prot database. Manually downloaded from [here](https://www.uniprot.org/uniprot/?) on 2021-02-09

md5sum:d96a92789f77e091f5555ffd4fb952e2

**indiv_crab_summary.csv**: A description of the variables for each individual crab. Some variables were determined prior to any analysis by me (such as shell condition, treatment group, etc), while others were determined as a result of the analysis (plasticity, number of immune transcripts expressed, etc).

**sample_ids.csv**: Shows the statistics for each sample - its matching crab, day of sample, number of reads per transcriptome it matches to, and more. 