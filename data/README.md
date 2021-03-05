# data

Contact: Aidan Coyle, afcoyle@uw.edu

Roberts Lab, UW-SAFS

Last edited README: 2021-03-04

## Folders

- **blast_db**: BLAST databases from all NCBI Alveolata nucleotide sequences and all Arthropoda nucleotide sequences. Sequences downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files. _Alveolata_ sequences downloaded 2021-02-17, _Arthropoda_ sequences downloaded 2021-02-24

- **indiv_libraries**: Trimmed individual libraries of infected crab. Downloaded from Gannet, available [here](https://gannet.fish.washington.edu/Atumefaciens/20200318_cbai_RNAseq_fastp_trimming/). Due to large file sizes, libraries were deleted from local machine after kallisto indices were created

- **ncbi_genomes**: Genomes of related species downloaded from NCBI database. Sequences downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files. [Link to _Chionoecetes opilio_ protein sequences](https://www.ncbi.nlm.nih.gov/genome/?term=txid41210%5BOrganism:exp%5D), and [link to _Amoebophrya sp. genome sequences](https://www.ncbi.nlm.nih.gov/genome/?term=txid1775427%5BOrganism:exp%5D). Genomes downloaded on 2021-03-02

- **pooled_libraries**: Trimmed pooled libraries of infected crab. Downloaded from Gannet, some from [here](https://gannet.fish.washington.edu/Atumefaciens/20200414_cbai_RNAseq_fastp_trimming/). Due to large file sizes, libraries were deleted from local machine after kallisto indices were created.

- **transcriptomes**: Transcriptomes created from libraries by Sam White and Grace Crandall prior to my involvement in this project. Available [here](https://robertslab.github.io/resources/Genomic-Resources/)

- **uniprot_taxa_seqs**: Downloads of all publicly-available nucleotide sequences within particular taxa. Downloaded from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cg) as FASTA files.

## Files:

**all_uniprot_info_inc_GOterms.tab**: The Swiss-Prot database. Manually downloaded from [here](https://www.uniprot.org/uniprot/?) on 2021-02-09

md5sum:d96a92789f77e091f5555ffd4fb952e2

**cbai_hemat_diamond_blastx_table_transcriptome_v2.0.txt**: BLASTx annotation of cbai_transcriptome_v2.0.fasta (by my notation, cbai_hemat_transcriptome_v2.0.fasta). Annotation description available [here](https://robertslab.github.io/sams-notebook/2020/05/02/Transcriptome-Assembly-C.bairdi-All-RNAseq-Data-Without-Taxonomic-Filters-with-Trinity-on-Mox.html), [direct file available here](https://gannet.fish.washington.edu/Atumefaciens/20200508_cbai_diamond_blastx_transcriptome-v2.0/20200507.C_bairdi.Trinity.blastx.outfmt6), and [original transcriptome available here](https://owl.fish.washington.edu/halfshell/genomic-databank/cbai_transcriptome_v2.0.fasta)

md5sum: ace82a75cb947574ac807d868427253c