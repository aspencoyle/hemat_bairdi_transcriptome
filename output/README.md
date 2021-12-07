# output

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

Last edited README: 2021-12-06

## accession_n_GOids

Contains a number of files that are the result of pairwise comparisons. The comparison is indicated in the filename. Files ending in GeneIDs.txt contain a single column of newline-separated accession IDs. Files ending in GOIDs.txt contain two columns - one of accession IDs, and the other of the GO terms that match those accession IDs.

## all_genes

Contains two files for each transcriptome. One (all_indiv_libraries.txt) contains the transcript and accession IDs found in each transcriptome, excluding transcripts without a matching accession ID. The other (all_gene_names.csv) contains data from UniProtKB/SwissProt on the characteristics of each accession ID found in the transcriptome.

## BLASTs

Contains the result of many BLASTs

- **alveolata_publicseqs/**: BLASTed against a database consisting of all available Alveolata nucleotide sequences that could be found on the NCBI taxonomy browser

- **amoebo_genome/**: BLASTed against a database consisting of the genome for _Amoebophrya sp._

- **arthropoda_publicseqs/**: BLASTed against a database consisting of all available Arthropoda nucleotide sequences that could be found on the NCBI taxonomy browser

- **input_seqs/**: In cases where a FASTA file needed to be created (for instance, when specifically BLASTing DEGs only), the input to a BLAST was placed here. 

- **opilio_genome/**: BLASTed against a database consisting of the genome for _Chionoecetes opilio_

If needed, further info on filenames is often located in a README inside that directory

## correlation

Shows the correlation between sample variables, both in the form of correlation charts and dot plots. Includes both physical characteristics (e.g. shell condition and carapace width), and expression characteristics (e.g. number of immune genes expressed and percent aligned)

## GO-MWU_output

Output files from running [GO-MWU](https://github.com/z0on/GO_MWU)

## immune_genes

For each transcriptome, has the following:
- Names of all genes, regardless of whether immune genes
- Names of all immune genes (in both .csv and .txt format)
- Transcript and accession IDs for all immune genes

## input_for_GO-MWU

GO-MWU requires two input files - one containing a CSV of accession IDs and p-values, and the other containing a tab-separated file of accession IDs and GO IDs with no repeats. This folder contains both files for each parwise comparison we ran through GO-MWU

## kallisto_libraries

Contains kallisto quantifications for each library 

## kallisto_matrices

Matrices containing kallisto quantifications for each pairwise comparison examined

## manual_clustering

A clustering script was run directly in R, and its results are shown here

## TPM_counts

For sequences that were differentially-expressed in a pairwise comparison, contains transcript per million (TPM) counts from all libraries.

## WGCNA_output

Output from the R package WGCNA, which performs weighted gene correlation network analyses