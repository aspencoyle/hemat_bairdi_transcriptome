# output

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

Last edited README: 2021-03-04

## accession_n_GOids

Contains a number of files that are the result of pairwise comparisons. The comparison is indicated in the filename. Files ending in GeneIDs.txt contain a single column of newline-separated accession IDs. Files ending in GOIDs.txt contain two columns - one of accession IDs, and the other of the GO terms that match those accession IDs.

## BLASTn

Contains the result of BLASTns

- **alveolata_publicseqs/**: BLASTed against a database consisting of all available Alveolata nucleotide sequences that could be found on the NCBI taxonomy browser

- **amoebo_genome/**: BLASTed against a database consisting of the genome for _Amoebophrya sp._

- **arthropoda_publicseqs/**: BLASTed against a database consisting of all available Arthropoda nucleotide sequences that could be found on the NCBI taxonomy browser

- **input_seqs/**: In cases where a FASTA file needed to be created (for instance, when specifically BLASTing DEGs only), the input to a BLAST was placed here. 

- **opilio_genome/**: BLASTed against a database consisting of the genome for _Chionoecetes opilio_

If needed, further info on filenames is often located in a README inside that directory

## GO-MWU_output

Output files from running [GO-MWU](https://github.com/z0on/GO_MWU)

## input_for_GO-MWU

GO-MWU requires two input files - one containing a CSV of accession IDs and p-values, and the other containing a tab-separated file of accession IDs and GO IDs with no repeats. This folder contains both files for each parwise comparison we ran through GO-MWU

## kallisto_libraries

Contains kallisto quantifications for each library 

## kallisto_matrices

Matrices containing kallisto quantifications for each pairwise comparison examined

## TPM_counts

For sequences that were differentially-expressed in a pairwise comparison, contains transcript per million (TPM) counts from all libraries.