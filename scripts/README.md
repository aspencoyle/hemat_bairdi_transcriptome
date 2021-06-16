Aidan Coyle, afcoyle@uw.edu
Updated 2021-06-15

### Notation

Scripts are organized with a two-digit ID number. The first corresponds to the analysis, and the second corresponds to the script number within the analysis.

Analysis 1: Putting mixed _C. bairdi_ / _Hematodinium_ transcriptome, AKA cbai_transcriptomev2.0 or cbai_hemat_transcriptomev2.0, through pipeline (kallisto -> DESeq2 -> GO-MWU)

Analysis 2: Misc analyses exploring output from Analysis 1, trying to determine what sequences are _Hematodinium_ and what sequences are _C. bairdi_.

Analysis 3: Putting a _Hematodinium_ only transcriptome, AKA hemat_transcriptomev1.6, through the same pipeline as Analysis 1.

Analysis 4: Putting a _C. bairdi_ - only transcriptome, AKA cbai_transcriptomev4.0, through the same pipeline as Analysis 1

Analysis 5: Running WGCNA to try to find correlation between either gene modules over time or between temperature treatments for various transcriptome alignments

Analysis 6: Visualizing WGCNA output by graphing expression of individual modules

Analysis 7: "Manual" gene clustering - an alternative to WGCNA, where genes are clustered by expression patterns into modules by setting distinct cut heights

Analysis 8: Using goseq to examine enrichment within individual modules

**1_1_download_libraries_run_kallisto.ipynb** 

Downloads individual and pooled libraries, downloads two transcriptomes of mixed _C. bairdi_ and _Hematodinium_, creates kallisto index files from transcriptomes, runs kallisto on libraries, and then combines output into matrix. Kallisto is quite computationally intensive, and as a result, this was done on a remote machine. Therefore, commands are copied and pasted, not run directly in the notebook.

**1_2_kallisto_to_deseq_to_accessionIDs.Rmd**

This script utilizes a series of custom-built functions for the following: 
- running a Trinity matrix of kallisto transcript counts through DESeq2
- turning DESeq2 output into newline-separated file of UniProt accessions
- creating Venn diagrams of transcript IDs and accession IDs (not needed for future steps - can skip if not wanted)

**1_3_uniprot_to_GO.Rmd**

This script is used to obtain one of the two inputs for GO-MWU (a tab-delimited 2-column table of gene IDs and GO terms). It also requires you to have a recently-downloaded copy of the SwissProt database.


**1_4_eliminate_duplicates.ipynb**

Takes the output from 13_uniprot_to_GO.R and eliminates duplicate rows, thus finalizing it for readiness in GO-MWU

**1_5_GO-MWU_prep.Rmd**

GO-MWU needs two input files - a tab-delimited 2-column table of gene IDs and GO terms with no header and one line per gene, and a 2-column CSV of gene IDs and GO terms with a header. We already created the former in 13_uniprot_to_GO.R and 14_eliminate_duplicates.ipynb. This script creates the latter file (the two-column CSV of gene IDs and GO terms with a header)

**1_6_running_GO-MWU/**

This folder contains: 

    1. an R script I created that runs GO-MWU (titled 16_running_GO-MWU)

    2. all files needed for running GO-MWU, available [here](https://github.com/z0on/GO_MWU)

    3. The specific input files for GO-MWU that we created (suboptimal structurally, but required for GO-MWU to run)

**2_1_obtaining_TPM_for_DEGs.Rmd**

Takes a file of genes that are differentially-expressed (output from DESeq2) and cross-references it with the transcriptome used to create the kallisto index. It then extracts the TPM (transcripts per million) counts from the kallisto libraries created earlier in the pipeline. This produces a single table containing the following:

- transcript and accession IDs for all inputted DEGs

- all TPM counts from individual kallisto libraries for those DEGs

**2_2_DEG_blast.ipynb**

This script takes a newline-separated file of accession IDs for all DEGs and BLASTs it (using BLASTn)against all Alveolata nucleotide sequences, and then again for all Arthropoda nucleotide sequences. The goal is to determine whether a particular DEG is more likely to be _C. bairdi_ or _Hematodinium sp._


**2_3_blastn_analysis.Rmd**

Takes the results from our BLASTs in 12_DEG_blast.ipynb and examines the e-values to determine whether sequences are more likely to have originated from _C. bairdi_ or _Hematodinium sp._

**2_4_ncbi_genome_blasts.ipynb**

Takes the closest available genome to _C. bairdi_ (_C. opilio_), and the closest available genome to _Hematodinium sp_. (_Amoebophrya sp._) and BLASTs cbai_transcriptome_v2.0 (AKA cbai_hemat_transcriptomev2.0) in a variety of ways

**Analysis 3**

These scripts perform the same processes as Analysis 1. However, while Analysis 1 analyzes cbai_transcriptome_v2.0 - which is not filtered by taxa - these scripts analyze hemat_transcriptome_v1.6, which takes all reads and filters to only include those from Alveolata (the taxa that contains _Hematodinium_). 

Steps are broadly the same, with the addition of 3_0_hemat1.6_indexcreation.ipynb. For cbai_transcriptomev2.0, a BLASTx index matching sequences to annotation IDs had already been created. However, for hemat_transcriptomev1.6, no such index existed. Therefore, one was created. To keep notation consistent between Analysis 3 and Analysis 1, this was put as Script 0 for this analysis.

**Analysis 4**

These scripts perform the same processes as Analysis 1. However, while Analysis 1 analyzes cbai_transcriptome_v2.0 - which is not filtered by taxa - these scripts analyze cbai_transcriptomev4.0, which takes all reads and filters to only include those from _Chionoecetes_.

Steps are broadly the same, with the addition of 40_cbai5.0_indexcreation.ipynb. For cbai_transcriptomev2.0, a BLASTx index matching sequences to annotation IDs had already been created. However, for cbai_transcriptomev4.0, no such index existed. Therefore, one was created. To keep notation consistent between Analysis 4 and Analysis 1, this was put as Script 0 for this analysis.

**Analysis 5**

Each script contains a different analysis using WGCNA. The first portion of the name refers to the transcriptome that each library is aligned to - either cbai_transcriptomev4.0 (in which case, only _Chionoecetes_ genes are compared), hemat_transcriptomev1.6 (for which only _Hematodinium_ genes are compared), or cbai_transcriptomev2.0, which is unfiltered by taxa.

The second portion of the name refers to the comparison being made. If it is simply AmbCrabs or ElevCrabs, then it examines expression among crabs in that treatment group (ambient or elevated) throughout the course of the experiment. If it is instead AmbVsElev, it examines expression between treatment groups (Ambient and Elevated), though also takes both day and crab into account.

**Analysis 6** 

Here, modules created in Analysis 5 are graphed to examine expression over time and between treatments

**Analysis 7**

In each file, expression over time is manually clustered, in an alternative to WGCNA, where genes are clustered by expression patterns into modules by setting distinct cut heights

**hematodinium_analysis_functions.R**

Contains a variety of functions created for the purpose of this project. Functions are called inside other R and Rmd scripts.