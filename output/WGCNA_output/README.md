# output/WGCNA_output

Contact: Aspen Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

Last edited README: 2021-12-06

## WGCNA_output

Output from the R package WGCNA, which performs weighted gene correlation network analyses

Contains one folder for each set of libraries analyzed. Reminder that:
- cbai_transcriptomev2.0: libraries aligned to complete transcriptome
- cbai_transcriptomev4.0: libraries aligned to Tanner crab transcriptome
- hemat_transcriptomev1.6: libraries aligned to _Hematodinium_transcriptome

Each contains several sub-folders. Most of these are from early runs of WGCNA which can be discarded due to improper setup.
- amb_crabs_no_filter: discard, didn't run with enough samples
- amb_vs_elev_days_02: same as above
- amb_vs_elev_no_filter: same as above
- elev_crabs_no_filter: same as above
- **all_crabs_no_filter**: The good one. Examines ALL individual libraries.

## [transcriptome_name]/all_crabs_no_filter/

Each contains several sub-directories, according to which variable was determined as "the" variable. This creates only very minor differences in how WGCNA was set up, since it inherently examines correlation to each variable, but for the sake of repetition, I ran for each.

The important one is **day_as_variable**

## [transcriptome_name]/all_crabs_no_filter/day_as_variable/

Each contains even more sub-directories from previous runs of WGCNA where it wasn't set up correctly. 
- including_cPCR: Both the cPCR and qPCR results for _Hematodinium_ infection were used in the analysis, creating issues

- unsigned: Used an unsigned network, which is generally suboptimal compared to a signed network

- hemat_as_binary: qPCR results for _Hematodinium_ infection were set as binary (high infection vs. low infection), rather than a more  fine-tune approach, which was using the specific qPCR SQ mean

In summary, **the important pathway is WGCNA_output/[transcriptome_name]/all_crabs_no_filter/day_as_variable/*.txt (or *.png)**
