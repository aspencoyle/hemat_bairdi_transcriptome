# hemat_bairdi_transcriptome

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

Last edited README: 2021-03-04

Examining differential expression of shared _Hematodinium_/_C. bairdi_ libraries over temperature and time.

All libraries originate from Grace Crandall's research. The same libraries and general notation is used here. 

# Tool Information and Software Versions

### kallisto: version 0.46.0

### Trinotate: version 3.2.1

### Jupyter Lab: version 2.1.5

### Jupyter Notebook: version 6.0.3

### Operating system:
```
Windows 10, x64 (OS build: 18363.1256)
Running WSL version 1
```


### R packages and info:
```
R version 4.0.4 (2021-02-15)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
[10] base     

other attached packages:
 [1] VennDiagram_1.6.20          futile.logger_1.4.3         ape_5.4-1                  
 [4] vsn_3.58.0                  forcats_0.5.1               stringr_1.4.0              
 [7] dplyr_1.0.4                 purrr_0.3.4                 readr_1.4.0                
[10] tidyr_1.1.3                 tibble_3.1.0                ggplot2_3.3.3              
[13] tidyverse_1.3.0             DESeq2_1.30.1               SummarizedExperiment_1.20.0
[16] Biobase_2.50.0              MatrixGenerics_1.2.1        matrixStats_0.58.0         
[19] GenomicRanges_1.42.0        GenomeInfoDb_1.26.2         IRanges_2.24.1             
[22] S4Vectors_0.28.1            BiocGenerics_0.36.0         apeglm_1.12.0  
```
## Folders:

data: unmodified data files

graphs: output from various R packages ran as part of analyses. Does include several tables

output: output from analyses

paper: publication

scripts: code to run analyses. Does include some output within the GO-MWU directory, as it is required for GO-MWU to run.

## Subfolders:

Many folders in output/ and graphs/ are split into subdirectories labeled cbaihemat_transcriptomev2.0 and hemat_transcriptomev1.6. These correspond to analyses performed using [cbai_transcriptomev2.0](https://robertslab.github.io/sams-notebook/2020/05/02/Transcriptome-Assembly-C.bairdi-All-RNAseq-Data-Without-Taxonomic-Filters-with-Trinity-on-Mox.html) and [hemat_transcriptomev1.6](https://robertslab.github.io/sams-notebook/2021/03/08/Transcriptome-Assembly-Hematodinium-Transcriptomes-v1.6-and-v1.7-with-Trinity-on-Mox.html)