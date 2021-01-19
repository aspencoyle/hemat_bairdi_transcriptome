# hemat_bairdii_transcriptome
 
First project upon starting grad school - examining differential expression of shared hematodinium/C. bairdii libraries at various temperatures and at different time points.

Data has same libraries and same general notation as Grace Crandall's research

Date started project: 2020/10/25 (approximation)
Date started Github repo: 2020/12/11

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

**Data Sources**

*Libraries*

[Individual and pooled libraries downloaded from here](https://owl.fish.washington.edu/nightingales/C_bairdi/). Also available - along with which libraries match to which treatments - at [this Google doc](https://docs.google.com/spreadsheets/d/1d17yg5F5gKKC66O8QkTIlPxljJeuX7ZsG46pkBr1lNQ/edit#gid=0)

*Transcriptome*

[Transcriptome v3.0 downloaded from here](https://owl.fish.washington.edu/halfshell/genomic-databank/cbai_transcriptome_v3.0.fasta). 

Further information for transcriptome 3.0 (including BLASTx annotation, GO terms annotation, etc) [available here](https://github.com/RobertsLab/resources/wiki/Genomic-Resources)

**Tool Information**

R: 
```
platform       x86_64-w64-mingw32          
arch           x86_64                      
os             mingw32                     
system         x86_64, mingw32             
status                                     
major          4                           
minor          0.3                         
year           2020                        
month          10                          
day            10                          
svn rev        79318                       
language       R                           
version.string R version 4.0.3 (2020-10-10)
nickname       Bunny-Wunnies Freak Out    
```

kallisto: version 0.46.0

Trinotate: version 3.2.1

Operating system:
```
Windows 10, x64 (OS build: 18363.1256)
Running WSL version 1
```

R packages:
```
apeglm: 1.12.0
BiocManager: 1.30.10
DESeq2: 1.30.0
tidyverse: 1.3.0
topGO: 2.42.0
VennDiagram: 1.6.20
vsn: 3.58.0
```