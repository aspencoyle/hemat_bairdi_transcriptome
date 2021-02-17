# hemat_bairdii_transcriptome
 
First project upon starting grad school - examining differential expression of shared hematodinium/C. bairdii libraries at various temperatures and at different time points.

Data has same libraries and same general notation as Grace Crandall's research

Date started project: 2020/10/25 (approximation)
Date started Github repo: 2020/12/11

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

**Data Sources**

*Libraries*


[Trimmed individual libraries downloaded from here](https://gannet.fish.washington.edu/Atumefaciens/20200318_cbai_RNAseq_fastp_trimming/) at 22:00 PST on 2021-02-02

[Pooled libraries downloaded from here](https://gannet.fish.washington.edu/Atumefaciens/20200414_cbai_RNAseq_fastp_trimming/) at 24:00 PST on 2021-02-02

Mapping between libraries and treatments is available [at this Google doc](https://docs.google.com/spreadsheets/d/1d17yg5F5gKKC66O8QkTIlPxljJeuX7ZsG46pkBr1lNQ/edit#gid=0)

*Transcriptome*

[Transcriptomes downloaded from here](https://owl.fish.washington.edu/halfshell/genomic-databank/). at 01:00 PST on 2021-02-03

Transcriptome checksums, along with additional information (including BLASTx annotation, GO terms annotation, etc) [available here](https://github.com/RobertsLab/resources/wiki/Genomic-Resources)

*Swiss-Prot Database*

[Swiss-Prot database manually downloaded from here](https://www.uniprot.org/uniprot/?) at 15:00 on 2021-02-09.
Downloaded an uncompressed .tab file that includes all GO terms

*Alveolata Database*
[All _Alveolata_ nucleotide sequences downloaded from here](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid33630[Organism:exp]) at 17:00 on 2021-02-16

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