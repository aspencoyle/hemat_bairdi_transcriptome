# DESeq2_output

Contact: Aidan Coyle, afcoyle@uw.edu
Roberts Lab, UW-SAFS

Last edited README: 2021-03-04

Contains folders with output from the R package DESeq2. Each folder represents one pairwise comparison of different libraries.

## Filenames within folders:

**AllGenes.txt**: Columns are as follows:
- baseMean: mean of normalized counts for all samples
- log2FoldChange: log2 fold change (MLE): condition treated vs. untreated
- lfcSE: standard error: condition treated vs untreated
- stat: Wald statistic: condition treated vs. untreated
- pvalue: Wald test p-value: condition treated vs. untreated
- padj: BH adjusted p-values

**AllGenes_wcols.txt**: Same as AllGenes.txt, but with column headers

**allres_MAplot.png**: MA plot of all sequences

**allres_shrunken_MAplot.png**: MA plot of all sequences shrunken using apeglm

**DEGlist.txt**: All differentially-expressed genes, as determined by the p-value set when running DESeq2. Columns are same as AllGenes.txt

**DEGlist_wcols.txt**: Same as DEGlist.txt, but with column headers

**dispersion_estimates.png**: Dispersion estimates 

**normalizedcts_v_log2foldchange**: Same as MA plot

**PCA_plot**: principal component analysis plot

**res05_MAplot**: MA plot, but only with sequences with an adjusted p-value <= 0.05