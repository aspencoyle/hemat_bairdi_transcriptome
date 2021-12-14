#############################
# Aidan Coyle, afcoyle@uw.edu
# Roberts lab, UW-SAFS
# 2021-03-15
#############################

# This file runs GO-MWU. Commands copied from the GO-MWU.R file
# in the GO-MWU Github repository, available at https://github.com/z0on/GO_MWU

# We will run GO-MWU twice - once for each comparison



library(ape)
library(tidyverse)

# I know, I know - suboptimal to have setwd(). 
# I assume that all scripts start from within your scripts directory, 
# so this just indicates that you need to move to the same subdirectory
# as all other GO-MWU files - both data files and analysis files.
# GO-MWU doesn't cooperate otherwise.
setwd("11_4_running_GO-MWU")

#### ALL LIBRARIES ------------------------------------

#### COMPLETE TRANSCRIPTOME MODULES ------------------------

#### GO-MWU Run 1: Complete Transcriptome, All Libraries, Blue Module --------------------------

# Edit these to match your data file names: 
input="cbaiv2.0_all_crabs_blue_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai2.0_all_crabs_blue_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO term at 10% FDR

# Not enough to graph, stopping analysis here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv2.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev2.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

#### GO-MWU Run 2: Complete Transcriptome, All Libraries, Cyan Module --------------------------

# Edit these to match your data file names: 
input="cbaiv2.0_all_crabs_cyan_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai2.0_all_crabs_cyan_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO terms at 10% FDR

# Stopping analysis here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv2.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev2.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))


#### GO-MWU Run 3: Complete Transcriptome, All Libraries, Salmon Module --------------------------

# Edit these to match your data file names: 
input="cbaiv2.0_all_crabs_salmon_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai2.0_all_crabs_salmon_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=5,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results

# 0 GO terms at 10% FDR

# Stopping analysis here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv2.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev2.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))


#### GO-MWU Run 4: Complete Transcriptome, All Libraries, Yellow Module --------------------------

# Edit these to match your data file names: 
input="cbaiv2.0_all_crabs_yellow_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai2.0_all_crabs_yellow_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO terms at 10% FDR

# Gives an error because not enough genes present in the module. 
# Stopping analysis here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv2.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev2.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

#### TANNER CRAB TRANSCRIPTOME MODULES ---------------------------

#### GO-MWU Run 5: Tanner crab Transcriptome, All Libraries, Brown Module --------------------------

# Edit these to match your data file names: 
input="cbaiv4.0_all_crabs_brown_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai4.0_all_crabs_brown_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO terms at 10% FDR

# Ending here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv4.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev4.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

#### GO-MWU Run 6: Tanner crab Transcriptome, All Libraries, Black Module --------------------------

# Edit these to match your data file names: 
input="cbaiv4.0_all_crabs_black_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai4.0_all_crabs_black_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO terms at 10% FDR

# Ending here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv4.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev4.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))


#### PARASITE TRANSCRIPTOME MODULES ------------------------

#### GO-MWU Run 8: Parasite Transcriptome, All Libraries, Brown Module --------------------------

# Edit these to match your data file names: 
input="hematv1.6_all_crabs_brown_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="hemat1.6_all_crabs_brown_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 7 GO terms at 10% FDR

grDevices::windows()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=0.05,  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs


# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("hematv1.6_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/hemat_transcriptomev1.6/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_hemat.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))


#### GO-MWU Run 9: Parasite Transcriptome, All Libraries, Pink Module --------------------------

# Edit these to match your data file names: 
input="hematv1.6_all_crabs_pink_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="hemat1.6_all_crabs_pink_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 2 GO terms at 10% FDR

grDevices::windows()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=0.05,  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs


# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("hematv1.6_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/hemat_transcriptomev1.6/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_hemat.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))


#### AMBIENT AND LOWERED TREATMENT LIBRARIES ----------------------------

#### GO-MWU Run 10: Complete Transcriptome, Amb + Low, Blue Module --------------------------

# Edit these to match your data file names: 
input="cbaiv2.0_amb_low_blue_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai2.0_amb_low_blue_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 1 GO terms at 10% FDR

# Can't graph, but note that we do have a significant GO term!

# Checked output file, looks like it's macromolecule modification (padj = 0.075)

# Stopping analysis here, because can't graph

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv2.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev2.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

#### GO-MWU Run 11: Tanner crab Transcriptome, Amb + Low, Red Module --------------------------

# Edit these to match your data file names: 
input="cbaiv4.0_amb_low_red_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="cbai4.0_amb_low_red_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 0 GO terms at 10% FDR

# Ending here

# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("cbaiv4.0_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/cbai_transcriptomev4.0/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_cbai.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

#### GO-MWU Run 8: Parasite Transcriptome, Amb + Low, Green Module --------------------------

# Edit these to match your data file names: 
input="hematv1.6_amb_low_green_module_kMEs.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="hemat1.6_amb_low_green_module_GOIDs_norepeats.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="C:/Users/acoyl/Documents/GradSchool/RobertsLab/Tools/perl/bin/perl.exe", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

# --------------- Results
# 8 GO terms at 10% FDR

grDevices::windows()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=0.05,  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.8 # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs


# Move the 3 files we created to a permanent folder, since GO-MWU automatically puts them in
# the same folder you run the script in

file_loc <- input %>%
  str_remove("hematv1.6_") %>%
  str_remove("_kMEs.csv")

filepath <- paste0("../../output/GO-MWU_output/WGCNA_modules/hemat_transcriptomev1.6/",
                   file_loc, "/")

# We're encountering some issues with a double-named filename (something like "dissim_BP_cbai2.0_amb0217_elev0_low0_vs_elev2_l2FC.csv_cbai2.0_amb0217_elev0_low0_vs_elev2_GOIDs_norepeats.txt")
# It's a problem with the function, but I don't want to touch the prebuilt GO-MWU stuff,
# so I'm just going to remove part of the name
files <- list.files(getwd(), pattern = "BP_")
newfiles <- gsub("\\.csv_hemat.*", ".txt", files)
file.rename(files, newfiles)

file.copy(list.files(getwd(), pattern = "BP_"), filepath)
file.remove(list.files(getwd(), pattern = "BP_"))

