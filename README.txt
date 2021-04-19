##Purpose:
#This script merges 16S and 18S ASV tables generated through the qiime2 pipeline modified by Jesse McNichol by normalizing raw ASV counts against 1) the percent of reads passing DADA2 and 2) the sequencing platform's bias against longer 18S sequences.

##Inputs:
#The script automatically imports raw counts of the 16S and 18S ASVs ("02-PROKs/10-exports/all-16S-seqs.with-tax.tsv" and "02-EUKs/15-exports/all-18S-seqs.with-PR2-tax.tsv", respectively).
#User must specify the bias against 18S sequences. We have found that a bias of 2x correctly normalizes our mixed mock communities of both 16S and 18S sequences (Yeh et al., 2018). 2x is therefore the default bias against 18S sequences in this script. However, samples from an environment with a higher proportion of long 18S sequences may have a larger bias against 18S.
#User must also specify the prefix of the output file names. 

##Outputs:
#Two .tsv files will be generated, the normalized counts of 16S and 18S sequences and the counts converted to proportions (out of 16S+18S total sequences).

##Useage: 
#$ Rscript normalize_16S_18S.R --bias <number> --outputfile <name>

##Dependent packages: 
#This script requires Rscript, base R, and optparse. 
#Dependent packages may be found within the conda .yml file in this repository. 

##Version: 04.19.2021
##Author: Colette Fletcher-Hoppe, with help from Jesse McNichol.
