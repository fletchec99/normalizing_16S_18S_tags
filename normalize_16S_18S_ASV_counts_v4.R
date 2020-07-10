#!/usr/bin/Rscript

#Purpose: normalize raw counts of 16S and 18S ASVs by doubling counts of euks and dividing counts of proks and euks by the % passing DADA2.
#Next, transform to proportions with Jesse's python script
#Note: this only works with Illumina HiSeq or MiSeq data, which have a 2-fold bias against 18S sequences (Yeh et al. 2018)
#Author: Colette Fletcher-Hoppe 
#Ver 6.26.20

#This script must be run from the base directory (the folder that contains 02-PROKs/ and 02-EUKs/)

#1. Calculate precent of reads that passed DADA2 denoising for both proks and euks -----
print("1. Calculating DADA2 stats")
proks_stats <- read.table("02-PROKs/03-DADA2d/stats.tsv", sep="\t", header=T, stringsAsFactors = F)
proks_stats$perc.passed=proks_stats$non.chimeric/proks_stats$input

euks_stats <- read.table("02-EUKs/08-DADA2d/stats.tsv", sep="\t", header=T, stringsAsFactors = F)
euks_stats$perc.passed=euks_stats$non.chimeric/euks_stats$input

#2. Normalize ASV counts (divide counts of ASVs/ percent passed for each sample, multiply euks ASV counts by 2)------
print("2. Normalizing ASV counts for proks and euks")

#a. Proks--------
#Read in ASV counts data
#Fix "#OTU ID" and read in the counts table with taxonomy 
temp <- "02-PROKs/10-exports/all-16S-seqs.with-tax.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
colnames <- colnames[-c(1:7)]
proks_data <- read.table(temp, col.names=c("OTU_ID",colnames), sep="\t", stringsAsFactors = F)
colnames(proks_data)=c("OTU_ID", colnames) #removes the "X" before colnames 

#Subset matrix without taxonomy (only numerical values)
proks_subs <- proks_data[-grep("taxonomy", colnames(proks_data))]

#Set up a new matrix for normalized data
proks_norm <- as.data.frame(matrix(nrow=nrow(proks_subs), ncol=ncol(proks_subs)))
proks_norm[,1]=proks_subs$OTU_ID
colnames(proks_norm)=colnames(proks_subs)

#Divide ASV count for each sample by percent passing, write into the new matrix
for(i in proks_stats$sample.id){
  proks_norm[,grep(i, colnames(proks_subs))]=proks_subs[,grep(i, colnames(proks_subs))]/proks_stats$perc.passed[grep(i, proks_stats$sample.id)]
}

#b. Now repeat normalization for the euks, and also double them-----
#Read in files (use file with PR2 taxonomy)
temp <- "02-EUKs/15-exports/all-18S-seqs.with-PR2-tax.tsv"
colnames <- scan(text=readLines(temp, 2), what="", quiet=T)
colnames <- colnames[-c(1:7)] #Keep taxonomy in for now
euks_data <- read.table(temp, col.names=c("OTU_ID",colnames), sep="\t", stringsAsFactors = F)
colnames(euks_data)=c("OTU_ID", colnames)

#Subset matrix without taxonomy 
euks_subs <- euks_data[,-grep("taxonomy", colnames(euks_data))]

#Set up a new matrix for normalized data
euks_norm <- as.data.frame(matrix(nrow=nrow(euks_subs), ncol=ncol(euks_subs)))
euks_norm[,1]=euks_subs$OTU_ID
colnames(euks_norm)=colnames(euks_subs)

#Divide ASV count by percent passing DADA2 for each sample and multiply by 2 to normalize ASV counts
for(i in euks_stats$sample.id){
  euks_norm[,grep(i, colnames(euks_subs))]=euks_subs[,grep(i, colnames(euks_subs))]*2/euks_stats$perc.passed[grep(i, euks_stats$sample.id)]
}

#3. Concatenate tables, and write out the normalized file -------
#First, we need to format a little, because this next part only works if colnames for proks and euks data are the same.

#a. Add in dummy columns for samples missing from one matrix or the other (i.e. the mocks) ------
print("3a. Add missing samples to proks and euks matrices")
#1. First, check if there are any samples missing from proks spreadsheet, and if so, write them in 
missing_from_proks <- c()
for(i in colnames(euks_norm)){
  if(i %in% colnames(proks_norm)){
  } else{
    missing_from_proks <- c(missing_from_proks, i)
  }
}

#Add in dummy columns for these samples--Must create a new matrix to do so
norm_proks <- as.data.frame(matrix(nrow=nrow(proks_norm), ncol=(ncol(proks_norm)+length(missing_from_proks))))
colnames(norm_proks)=c(colnames(proks_norm), missing_from_proks)
norm_proks[,1:ncol(proks_norm)]=proks_norm
#Now add in dummy columns
for(i in c(1:length(missing_from_proks))){
  norm_proks[,i+ncol(proks_norm)]=0
}

#2. Then, check if there are columns missing in euks spreadsheet and if so, write them in
missing_from_euks <- c()
for(i in colnames(proks_norm)){
  if(i %in% colnames(euks_norm)){
  } else {
    missing_from_euks <- c(missing_from_euks, i)
  }
}

#Add dummy columns for these samples into normalized euks data
norm_euks <- as.data.frame(matrix(nrow=nrow(euks_norm), ncol=(ncol(euks_norm)+length(missing_from_euks))))
colnames(norm_euks)=c(colnames(euks_norm), missing_from_euks)
norm_euks[,1:ncol(euks_norm)]=euks_norm
#Now add in dummy columns 
for(i in c(1:length(missing_from_euks))){
  norm_euks[,i+ncol(euks_norm)]=0
}

#b. Make sure normalized proks and euks data are in the same order, but OTU ID comes first and taxonomy comes last------
print("3b. Re-ordering samples")
#Proks first
proks_ordered <- as.data.frame(matrix(nrow=nrow(norm_proks), ncol=ncol(norm_proks)))
proks_ordered[,1]=norm_proks$OTU_ID 
colnames(proks_ordered)[1]="OTU_ID"
proks_ordered[,2:ncol(proks_ordered)]=norm_proks[,order(colnames(norm_proks))][-grep("OTU_ID", colnames(norm_proks)[order(colnames(norm_proks))])]
colnames(proks_ordered)[2:ncol(proks_ordered)]=colnames(norm_proks)[order(colnames(norm_proks))][-grep("OTU_ID", colnames(norm_proks)[order(colnames(norm_proks))])]
#Write in taxonomy 
proks_ordered$taxonomy=proks_data$taxonomy

#Now euks second
euks_ordered <- as.data.frame(matrix(nrow=nrow(norm_euks), ncol=ncol(norm_euks)))
euks_ordered[,1]=norm_euks$OTU_ID
colnames(euks_ordered)[1]="OTU_ID"
euks_ordered[,2:ncol(euks_ordered)]=norm_euks[,order(colnames(norm_euks))][-grep("OTU_ID", colnames(norm_euks)[order(colnames(norm_euks))])]
colnames(euks_ordered)[2:ncol(euks_ordered)]=colnames(norm_euks)[order(colnames(norm_euks))][-grep("OTU_ID", colnames(norm_euks)[order(colnames(norm_euks))])]
euks_ordered$taxonomy=euks_data$taxonomy

#c. Merge normalized proks and euks data-----
print("3c. Merging proks and euks data")
norm_ordered_proks_euks <- rbind(proks_ordered, euks_ordered, stringsAsFactors=F)

#d. Tweak it slightly and write it out!-----
print("3d. Writing out the final file!")
#Write out the file and read it back in with no colnames
write.table("temp_18S_16S.tsv", x=norm_ordered_proks_euks, sep="\t", row.names=F, quote=F)
final_data <- read.table("temp_18S_16S.tsv", stringsAsFactors = F, header=F, sep="\t", row.names=NULL)
#Set colnames (new header) to be "#Constructed from biom file"
colnames(final_data)[1]="#Constructed from biom file"
#Set all other colnames to be blank
colnames(final_data)[2:length(colnames(final_data))]=""
#Set data cell 1 (first colname) to be "#OTU ID"
final_data[1,1]="#OTU ID"
#Write out file
write.table("normalized_merged_16S_18S.tsv", x=final_data, sep="\t", row.names=F, quote=F)
#Remove temp file
file.remove("temp_18S_16S.tsv")

print("Script complete! Please find your normalized, merged 16S+18S abundance counts in this directory; the file is called 'normalized_merged_16S_18S.tsv'. Next run Jesse's 'transform-ESV-tsv-to-proportions.py' script on this file.")
