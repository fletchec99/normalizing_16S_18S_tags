#Purpose: normalize raw counts of 16S and 18S ASVs by dividing counts of proks and euks by the % passing DADA2, then converting to proportions and multiplying counts of euks by the bias against them (as user specifies).
#Required packages: Rscript, args 
#Output: This returns a file with ASV relative abundances out of (16S + 18S).
#Note: we recommend assuming a 2-fold bias against 18S sequences, which has been found with Illumina HiSeq or MiSeq data (Yeh et al. 2018)
#This script must be run from the base directory (the folder that contains 02-PROKs/ and 02-EUKs/)
#Author: Colette Fletcher-Hoppe 
#Ver 12.7.20

#Make sure input arguments (bias against 18S) is provided
args=as.numeric(commandArgs(TRUE))

if(length(args)==0){
  stop("Please specify the bias against 18S sequences, e.g. '$ ./normalize_16S_18S_ASV_proportions.R <bias>'. For Illumina HiSeq or MiSeq data, we recommend assuming a 2-fold bias against 18S sequences (Yeh et al. 2018); i.e. specify '$ ./normalize_16S_18S_ASV_proportions.R 2'.n", call=FALSE)
}

if(length(args)>1){
  stop("Please only specify one number for the bias against 18S sequences.n", call=FALSE)
}

#1. Calculate percent of reads that passed DADA2 denoising for both proks and euks -----
print("1. Calculating DADA2 stats")
proks_stats <- read.table("02-PROKs/03-DADA2d/stats.tsv", sep="\t", header=T, stringsAsFactors = F)
proks_stats$perc.passed=proks_stats$non.chimeric/proks_stats$input

euks_stats <- read.table("02-EUKs/08-DADA2d/stats.tsv", sep="\t", header=T, stringsAsFactors = F)
euks_stats$perc.passed=euks_stats$non.chimeric/euks_stats$input

#2. Normalize ASV counts (divide counts of ASVs/ percent passed for each sample, multiply euks ASV counts by the bias you specified)------
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

#b. Now repeat normalization for the euks, and also multiply by the bias you specified-----
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

#Divide ASV count by percent passing DADA2 for each sample and multiply by the given bias to normalize ASV counts
for(i in euks_stats$sample.id){
  euks_norm[,grep(i, colnames(euks_subs))]=args[1]*euks_subs[,grep(i, colnames(euks_subs))]/euks_stats$perc.passed[grep(i, euks_stats$sample.id)]
}


#3. Combine proks and euks tables of normalized sequencing counts-------
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

#b. Make sure normalized proks and euks data are in the same order, but OTU ID comes first------
print("3b. Re-ordering samples")
#Proks first
proks_ordered <- as.data.frame(matrix(nrow=nrow(norm_proks), ncol=ncol(norm_proks)))
proks_ordered[,1]=norm_proks$OTU_ID 
colnames(proks_ordered)[1]="OTU_ID"
proks_ordered[,2:ncol(proks_ordered)]=norm_proks[,order(colnames(norm_proks))][-grep("OTU_ID", colnames(norm_proks)[order(colnames(norm_proks))])]
colnames(proks_ordered)[2:ncol(proks_ordered)]=colnames(norm_proks)[order(colnames(norm_proks))][-grep("OTU_ID", colnames(norm_proks)[order(colnames(norm_proks))])]

#Now euks second
euks_ordered <- as.data.frame(matrix(nrow=nrow(norm_euks), ncol=ncol(norm_euks)))
euks_ordered[,1]=norm_euks$OTU_ID
colnames(euks_ordered)[1]="OTU_ID"
euks_ordered[,2:ncol(euks_ordered)]=norm_euks[,order(colnames(norm_euks))][-grep("OTU_ID", colnames(norm_euks)[order(colnames(norm_euks))])]
colnames(euks_ordered)[2:ncol(euks_ordered)]=colnames(norm_euks)[order(colnames(norm_euks))][-grep("OTU_ID", colnames(norm_euks)[order(colnames(norm_euks))])]

#c. Merge normalized proks and euks data-----
print("3c. Merging proks and euks data")
norm_ordered_proks_euks <- rbind(proks_ordered, euks_ordered, stringsAsFactors=F)

#4. Convert normalized sequence counts to proportions-----
print("4. Converting normalized ASV counts to proportions")

#Sum up the ASV counts in each sample (16S + 18S)
colsums <- c()
for(i in c(2:ncol(norm_ordered_proks_euks))){
  colsums <- c(colsums, sum(norm_ordered_proks_euks[,i]))
}

#Set up a new dataframe for relative abundance data
relabun_proks_euks=as.data.frame(matrix(nrow=nrow(norm_ordered_proks_euks), ncol=ncol(norm_ordered_proks_euks)))
colnames(relabun_proks_euks)=colnames(norm_ordered_proks_euks)
relabun_proks_euks$OTU_ID=norm_ordered_proks_euks$OTU_ID

#Write in normalized data
for(i in c(2:ncol(relabun_proks_euks))){
  relabun_proks_euks[,i]=norm_ordered_proks_euks[,i]/colsums[(i-1)]
}

#5. Add in taxonomy and write it out!-----
print("5. Writing out the final file!")

#Write in taxonomy 
#First write in a dummy taxonomy
relabun_proks_euks$taxonomy="SAR--11"
#Then write over it with real taxonomy from proks_data and euks_data
relabun_proks_euks$taxonomy[c(1:nrow(proks_ordered))]=proks_data$taxonomy 
relabun_proks_euks$taxonomy[c((nrow(proks_ordered)+1):nrow(relabun_proks_euks))]=euks_data$taxonomy

#Write out file
write.table("ASV_proportions_16S_18S.tsv", x=relabun_proks_euks, sep="\t", row.names=F, quote=F)

print("Script complete! Please find the relative abundances of your 16S and 18S ASVs (out of (16S + 18S) sequences) in this directory, the file is called 'ASV_proportions_16S_18S.tsv.' Thank you!")
