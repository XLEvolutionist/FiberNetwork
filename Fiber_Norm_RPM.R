##
#A script to normalize Fiber expression data
#Simon Renny-Byfield, Iowa State University, Aug 2014
##

#set the working directory
setwd("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis")

#load in the data
Counts.A1<-read.table("A1_counts.txt",header = T, sep = "\t")

Counts.AD1<-read.table("AD1_counts.txt",header = T, sep = "\t")


#set up a data frame with sample info
sampleInfo<-read.csv("SampleInfo.csv",
                     header = T)

#re-arrange the data a bit
row.names<-Counts.AD1[,1]
counts<-cbind(Counts.A1[,-1],Counts.AD1[,-1])
rownames(counts)<-row.names
counts<-counts[,order(match(colnames(counts),sampleInfo$sample))]

extra.Dat<-read.table("extra.Data.txt", header = T , sep = "\t")

#replace some samples with Josh's data
#Maxxa 20 DPA
counts[,31]<-extra.Dat[,7]
#Yuc 10 DPA
counts[,45]<-extra.Dat[,11]

#find the number of mapped reads in each library
colSums(counts)
#remove teh rows with zero across the board
rowTotals<-rowSums(counts)
#counts<-counts[rowTotals>0,]

#modify the order of columns ot match sampleInfo
counts<-counts[,order(match(colnames(counts),sampleInfo$sample))]

#save the raw counts into an R object
save(counts, file="Fiber_counts.R")

#turn the counts into RPM
counts_per_RPM<-sweep(counts,2,colSums(counts), FUN="/")
counts_per_RPM<-(counts_per_RPM)*1e6

exp.Dat<-counts_per_RPM
exp.Dat<-as.matrix(exp.Dat)

#save the count data
save(exp.Dat, file="FiberTranscriptome_RPM.R")
