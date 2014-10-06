rm(list=ls(all=TRUE))
library(WGCNA)
library(igraph)
library(compiler)
library(MASS)

source("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/MakeGrapObject.R")
source("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/plotg.R")
power<-20

#####################################
#Start the analysis
#####################################

#set the wd
wdir<-"/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/WGCNA\ networks"
#wdir<-getwd()
#setwd(paste(wdir,power,sep=""))
setwd(wdir)

#grab the relavent files and load in the data
files <- list.files(getwd())
Rfiles <- files[grep(glob2rx("*.R"), files)]

for ( i in Rfiles) {
  load(i)
}#for

###############################################################
#Calculate adjacency
###############################################################

#for At dom graph
matAt<-adjacency(Expr[[1]]$data)
write.matrix(matAt, file="At.adj.matrix.txt")
#for Dt dom graph...
matDt<-adjacency(Expr[[2]]$data)
write.matrix(matDt, file="Dt.adj.matrix.txt")


##
#Use a perl program to trim the data, removing weak edges and duplicates etc
##

system(command="perl /Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/edgeList.pl At.adj.matrix.txt At.adj.edges.txt")
system(command="perl /Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/edgeList.pl Dt.adj.matrix.txt Dt.adj.edges.txt")

#load the data back in
edgeTestAt<-read.table("At.adj.edges.txt", sep = "\t", header = FALSE)
edgeTestDt<-read.table("Dt.adj.edges.txt", sep = "\t", header = FALSE)

######################################################
#Calculate the networks
######################################################

######
# Create graph
######

#define the colors
colorsAt<-labels2colors(At.net$colors)
names(colorsWild)<-colnames(Expr[[2]]$data)

#make the grpah object
gAt<-MakeGraphObject(edgeTestAt, geneNames=colnames(Expr[[2]]$data), geneCols=colorsAt,minDegree=5)
gDt<-MakeGraphObject(edgeTestDt, geneNames=colnames(Expr[[2]]$data), geneCols=colorsAt,minDegree=5)
#plotg(gWild,lwd = 0.1 , cex = 0.65, bg = "white",seg.col="darkgrey")
#plotg(gDom,lwd = 0.1 , cex = 0.65, bg = "white",seg.col="darkgrey")

#print to pdf
pdf("At_vs_Dt_igraph.pdf", height = 12, width =24)
par(mfrow=c(1,2))
plotg(gAt, lwd = 0.1 , cex = 0.9, bg = "white",seg.col="grey",col=gAt$C)
plotg(gDt, lwd = 0.1 , cex = 0.9, bg = "white",seg.col="grey",col=gDt$C)
dev.off()

####
# END
####
