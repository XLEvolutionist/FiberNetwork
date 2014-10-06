rm(list=ls(all=TRUE))
library(ggplot2)
library(WGCNA)
library(plotrix)
library(RColorBrewer)
library(vioplot)
library(gtools)

#source useful functions
source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/moduleCompare.R")
source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/smoothConnectivityPlot.R")
source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/labeledHeatMap2.R")

##
# Set up some parameters
##

#the minimum num of genes befroe considering a collection as a module
minModuleSize = 100
#merge any modules with greater than 1-mergeCutHeight similarity
mergeCutHeight = 0.20
#the power used to calculate dissimilarity
power = 20
#whether to consider just absolute correlation ("unsigned") or signed ("signed") correlation
TOMType = "signed"
#the minimum level of expression in RPM
minNormCount<-2
#the number of top hub genes to report in each module.
numGenes<-50

#set the wd
wdir<-setwd(paste("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/RPM_sft_", power , sep=""))
setwd(paste("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/RPM_sft_", power , sep=""))

#grab the relavent files and load in the data
files <- list.files()
Rfiles <- files[grep(glob2rx("*.R"), files)]
Rfiles2 <- files[grep(glob2rx("*.RData"), files)]

for ( i in c(Rfiles,Rfiles2) ) {
  load(i)
}
#load in the sampleInfo
sampleInfo<-read.csv("SampleInfo.csv")

#compare to co-expression networks
pdf(file="compare_modules.pdf")
moduleCompare(x=AD1.wild.net,y=AD1.dom.net,xlab="AD1 wild", ylab = "AD1 dom")
dev.off()

#generate some colors for the genes, using AD1 domesticated and reference
mergedColorsAD1.dom = labels2colors(AD1.dom.net$colors)
colorsAD1<-names(table(mergedColorsAD1.dom))
#generate some colors for genes using AD1 wild and reference
mergedColorsAD1.wild = labels2colors(AD1.wild.net$colors)
colorsAD1.wild<-names(table(mergedColorsAD1.wild))
#gene colors using module membership in A1 dom
mergedColorsA1.dom = labels2colors(A1.dom.net$colors)
colorsA1.dom<-names(table(mergedColorsA1.dom))
#gene colors using module membership in A1 wild
mergedColorsA1.wild = labels2colors(A1.wild.net$colors)
colorsA1.wild<-names(table(mergedColorsA1.wild))

##
# now look at module mebership (kME) for AD1 dom and AD1 wild,
# using the modules definition of AD1 wild
##

PCs.A1.Dom    = moduleEigengenes(multiExpr[[1]]$data,  colors=mergedColorsA1.wild) 
ME_A1.Dom    = PCs.A1.Dom$eigengenes

PCs.A1.wild    = moduleEigengenes(multiExpr[[2]]$data,  colors=mergedColorsA1.wild) 
ME_A1.wild    = PCs.A1.wild$eigengenes

PCs.AD1.Dom    = moduleEigengenes(multiExpr[[3]]$data,  colors=mergedColorsAD1.wild) 
ME_AD1.Dom    = PCs.AD1.Dom$eigengenes

PCs.AD1.wild    = moduleEigengenes(multiExpr[[4]]$data,  colors=mergedColorsAD1.wild) 
ME_AD1.wild    = PCs.AD1.wild$eigengenes

#calculte kME for all genes, but with AD1 wild as a reference for module membership
geneModuleMembershipA1.Dom <- signedKME( multiExpr[[1]], ME_A1.wild )
geneModuleMembershipA1.wild <- signedKME( multiExpr[[2]], ME_A1.wild )
geneModuleMembershipAD1.Dom <- signedKME( multiExpr[[3]], ME_AD1.wild )
geneModuleMembershipAD1.wild <- signedKME( multiExpr[[4]], ME_AD1.wild )


#calcualte mean connectivity per gene in each co-expression network
soft.Connectivity.A1.Dom<-softConnectivity(multiExpr[[1]]$data,type=TOMType, power = power)
soft.Connectivity.A1.wild<-softConnectivity(multiExpr[[2]]$data,type=TOMType,power = power)
soft.Connectivity.AD1.Dom<-softConnectivity(multiExpr[[3]]$data,type=TOMType,power = power)
soft.Connectivity.AD1.wild<-softConnectivity(multiExpr[[4]]$data,type=TOMType,power = power)

###
#plot MF vs LF for connectivity
###
names(soft.Connectivity.A1.Dom)<-colnames(multiExpr[[1]]$data)
smoothConnectivityPlot(y=soft.Connectivity.AD1.wild[match(LF.dge.names,colnames(multiExpr[[1]]$data))],
                          x=soft.Connectivity.AD1.wild[match(MF.dge.names,colnames(multiExpr[[1]]$data))], nrpoints=9000,colors=mergedColorsAD1.wild,
                       xlab="MF kME", ylab="LF kME")


#a useful color pallette for plotting later
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))




#draw a plot of connectivity for each gene within wild and domesticated networks
pdf("Connectivity.pdf", height = 10 , width = 9)
par(mfrow = c(2,2))
smoothConnectivityPlot(y=soft.Connectivity.AD1.Dom,x=soft.Connectivity.AD1.wild, nrpoints=9000,colors=mergedColorsAD1.wild,
                        xlab="AD1 wild", ylab="AD1 domesticated")
smoothConnectivityPlot(y=soft.Connectivity.A1.Dom,x=soft.Connectivity.A1.wild, nrpoints=9000,colors=mergedColorsA1.wild,
                        xlab = "A1 wild", ylab = "A1 domesticated")
smoothConnectivityPlot(y=soft.Connectivity.AD1.Dom,x=soft.Connectivity.AD1.wild, nrpoints=9000,colors=mergedColorsAD1.wild,
                       xlab="AD1 wild", ylab="AD1 domesticated",points=FALSE)
smoothConnectivityPlot(y=soft.Connectivity.A1.Dom,x=soft.Connectivity.A1.wild, nrpoints=9000,colors=mergedColorsA1.wild,
                       xlab = "A1 wild", ylab = "A1 domesticated",points=FALSE)
dev.off()


MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembershipAD1.Dom),dim(multiExpr[[3]]$data)[1]); 

colnames(MMPvalue1)=paste("AD1_Dom_",colnames(geneModuleMembershipAD1.Dom),".pval",sep="");
colnames(geneModuleMembershipAD1.Dom)=paste("AD1_Dom_",colnames(geneModuleMembershipAD1.Dom),".cor",sep=""); 

Gene <- rownames(t(multiExpr[[3]]$data))
kMEtable1<-NULL
kMEtable1<-cbind(Gene,mergedColorsAD1.dom)


kMEtable1 <- cbind(kMEtable1, geneModuleMembershipAD1.Dom, MMPvalue1)

colnames(kMEtable1)=c("Gene","Module",colorsAD1, colnames(MMPvalue1))

write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembershipAD1.wild),dim(multiExpr[[3]]$data)[1]); 

colnames(geneModuleMembershipAD1.wild)=paste("AD1_wild",colnames(geneModuleMembershipAD1.wild),".cor",sep=""); 
colnames(MMPvalue2)=paste("AD1_wild",colnames(geneModuleMembershipAD1.wild),".pval",sep="");

kMEtable2  = cbind(Gene,mergedColorsAD1.dom)
kMEtable2 = cbind(kMEtable2, geneModuleMembershipAD1.wild, MMPvalue2)
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

#print out some kME plots, with kME for all genes in all modules...
pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=11)
par(mfrow=c(2,3))
for (c in 1:length(colorsAD1.wild)){
  verboseScatterplot(geneModuleMembershipAD1.wild[,c],geneModuleMembershipAD1.Dom[,c],main=colnames(geneModuleMembershipAD1.wild)[c],
                     xlab="kME in wild",ylab="kME in dom", cex.main = 0.9, cex.axis = 0.99, cex = 0.95)
  smoothScatter(geneModuleMembershipAD1.wild[,c],geneModuleMembershipAD1.Dom[,c],xlab="kME in wild",ylab="kME in dom", 
                      bandwidth = 0.01,colramp=colorRampPalette(my.cols),main=colnames(geneModuleMembershipAD1.wild)[c])
  vioplot(geneModuleMembershipAD1.wild[,c],geneModuleMembershipAD1.Dom[,c],names=c("kME wild","kME domesticated"), 
              col = colorsAD1.wild[c])
}; dev.off()

#now kME for within modules
pdf("inModule_kMEtable2_vs_kMEtable1_signed.pdf",,height=8,width=11)
par(mfrow=c(2,3))
for (c in 1:length(colorsAD1.wild)){
  inMod = mergedColorsAD1.wild== colorsAD1.wild[c]
  verboseScatterplot(geneModuleMembershipAD1.wild[inMod,c],geneModuleMembershipAD1.Dom[inMod,c],main=colorsAD1.wild[c],
                     xlab="kME in wild",ylab="kME in dom", cex.main = 0.9, cex.axis = 0.99, cex = 0.95)
  smoothScatter(geneModuleMembershipAD1.wild[inMod,c],geneModuleMembershipAD1.Dom[inMod,c],xlab="kME in wild",ylab="kME in dom", 
                bandwidth = 0.01,colramp=colorRampPalette(my.cols),main=colorsAD1.wild[c])
  vioplot(geneModuleMembershipAD1.wild[inMod,c],geneModuleMembershipAD1.Dom[inMod,c], names=c("kME wild","kME domesticated"),
              col=colorsAD1.wild[c])
}; dev.off()

#now look at Module Eigen gene value in wild vs domesticated across all DPA
pdf("eigenGene_expression.pdf",,height=8,width=11)
par(mfrow=c(2,2))
for (mod in 1:length(colorsAD1.wild)){
  eigen<-NULL
  boxLabel<-NULL
  #for each developmental stage
  for (dev in sort(unique(sampleInfo$DPA))) {
      if ( dev == 5 ) { dev2<-"05"}
      else { dev2<-dev }
      eigen<-c(eigen,ME_AD1.Dom[sampleInfo$DPA[25:36] == dev,mod],ME_AD1.wild[sampleInfo$DPA[37:48] == dev,mod])
      boxLabel<-c(boxLabel,c(rep(paste(dev2,"dom"),3),rep(paste(dev2,"wild"),3)))
  }
  verboseBoxplot(eigen, g=boxLabel, notch = FALSE,col=colorsAD1.wild[mod], 
                 ylab = "eigenGene value", xlab = "", cex.axis = 0.9, main = colorsAD1.wild[mod])
}; dev.off()

#find the genes with the most negative kME, hub genes?
NegNames<-NULL
NegGenesKME = NULL
for (c in 1:length(colorsAD1)){
  genes1 <-  Gene[order(geneModuleMembershipAD1.Dom[,c])][1:numGenes]
  data1 <-  geneModuleMembershipAD1.Dom[order(geneModuleMembershipAD1.Dom[,c]),c][1:numGenes]
  genes2 <-  Gene[order(geneModuleMembershipAD1.wild[,c])][1:numGenes]
  data2 <-  geneModuleMembershipAD1.wild[order(geneModuleMembershipAD1.wild[,c]),c][1:numGenes]
  NegGenesKME = cbind(NegGenesKME,genes1,data1,genes2,data2)
  NegNames <- c(NegNames,c(paste("genes1",colorsAD1[c],sep=""),paste("data1",colorsAD1[c],sep=""),
                           paste("genes2",colorsAD1[c],sep=""),paste("data2",colorsAD1[c],sep="")))
}; colnames(NegGenesKME)<-NegNames
NegGenesKME

save(NegGenesKME,file="ref_AD1_dom_vs_AD1_wild_neg_kME.Data.R")

#find the genes with the most postive kME, hub genes?
PosNames<-NULL
PosGenesKME = NULL
for (c in 1:length(colorsAD1)){
  genes1 <-  Gene[rev(order(geneModuleMembershipAD1.Dom[,c]))][1:numGenes]
  data1 <-  geneModuleMembershipAD1.Dom[rev(order(geneModuleMembershipAD1.Dom[,c])),c][1:numGenes]
  genes2 <-  Gene[rev(order(geneModuleMembershipAD1.wild[,c]))][1:numGenes]
  data2 <-  geneModuleMembershipAD1.wild[rev(order(geneModuleMembershipAD1.wild[,c])),c][1:numGenes]
  PosGenesKME = cbind(PosGenesKME,genes1,data1,genes2,data2)
  PosNames <- c(PosNames,c(paste("genes1",colorsAD1[c],sep=""),paste("data1",colorsAD1[c],sep=""),
                           paste("genes2",colorsAD1[c],sep=""),paste("data2",colorsAD1[c],sep="")))
}; colnames(PosGenesKME)<-PosNames
PosGenesKME

save(PosGenesKME,file="ref_AD1_dom_vs_AD1_wild_pos_kME.Data.R")
