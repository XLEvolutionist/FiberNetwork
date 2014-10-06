##
#A script to pwrform network analysis using WGCNA, compare Homeologous networks in domesticated cotton
##

#Simon Renny-Byfield, Iowa State University, September 2014
#updated 3/10/14

#load in the libraries
library(edgeR)
library(vegan)
library(ggplot2)
library(WGCNA)
library(plotrix)
library(RColorBrewer)
library(vioplot)
library(clusterRepro)
library(impute)
library(gplots)

options(stringsAsFactors = FALSE);
power<-16
minModuleSize<-30
mergeCutHeight<-0.20
NumGenes<-50
minNormCount<-3

source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/moduleCompare.R")
source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/smoothConnectivityPlot.R")
source(file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/scripts/labeledHeatMap2.R")

setwd("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/WGCNA\ networks")
#load in the Homeologous expression data, add your own data in here
load("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/HomeoExprTrimmed.RData")
#load in the A1 data, choose wild
load("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/FiberTranscriptome_RPM.R")
A1Expr<-exp.Dat[,13:24]
A1DomExpr<-exp.Dat[,1:12]
#load in the sampleInfo
sampleInfo<-read.csv("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/SampleInfo.csv")

#now seperate datasets for At and Dt
grab<-seq(from=1,to=dim(HomeoSort$normalized)[2], by = 3)
AtExpr<-HomeoSort$normalized[,grab]
DtExpr<-HomeoSort$normalized[,grab+1]

#remove any genes where there is zoro Dt or At expression accross rthe data set
#this is neccesary id you want to make seperate modules with At and Dt data

#AtExlc<-rownames(AtExpr)[rowSums(AtExpr)==0]
#DtExlc<-rownames(DtExpr)[rowSums(DtExpr)==0]
#A1Exlc<-rownames(A1Expr)[rowSums(A1Expr)==0]
#excl<-unique(c(AtExlc,DtExlc,A1Exlc))
#AtExpr<-AtExpr[!rownames(AtExpr) %in% excl,]
#DtExpr<-DtExpr[!rownames(DtExpr) %in% excl,]
#A1Expr<-A1Expr[!rownames(A1Expr) %in% excl,]

index<-apply(cbind(AtExpr,DtExpr,A1Expr,A1DomExpr), 1, function(x){any(x<minNormCount)})
#index<-apply(cbind(AtExpr,DtExpr), 1, function(x){any(x<minNormCount)})

AtExpr<-AtExpr[index == FALSE,]
DtExpr<-DtExpr[index == FALSE,]
A1Expr<-A1Expr[index == FALSE,]
A1DomExpr<-A1DomExpr[index == FALSE,]

nSets = 3
Expr<-vector(mode="list",length=nSets)
Expr[[1]]<-list(data=as.data.frame(t(AtExpr)))
Expr[[2]]<-list(data=as.data.frame(t(DtExpr)))
Expr[[3]]<-list(data=as.data.frame(t(A1Expr)))
#Expr[[4]]<-list(data=as.data.frame(t(A1DomExpr)))
SetLabels<-c("At","Dt","A1")
names(Expr)<-SetLabels
# cluster samples
sampleTree = flashClust(dist(t(cbind(AtExpr,DtExpr,A1Expr)), method = "euclidean"))
pdf("sampleTree.pdf")
plot(sampleTree)
dev.off()

#powers = c(c(1:10), seq(from = 12, to=30, by=2))
## Call the network topology analysis function
#sft = pickSoftThreshold(t(cbind(AtExpr,DtExpr,A1Expr)), powerVector = powers, verbose = 5)
#cex1<-1
## Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file="softThresholding.pdf", height = 6 , width = 7)
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#     main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels=powers,cex=cex1,col="red");
#abline(h=power/10,col="red")
#dev.off()

checkSets(Expr)
#first construct a consensus module set
netCon = blockwiseConsensusModules(Expr, power = power,
                                   TOMType = TOMType, minModuleSize = minModuleSize,
                                   reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                                   numericLabels = TRUE, pamRespectsDendro = FALSE,
                                   saveTOMs = TRUE,
                                   #saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                                   verbose = 3,
                                   maxBlockSize=25000)

consMEs = netCon$multiMEs;
moduleLabels = netCon$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = netCon$dendrograms[[1]];

pdf("consensus_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
#save the consesnsus modules and plot
save(netCon,file="Consensus_Modules.data.R")

# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(Expr, universalColors = moduleColors);
setLabels<-c("At", "Dt" , "A1 wild", "A1 dom")
pdf("eigenGeneComparison.pdf", width = 10, height = 10)
plotEigengeneNetworks(consMEsC, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      xLabelsAngle = 90)
dev.off()

####
#Calculate modules for At
####

#load("At.modules.RData")

At.net = blockwiseModules(Expr[[1]]$data, power = power,
                              TOMType = TOMType, minModuleSize = minModuleSize,
                              reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                              verbose = 3,
                              maxBlockSize=23000)
mergedColorsAt = labels2colors(At.net$colors)
#Plot the dendrogram and the module colors underneath
pdf("At_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(At.net$dendrograms[[1]], mergedColorsAt,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="At Cluster Dendrogram")
dev.off()
table(At.net$colors)


save(At.net, file="At.modules.RData")
####
#Calcualte the same for Dt
####
#load("Dt.modules.RData")

Dt.net = blockwiseModules(Expr[[2]]$data, power = power,
                          TOMType = TOMType, minModuleSize = minModuleSize,
                          reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                          verbose = 3,
                          maxBlockSize=23000)

mergedColorsDt = labels2colors(Dt.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("Dt_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(Dt.net$dendrograms[[1]], mergedColorsDt,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="Dt Cluster Dendrogram")
dev.off()
table(Dt.net$colors)


save(Dt.net, file="Dt.modules.RData")

####
#Calcualte the same for A1
####
#load("A1.modules.RData")

A1.net = blockwiseModules(Expr[[3]]$data, power = power,
                          TOMType = TOMType, minModuleSize = minModuleSize,
                          reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                          verbose = 3,
                          maxBlockSize=23000)

mergedColorsA1 = labels2colors(A1.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("A1_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(A1.net$dendrograms[[1]], mergedColorsA1,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="A1 Cluster Dendrogram")
dev.off()
table(A1.net$colors)


save(A1.net, file="A1.modules.RData")

####
#Compare the modules
####

#compare to co-expression networks
pdf(file="compare_modules.pdf")
moduleCompare(x=At.net,y=Dt.net,xlab="At", ylab = "Dt",cex.text=0.7)
dev.off()

#compare to co-expression networks
pdf(file="compare_modulesAtvsA1.pdf")
moduleCompare(x=At.net,y=A1.net,xlab="At", ylab = "A1",cex.text=0.7)
dev.off()

####
#Look at connectivity
####

##
# now look at module mebership (kME) for AD1 dom and AD1 wild,
# using the modules definition of AD1 wild
##

PCs.At    = moduleEigengenes(Expr[[1]]$data,  colors=mergedColorsAt) 
ME_At    = PCs.At$eigengenes

PCs.Dt    = moduleEigengenes(Expr[[2]]$data,  colors=mergedColorsAt) 
ME_Dt    = PCs.Dt$eigengenes

#now for the consensus modules
PCs.At.cons    = moduleEigengenes(Expr[[1]]$data,  colors=moduleColors) 
ME_At.cons    = PCs.At.cons$eigengenes

PCs.Dt.cons    = moduleEigengenes(Expr[[2]]$data,  colors=moduleColors) 
ME_Dt.cons    = PCs.Dt.cons$eigengenes

#calculte kME for all genes, but with AD1 wild as a reference for module membership
geneModuleMembershipAt <- signedKME( Expr[[1]], ME_At  )
geneModuleMembershipDt <- signedKME( Expr[[2]], ME_At  )

PCs.Dt.Ref   = moduleEigengenes(Expr[[2]]$data,  colors=mergedColorsDt) 
ME_Dt.Ref    = PCs.Dt.Ref$eigengenes

geneModuleMembershipAt.DtRef <- signedKME( Expr[[1]], ME_Dt.Ref  )
geneModuleMembershipDt.DtRef <- signedKME( Expr[[2]], ME_Dt.Ref  )



#calcualte mean connectivity per gene in each co-expression network
soft.Connectivity.At<-softConnectivity(Expr[[1]]$data,type=TOMType, power = power)
soft.Connectivity.Dt<-softConnectivity(Expr[[2]]$data,type=TOMType, power = power)

pdf("connectivty.pdf", height = 20, width = 10)
par(mfrow=c(2,1))
smoothConnectivityPlot(y=soft.Connectivity.At,x=soft.Connectivity.Dt, nrpoints=9000,colors=mergedColorsAt,
                       xlab="Dt", ylab="At")
smoothConnectivityPlot(y=soft.Connectivity.At,x=soft.Connectivity.Dt, nrpoints=9000,colors=mergedColorsAt,
                       xlab = "Dt", ylab = "At", points = FALSE)
dev.off()

colorsAt<-names(table(mergedColorsAt))
colorsDt<-names(table(mergedColorsDt))

#a useful color pallette for plotting later
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

#print out some kME plots, with kME for all genes in all modules...
pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=11)
par(mfrow=c(2,3))

for (c in 1:length(colorsAt)){
  verboseScatterplot(geneModuleMembershipAt[,c],geneModuleMembershipDt[,c],main=colnames(geneModuleMembershipAt)[c],
                     xlab="kME in At",ylab="kME in Dt", cex.main = 0.9, cex.axis = 0.99, cex = 0.95)
  smoothScatter(geneModuleMembershipAt[,c],geneModuleMembershipDt[,c],xlab="kME in At",ylab="kME in Dt", 
                bandwidth = 0.01,colramp=colorRampPalette(my.cols),main=colnames(geneModuleMembershipAt)[c])
  vioplot(geneModuleMembershipAt[,c],geneModuleMembershipDt[,c],names=c("kME At","kME Dt"), 
          col = colorsAt[c])
}; dev.off()

#now do module by module

#now kME for within modules
pdf("inModule_kMEtable2_vs_kMEtable1_signed.pdf",,height=8,width=11)
par(mfrow=c(2,3))
for (c in 1:length(colorsAt)){
  inMod = mergedColorsAt == colorsAt[c]
  verboseScatterplot(geneModuleMembershipAt[inMod,c],geneModuleMembershipDt[inMod,c],main=colorsAt[c],
                     xlab="kME in wild",ylab="kME in dom", cex.main = 0.9, cex.axis = 0.99, cex = 0.95)
  smoothScatter(geneModuleMembershipAt[inMod,c],geneModuleMembershipDt[inMod,c],main=colorsAt[c],xlab="kME in At",ylab="kME in Dt", 
                bandwidth = 0.01,colramp=colorRampPalette(my.cols))
  vioplot(geneModuleMembershipAt[inMod,c],geneModuleMembershipDt[inMod,c], names=c("kME Dt","kME At"),
          col=colorsAt[c])
}; dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1_signedDtRef.pdf",,height=8,width=11)
par(mfrow=c(2,3))
for (c in 1:length(colorsDt)){
  inMod = mergedColorsDt == colorsDt[c]
  verboseScatterplot(geneModuleMembershipAt.DtRef[inMod,c],geneModuleMembershipDt.DtRef[inMod,c],main=colorsDt[c],
                     xlab="kME in At",ylab="kME in Dt", cex.main = 0.9, cex.axis = 0.99, cex = 0.95)
  smoothScatter(geneModuleMembershipAt.DtRef[inMod,c],geneModuleMembershipDt.DtRef[inMod,c],main=colorsDt[c],xlab="kME in At",ylab="kME in Dt", 
                bandwidth = 0.01,colramp=colorRampPalette(my.cols))
  vioplot(geneModuleMembershipAt.DtRef[inMod,c],geneModuleMembershipDt.DtRef[inMod,c], names=c("kME At","kME Dt"),
          col=colorsDt[c])
}; dev.off()




#now look at Module Eigen gene value in wild vs domesticated across all DPA
pdf("eigenGene_expression.pdf",,height=8,width=11)
par(mfrow=c(2,2))
for (mod in 1:length(colorsAt)){
  eigen<-NULL
  boxLabel<-NULL
  #for each developmental stage
  for (dev in sort(unique(sampleInfo$DPA))) {
    if ( dev == 5 ) { dev2<-"05"}
    else { dev2<-dev }
    eigen<-c(eigen,ME_At[sampleInfo$DPA[25:36] == dev,mod],ME_Dt[sampleInfo$DPA[25:36] == dev,mod])
    boxLabel<-c(boxLabel,c(rep(paste(dev2,"At"),3),rep(paste(dev2,"Dt"),3)))
  }
  verboseBoxplot(eigen, g=boxLabel, notch = FALSE,col=colorsAt[mod], 
                 ylab = "eigenGene value", xlab = "", cex.axis = 0.9, main = colorsAt[mod])
  #now try to generate a heat map for genes in each module
  #grab the expression data
   
}; dev.off()

pdf("heatmaps.pdf", height = 10, width = 10)
par(mfrow=c(2,2))
for (mod in 1:length(colorsAt)){
  At<-t(Expr[[1]]$data)
  Dt<-t(Expr[[2]]$data)
  At<-At[mergedColorsAt==colorsAt[mod],]
  Dt<-Dt[mergedColorsAt==colorsAt[mod],]
  order<-c(4,8,12,1,5,9,2,6,10,3,7,11)
  At<-At[,order]
  Dt<-Dt[,order]
  moduleHeatData<-rbind(At,Dt)
  #figure out the 
  group<-c("At 5 DPA","Dt 5 DPA","At 10 DPA","Dt 10 DPA", "At 15 DPA", "Dt 15 DPA", "At 20 DPA", "Dt 20 DPA")
  m<-cbind(rowMeans(moduleHeatData[,1:3]),rowMeans(moduleHeatData[,4:6]),rowMeans(moduleHeatData[,7:9]),rowMeans(moduleHeatData[,10:12]))
  m<-matrix(m,nrow=nrow(m)/2)
  heatmap(m,Colv = NA,Rowv=NA,labRow=NA,labCol=group, col = redgreen(10),main=colorsAt[mod])
}#for
dev.off()

#now look at Module Eigen gene value in wild vs domesticated across all DPA
colorsCons<-names(table(moduleColors))

pdf("eigenGeneCons_expression.pdf",,height=8,width=11)
par(mfrow=c(2,2))
for (mod in 1:length(colorsCons)){
  print(colorsCons[mod])
  eigen<-NULL
  boxLabel<-NULL
  #for each developmental stage
  for (dev in sort(unique(sampleInfo$DPA))) {
    if ( dev == 5 ) { dev2<-"05"}
    else { dev2<-dev }
    eigen<-c(eigen,ME_At.cons[sampleInfo$DPA[25:36] == dev,mod],ME_Dt.cons[sampleInfo$DPA[25:36] == dev,mod])
    boxLabel<-c(boxLabel,c(rep(paste(dev2,"At"),3),rep(paste(dev2,"Dt"),3)))
  }
  verboseBoxplot(eigen, g=boxLabel, notch = FALSE,col=colorsCons[mod], 
                 ylab = "eigenGene value", xlab = "", cex.axis = 0.9, main = colorsCons[mod])
}; dev.off()

#####
#Find the most negatively upregulated genes in each sample
#####
Gene <- rownames(t(Expr[[2]]$data))

#find the genes with the most negative kME, hub genes?
NegNames<-NULL
NegGenesKME = NULL
for (c in 1:length(colorsDt)){
  genes1 <-  Gene[order(geneModuleMembershipDt.DtRef[,c])][1:numGenes]
  data1 <-  geneModuleMembershipDt.DtRef[order(geneModuleMembershipDt.DtRef[,c]),c][1:numGenes]
  genes2 <-  Gene[order(geneModuleMembershipAt.DtRef[,c])][1:numGenes]
  data2 <-  geneModuleMembershipAt.DtRef[order(geneModuleMembershipAt.DtRef[,c]),c][1:numGenes]
  NegGenesKME = cbind(NegGenesKME,genes1,data1,genes2,data2)
  NegNames <- c(NegNames,c(paste("genes1",colorsDt[c],sep=""),paste("data1",colorsDt[c],sep=""),
                           paste("genes2",colorsDt[c],sep=""),paste("data2",colorsDt[c],sep="")))
}; colnames(NegGenesKME)<-NegNames
NegGenesKME

save(NegGenesKME,file="ref_At_vs_Dt_neg_kME_DtRef.Data.R")

#find the genes with the most postive kME, hub genes?
PosNames<-NULL
PosGenesKME = NULL
for (c in 1:length(colorsDt)){
  genes1 <-  Gene[rev(order(geneModuleMembershipAt.DtRef[,c]))][1:numGenes]
  data1 <-  geneModuleMembershipDt.DtRef[rev(order(geneModuleMembershipDt.DtRef[,c])),c][1:numGenes]
  genes2 <-  Gene[rev(order(geneModuleMembershipDt.DtRef[,c]))][1:numGenes]
  data2 <-  geneModuleMembershipAt.DtRef[rev(order(geneModuleMembershipAt.DtRef[,c])),c][1:numGenes]
  PosGenesKME = cbind(PosGenesKME,genes1,data1,genes2,data2)
  PosNames <- c(PosNames,c(paste("genes1",colorsAD1[c],sep=""),paste("data1",colorsAD1[c],sep=""),
                           paste("genes2",colorsAD1[c],sep=""),paste("data2",colorsAD1[c],sep="")))
}; colnames(PosGenesKME)<-PosNames
PosGenesKME

save(PosGenesKME,file="ref_AD1_dom_vs_AD1_wild_pos_kME.Data.R")

######
#Calculate module preservation
#####

modules<-list(mergedColorsAt,mergedColorsDt,mergedColorsA1)
names(modules)<-SetLabels

modPres<-modulePreservation(Expr,modules,nPermutations=10,networkType="signed",referenceNetworks=3)

eigengenes<-list()
for (set in 1:nSets) {
  eigengenes[[set]] = multiSetMEs(Expr, universalColors=modules[[set]], excludeGrey=TRUE)
  for ( ss in 1:nSets ) {
    rownames(eigengenes[[set]][[ss]]$data) = rownames(Expr[[set]]$data)
  }
}

cr <-list()

set.seed(20)

for ( ref in 1:nSets) {
  cr[[ref]] = list()
  for (test in 1:nSets ) {
    printFlush(system.time({
      cr[[ref]][[test]] <- clusterRepro(Centroids= as.matrix(eigengenes[[ref]][[test]]$data),
                                        New.data = as.matrix(Expr[[test]]$data),
                                        Number.of.permutations = 10)
    }));
    collectGarbage()
  }
}

