##
#A script to pwrform network analysis using WGCNA
#Simon Renny-Byfield, Iowa State University, Aug 2014
##
library(edgeR)
library(vegan)
library(ggplot2)
library(WGCNA)
library(plotrix)
library(RColorBrewer)

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

#output the params in a file

#set the wd
setwd(paste("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/RPM_sft_", power , sep=""))

fileConn<-file("params.txt")
writeLines(c(paste("minModuleSize=",minModuleSize,sep=""),c(paste("mergeCutHeight=",mergeCutHeight,sep="")),
             paste("power=",power,sep=""),paste("TOMType=",TOMType,sep=""),paste("minNormCount=",minNormCount,sep="")),fileConn)
close(fileConn)

#load in the sample data
sampleInfo<-read.csv("SampleInfo.csv")
#load in the normalized data
load(file="FiberTranscriptome_RPM.R")
#eliminate those genes with which we do not have good expression data for
#remove gene if it falls below a count of at least 0.75 RPMK in all samples

index<-apply(exp.Dat, 1, function(x){any(x<minNormCount)})
#remove those genes falling below the threshold
exp.Dat<-exp.Dat[!index,]

exp.Dat<-as.data.frame(t(exp.Dat))
#exp.Dat<-log(exp.Dat+1)

goodGenes = goodSamplesGenes(exp.Dat, verbose = 3);
goodGenes$allOK
#remove the genes if WGCNA thinks they are bad

traitRows = match(sampleInfo$sample, rownames(exp.Dat));
A1Traits = sampleInfo[traitRows, ];
#A1Traits<-na.omit(A1Traits)
rownames(A1Traits)<-A1Traits[,1]
A1Traits<-A1Traits[,-1]
A1Traits[,2]<-as.character(A1Traits[,2])
collectGarbage();

cex1=1.3

# Re-cluster samples
sampleTree2 = flashClust(dist(exp.Dat, method = "euclidean"))
# Convert traits to a color representation: white means low, red means high, grey means missing entry
#turn the factor DPA into a colors
colNums<-as.numeric(factor(A1Traits$DPA))

traitColors = labels2colors(A1Traits)
# Plot the sample dendrogram and the colors underneath.

pdf(file="Sample_Dendrogram.pdf", height = 6 , width = 7)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(A1Traits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
par(mfrow=c(1,1), mgp= c(3, 1, 0))
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(exp.Dat, powerVector = powers, verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file="softThresholding.pdf", height = 6 , width = 7)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
    #abline(h=power/10,col="red")
dev.off()


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

####
# Make a multiple set dataset, so we can calcualte consensus modules
####
nSets = 4
multiExpr<-vector(mode="list",length=nSets)
#load in each data set [1] A1 dom , [2] A1 wild, [3] AD1 dom and [4] AD1 wild

multiExpr[[1]]<-list(data=as.data.frame(exp.Dat[1:12,]))
multiExpr[[2]]<-list(data=as.data.frame(exp.Dat[13:24,]))
multiExpr[[3]]<-list(data=as.data.frame(exp.Dat[25:36,]))
multiExpr[[4]]<-list(data=as.data.frame(exp.Dat[37:48,]))
#save the data
save(multiExpr, file = "multiExprData.R")

netCon = blockwiseConsensusModules(multiExpr, power = power,
                       TOMType = TOMType, minModuleSize = minModuleSize,
                       reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                       verbose = 3,
                       maxBlockSize=15000)

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
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
setLabels<-c("A1 dom" , "A1 wild", "AD1 dom", "AD1 wild")
pdf("eigenGeneComparison.pdf", width = 10, height = 10)
plotEigengeneNetworks(consMEsC, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                       xLabelsAngle = 90)
dev.off()
#####
# Calculate modules for A1 dom
#####

A1.dom.net = blockwiseModules(multiExpr[[1]]$data, power = power,
                       TOMType = TOMType, minModuleSize = minModuleSize,
                       reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                       verbose = 3,
                       maxBlockSize=15000)

mergedColorsA1.dom = labels2colors(A1.dom.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("A1_dom_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(A1.dom.net$dendrograms[[1]], mergedColorsA1.dom[A1.dom.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="A1.dom Cluster Dendrogram")
dev.off()
table(A1.dom.net$colors)

#save the A1 dom modules and plot
save(A1.dom.net, file="A1.dom.modules.data.R")

#print a pdf of the cladogram and module assignments

#####
# Calculate modules for A1 wild
#####

A1.wild.net = blockwiseModules(multiExpr[[2]]$data, power = power,
                              TOMType = TOMType, minModuleSize = minModuleSize,
                              reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                              verbose = 3,
                              maxBlockSize=15000)

mergedColorsA1.wild = labels2colors(A1.wild.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("A1_wild_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(A1.wild.net$dendrograms[[1]], mergedColorsA1.wild[A1.wild.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="A1.wild Cluster Dendrogram")
dev.off()

table(A1.wild.net$colors)

#save the A1 wild modules and plot
save(A1.wild.net, file="A1.wild.modules.data.R")

#####
# Calculate modules for AD1 dom
#####

AD1.dom.net = blockwiseModules(multiExpr[[3]]$data, power = power,
                               TOMType = TOMType, minModuleSize = minModuleSize,
                               reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                               verbose = 3,
                               maxBlockSize=15000)

mergedColorsAD1.dom = labels2colors(AD1.dom.net$colors)
# Plot the dendrogram and the module colors underneath

pdf("AD1_dom_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(AD1.dom.net$dendrograms[[1]], mergedColorsAD1.dom[AD1.dom.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="AD1.dom Cluster Dendrogram")
dev.off()
table(AD1.dom.net$colors)

#save the AD1 dom modules and plot
save(AD1.dom.net, file="AD1.dom.modules.data.R")

#####
# Calculate modules for AD1 wild
#####

AD1.wild.net = blockwiseModules(multiExpr[[4]]$data, power = power,
                               TOMType = TOMType, minModuleSize = minModuleSize,
                               reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                               verbose = 3,
                               maxBlockSize=15000)

mergedColorsAD1.wild = labels2colors(AD1.wild.net$colors)
# Plot the dendrogram and the module colors underneath
#print a pdf of the cladogram and module assignments
pdf("AD1_wild_cladogram.pdf", height = 6 , width = 7)
plAD1Wild<-plotDendroAndColors(AD1.wild.net$dendrograms[[1]], mergedColorsAD1.wild[AD1.wild.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="AD1.wild Cluster Dendrogram")
dev.off()
table(AD1.wild.net$colors)

save(AD1.wild.net, file ="AD1.wild.modules.data.R")

#####
# Calculate modules for AD1
#####
#Make an expression data set
ADExpr<-rbind(multiExpr[[3]]$data,multiExpr[[4]]$data)

AD1.net = blockwiseModules(ADExpr, power = power,
                                TOMType = TOMType, minModuleSize = minModuleSize,
                                reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                                numericLabels = TRUE, pamRespectsDendro = FALSE,
                                saveTOMs = TRUE,
                                saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                                verbose = 3,
                                maxBlockSize=15000)

mergedColorsAD1 = labels2colors(AD1.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("AD1_cladogram.pdf", height = 6 , width = 7)
plAD1<-plotDendroAndColors(AD1.net$dendrograms[[1]], mergedColorsAD1[AD1.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="AD1.wild Cluster Dendrogram")

dev.off()

table(AD1.net$colors)

save(AD1.net, file ="AD1.modules.data.R")

#####
# Calculate modules for A1
#####

A1Expr<-rbind(multiExpr[[1]]$data,multiExpr[[2]]$data)

A1.net = blockwiseModules(A1Expr, power = power,
                           TOMType = TOMType, minModuleSize = minModuleSize,
                           reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/WGCNA_analysis/Modules.ALL.Samples",
                           verbose = 3,
                           maxBlockSize=15000)

mergedColorsA1 = labels2colors(A1.net$colors)
# Plot the dendrogram and the module colors underneath
pdf("A1_cladogram.pdf", height = 6 , width = 7)
plotDendroAndColors(A1.net$dendrograms[[1]], mergedColorsA1[A1.net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main="AD1.wild Cluster Dendrogram")
dev.off()

table(A1.net$colors)

save(A1.net, file ="A1.modules.data.R")

