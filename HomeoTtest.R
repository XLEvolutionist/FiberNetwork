rm(list=ls(all=TRUE))
library(caroline)
library(data.table)
library(ggplot2)
library(VennDiagram)
#Perform student's t-test on homeolog expression across 5 to 20 DPA.

#load in the expression data, normalized
NormData<-load("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/HomeoExprTrimmed.RData")
#check that the data look ok
head(HomeoSort$normalized)
#split the data up by expression
grab<-seq(from=1,to=dim(HomeoSort$normalized)[2], by = 3)
AtHomeo<-HomeoSort$normalized[,grab]
DtHomeo<-HomeoSort$normalized[,grab+1]
#make a vector of DPA labels
DPA<-c(rep(c("10","15","20","05"),3))

#sort the data so that it is in order of DPA
AtHomeo<-AtHomeo[,order(DPA)]
DtHomeo<-DtHomeo[,order(DPA)]
#combine the data agin
HomeoTtest<-cbind(AtHomeo,DtHomeo)

#remove genes that aren't "good" in terms of expression
minNormCount<-0.1
index<-apply(cbind(HomeoTtest), 1, function(x){any(x<minNormCount)})
#now remove those offending genes
HomeoTtest<-HomeoTtest[index == FALSE,]

#a function to test DGE between homeologs
HomeoDGE<-function(range=c(1:3),skip=12,expr, adjust = TRUE) { 
  cols2comp<-range
  cols2comp2<-cols2comp+12
  pval5DPA <- apply(expr,1, function(x)  t.test(x= x[cols2comp], y=x[cols2comp2]))
  #make a table of useful data
  HomeoDGEtable<-NULL
  for ( i in 1:length(pval5DPA)) {
    row<-c(pval5DPA[[i]]$estimate[1],pval5DPA[[i]]$estimate[2],pval5DPA[[i]]$p.value)
    HomeoDGEtable<-rbind(HomeoDGEtable,row)
    
  }#for
  rownames(HomeoDGEtable)<-rownames(expr)
  colnames(HomeoDGEtable)<-c("mean At","mean Dt","pvalue")
  if (adjust == TRUE ) {
    HomeoDGEtable[,"pvalue"]<-p.adjust(HomeoDGEtable[,"pvalue"], method="BH")
    colnames(HomeoDGEtable)[3]<-"adj pvalue"
  }#if
  return(HomeoDGEtable)
}#HomeoDGE

Homeo5DPA<-HomeoDGE(range=c(1:3), skip = 12 , HomeoTtest )
Homeo10DPA<-HomeoDGE(range=c(4:6), skip = 12 , HomeoTtest )
Homeo15DPA<-HomeoDGE(range=c(7:9), skip = 12 , HomeoTtest )
Homeo20DPA<-HomeoDGE(range=c(10:12), skip = 12 , HomeoTtest )

#save(Homeo5DPA,Homeo10DPA,Homeo15DPA,Homeo20DPA, file = "HomeoDGE.RData")
load("/Users/srbyfield/Desktop/WGCNA\ networks/HomeoDGE.RData")

#find the genes homoeologs that are DGE at 5 DPA
Homeo5DPA<-data.frame(Homeo5DPA, row.names=rownames(Homeo5DPA))
Homeo5DPADGE<-subset(Homeo5DPA, adj.pvalue < 0.05)
write.delim(cbind(gene=rownames(Homeo5DPADGE),Homeo5DPADGE), file ="/Users/srbyfield/Desktop/WGCNA\ networks/DGE_Homeologs/Homoe5DPDGE.txt")
#find th genes that are DGE at 10 DPA
Homeo10DPA<-data.frame(Homeo10DPA, row.names=rownames(Homeo10DPA))
Homeo10DPADGE<-subset(Homeo10DPA, adj.pvalue < 0.05)
write.delim(cbind(gene=rownames(Homeo10DPADGE),Homeo10DPADGE), file ="/Users/srbyfield/Desktop/WGCNA\ networks/DGE_Homeologs/Homoe10DPDGE.txt")

#find th genes that are DGE at 15 DPA
Homeo15DPA<-data.frame(Homeo15DPA, row.names=rownames(Homeo15DPA))
Homeo15DPADGE<-subset(Homeo15DPA, adj.pvalue < 0.05)
write.delim(cbind(gene=rownames(Homeo15DPADGE),Homeo15DPADGE), file ="/Users/srbyfield/Desktop/WGCNA\ networks/DGE_Homeologs/Homoe15DPDGE.txt")

#find th genes that are DGE at 15 DPA
Homeo20DPA<-data.frame(Homeo20DPA, row.names=rownames(Homeo20DPA))
Homeo20DPADGE<-subset(Homeo20DPA, adj.pvalue < 0.05)
write.delim(cbind(gene=rownames(Homeo20DPADGE),Homeo20DPADGE), file ="/Users/srbyfield/Desktop/WGCNA\ networks/DGE_Homeologs/Homoe20DPDGE.txt")

#what about up and down-regulated at each stage.
#define a useful function
upDown<-function(x) {
  y<-subset( x, x[,1] > x[,2] )
  z<-subset( x, x[,2] > x[,1] )
  l<-list("AtUp" =y , "DtUp"=z)
  return (l)
}#updown
#grab upregulated At and seperatately upregulated Dt genes
UpDown5<-upDown(Homeo5DPADGE)
UpDown10<-upDown(Homeo10DPADGE)
UpDown15<-upDown(Homeo15DPADGE)
UpDown20<-upDown(Homeo20DPADGE)

#use a binomila test to see if the difference between homoelogs would be expected at random
#at stage 5 DPA
pval_5<-binom.test(dim(UpDown5$DtUp)[1], n=(dim(UpDown5$DtUp)[1]+dim(UpDown5$AtUp)[1]), p=0.5)$p.value 
pval_10<-binom.test(dim(UpDown10$DtUp)[1], n=(dim(UpDown10$DtUp)[1]+dim(UpDown10$AtUp)[1]), p=0.5)$p.value 
pval_15<-binom.test(dim(UpDown15$DtUp)[1], n=(dim(UpDown15$DtUp)[1]+dim(UpDown15$AtUp)[1]), p=0.5)$p.value 
pval_20<-binom.test(dim(UpDown20$DtUp)[1], n=(dim(UpDown20$DtUp)[1]+dim(UpDown20$AtUp)[1]), p=0.5)$p.value 

#make a data.frame for ggplot
dpas<-c(5,5,10,10,15,15,20,20)
homeo<-rep(c("At","Dt"),4)

AtDtBias<-data.frame("dpa"=dpas,
          "Homeo"=homeo,"upreg"=c(dim(UpDown5$AtUp)[1],dim(UpDown5$DtUp)[1],dim(UpDown10$AtUp)[1],
                    dim(UpDown10$DtUp)[1],dim(UpDown15$AtUp)[1],dim(UpDown15$DtUp)[1],
                            dim(UpDown20$AtUp)[1],dim(UpDown20$DtUp)[1]))

summary(AtDtBias)
aggregate(AtDtBias[,3], by = list(AtDtBias$Homeo), sum)
biasBarplot<-ggplot(AtDtBias, aes(y = upreg, x=Homeo)) +
  geom_bar(stat="identity", colour = "royalblue", fill="royalblue")+
  facet_wrap(~dpa, ncol=4)+
  ylab("upregulated genes") +
  theme(axis.title.y = element_text(vjust = 1.8))+
  theme(axis.title.x = element_text(vjust = -0.3))+
  xlab("Homoeolog") +
  theme(axis.title.x = element_text(face="bold", size=20)) +
  theme(axis.text.x = element_text(face="bold", size=20)) +
  theme(axis.title.y = element_text(face="bold", size=20)) +
  theme(axis.text.y = element_text(face="bold", size=20)) +
  theme(panel.background =  element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour="grey50")) +
  theme(strip.text.x = element_text(size=20, face="bold"))
tiff(file="test.tiff")
biasBarplot
dev.off()
