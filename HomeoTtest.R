rm(list=ls(all=TRUE))

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

save(Homeo5DPA,Homeo10DPA,Homeo15DPA,Homeo20DPA, file = "HomeoDGE.RData")

