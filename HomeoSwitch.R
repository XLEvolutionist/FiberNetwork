#A script to examine homeolog expression over time.
#Specifically, this assess swithces in homeolog use
#over a time series of cotton fiber development.

#Simon Renny-Byfield, Iowa State University, September 2014

######
##Don't forget to swap out the Maxxa 20 DPA data for Josh's data
######

####
#Functions
####
source("/Users/simonrenny-byfield/scripts/HomeoNorm.R")

#function to find if each row (gene) can distinguish 
#At and Dt well enough

HomeoCheck<-function (x, cut=20 ) {
  #print(x)
  ColSeq<-seq(from=1, to=length(x), by=3)
  for ( i in ColSeq ) {
    AtDt<-as.numeric((x[i] + x[i+1]))
    N<-as.numeric(x[i+2])
    #print(AtDt)
    #print(N)
    if ( is.finite(AtDt) & AtDt > 0) {
      N<-x[i+2]
      perc<-(AtDt/(AtDt+N))*100
      #print(perc)
      if ( perc < cut ) {
      return(FALSE)
      #stop("")
      }#if
    }#if
    else {
      return(FALSE)
    }#else
  }#for
  return(TRUE)
}#HomeoCheck


#Load in count data
HomeoData<-read.table("/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/AD1_domesticated_homeo_counts.txt", header = T , sep = "\t")
#Normalize the data
HomeoDataNorm<-HomeoNorm(HomeoData)
LibSizes<-HomeoDataNorm$sizes

HomeoDataNorm<-HomeoDataNorm$normalized

#screen the genes to make sure those with insuffiient expression, or
#too many Ns, are not included in the analysis
minRPM<-0.2
#rowmeans have to be over minRPM
HomeoDataNorm<-HomeoDataNorm[rowMeans(HomeoDataNorm)> minRPM,]

#now where the At and Dt counts are at least 20% of the total counts
#At + Dt + N
keep<-apply(HomeoDataNorm,1,function(x) HomeoCheck(x))
HomeoTrim<-HomeoDataNorm[keep,]
#sort the columns of the matrix in order to get Maxxa together
HomeoSort<-HomeoTrim[,sort(colnames(HomeoTrim))]

save(HomeoSort, file="/Users/simonrenny-byfield/cotton/diploid_domestication/FiberTranscriptomeProject/homeoData/HomeoExprTrimmed.RData")
#set up a string contaning sample info
samples<-rep(c("10A","10D","10N","15A","15D","15N","20A","20D","20N","5A","5D","5N"),3)

HomeoMeans<-aggregate(t(HomeoSort),by=list(samples), FUN=mean)
rownames(HomeoMeans)<-HomeoMeans[,1]
HomeoMeans<-HomeoMeans[,-1]

HomeoMeans<-t(HomeoMeans[c(10:12,1:9),])

#sepearate At and Dt counts
grab<-c(1,4,7,10)
AtMeans<-HomeoMeans[,grab]
DtMeans<-HomeoMeans[,grab+1]
totalAtDt<-DtMeans+AtMeans

fractionDt<-DtMeans/totalAtDt

#trim the data so that only those with min values of less than 30, 
#and max values of greater than 70

#find the min and max Dt percentage
min<-apply(fractionDt, 1, min) 
max<-apply(fractionDt, 1, max) 

swithcData<-fractionDt[(min <= 0.3) & (max >= 0.7),]


plot(rep(seq(from=5,to=20,by=5),dim(HomeoMeans)[1]),t(fractionDt), type="n", xlab="development DPA", ylab = "fraction Dt contribution", cex.lab = 1.4,cex.axis=1.4 )
lines(rep(seq(from=5,to=20,by=5),dim(HomeoMeans)[1]),t(fractionDt), lwd= 0.015,col="black")


plot(rep(seq(from=5,to=20,by=5),dim(HomeoMeans)[1]),t(fractionDt), type="n", xlab="development DPA", ylab = "fraction Dt contribution", cex.lab = 1.4,cex.axis=1.4  )
lines(rep(seq(from=5,to=20,by=5),dim(swithcData)[1]),t(swithcData), lwd= 0.15, col = "red")

#use clustering to cluster genes that switch
testDist<-dist(swithcData,method = "euclidean")
clust<-hclust(testDist)
#cut the tree at hiehgt 0.8
groups<-cutree(clust, h=.8)
#try those ALL genes
AlltestDist<-dist(fractionDt,method = "euclidean")
Allclust<-hclust(AlltestDist)
#cut the tree at hiehgt 0.8
groups<-cutree(Allclust, h=.8)

pdf("clusterHomeo.pdf", height = 10, width = 15)
par(mfrow=c(2,2))
for ( i in unique(groups)) {
 print(i)
 plot(rep(seq(from=5,to=20,by=5),dim(HomeoMeans)[1]),t(fractionDt), type="n", xlab="development DPA", ylab = "fraction Dt contribution", cex.lab = 1.6,cex.axis=1.6  )
 aveLine<-rowMeans(t(swithcData[names(groups[groups == i]),]))
 lines(rep(seq(from=5,to=20,by=5),dim(swithcData[names(groups[groups == i]),])[1]),t(swithcData[names(groups[groups == i]),]), lwd= 0.05, col = "red")
 lines(seq(from=5,to=20,by=5),aveLine, lwd= 4, col = "black", lty = "dashed")
 text(7.5,0.9,paste("N=",dim(swithcData[names(groups[groups == i]),])[1]), cex = 1.5)
}
dev.off()

