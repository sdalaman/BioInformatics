###############  CLUSTERING  

numOfGenes = 500


# Looking at Total within-cluster sum of squares for different k
# to select best k
bestKmeans=function(data,maxk,heading){ 
  set.seed(10)
  maxk=maxk
  err=numeric(maxk)
  errchn=numeric(maxk)
  errchn2=numeric(maxk)
  last=0
  last2=0
  for (k in 1:maxk){
    km.out=kmeans(data,k,nstart=150)
    err[k]=km.out$tot.withinss
    errchn[k]=km.out$tot.withinss - last
    errchn2[k]=errchn[k] - last2
    last=km.out$tot.withinss
    last2=errchn[k]
    cat("# of centers : " ,k, "\n")
  }
  matplot(1:maxk,cbind(err,errchn,errchn2),pch=19,col=c("red","blue","black"),
          type="b",ylab="Mean Squared Error",xlab="# of clusters",main=heading)
  legend("topright",legend=c("Err","Errchn","Errchn2"),pch=19,col=c("red","blue","black"))
}

bestKmeans(dataTopAll[1:numOfGenes,],10,"Combined")
bestKmeans(dataTopDiabetes[1:numOfGenes,],10,"Diabetes")
bestKmeans(dataTopLeukemia[1:numOfGenes,],10,"Leukemia")
bestKmeans(dataTopParkinson[1:numOfGenes,],10,"Parkinson")
bestKmeans(dataTopAsthma[1:numOfGenes,],10,"Asthma")
bestKmeans(dataTopPancreatic[1:numOfGenes,],10,"Pacreatic Cancer")
bestKmeans(dataTopLung[1:numOfGenes,],10,"Lung Cancer")
bestKmeans(dataTopBreast[1:numOfGenes,],10,"Breast Cancer")

nclusAll = 3
nclusDiabetes = 3
nclusLeukemia = 3
nclusParkinson = 3
nclusAsthma = 3
nclusPancreatic = 3
nclusLung = 3
nclusBreast = 3


# K-Means Clustering with  clusters
numOfGenes = 500
kMeansModelAll <- kmeans(dataTopAll[1:numOfGenes,], nstart=150,nclusAll)
kMeansModelDiabetes <- kmeans(dataTopDiabetes[1:numOfGenes,], nstart=150,nclusDiabetes)
kMeansModelLeukemia <- kmeans(dataTopLeukemia[1:numOfGenes,], nstart=150,nclusLeukemia)
kMeansModelParkinson <- kmeans(dataTopParkinson[1:numOfGenes,], nstart=150,nclusParkinson)
kMeansModelAsthma <- kmeans(dataTopAsthma[1:numOfGenes,], nstart=150,nclusAsthma)
kMeansModelPancreatic <- kmeans(dataTopPancreatic[1:numOfGenes,], nstart=150,nclusPancreatic)
kMeansModelLung <- kmeans(dataTopLung[1:numOfGenes,], nstart=150,nclusLung)
kMeansModelBreast <- kmeans(dataTopBreast[1:numOfGenes,], nstart=150,nclusBreast)

# Cluster Plot against 1st 2 principal components
# vary parameters for most readable graph
library(cluster) 

clusplot(dataTopAll[1:numOfGenes,], kMeansModelAll$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Combined Samples)")
clusplot(dataTopDiabetes[1:numOfGenes,], kMeansModelDiabetes$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Diabetes)")
clusplot(dataTopLeukemia[1:numOfGenes,], kMeansModelLeukemia$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Leukemia)")
clusplot(dataTopParkinson[1:numOfGenes,], kMeansModelParkinson$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Parkinson)")
clusplot(dataTopAsthma[1:numOfGenes,], kMeansModelAsthma$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Asthma)")
clusplot(dataTopPancreatic[1:numOfGenes,], kMeansModelPancreatic$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Pancreatic Cancer)")
clusplot(dataTopLung[1:numOfGenes,], kMeansModelLung$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Lung Cancer)")
clusplot(dataTopBreast[1:numOfGenes,], kMeansModelBreast$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0,main="PCA1 VS PC2 Cluster Plot(Breast Cancer)")


# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(dataTopAll[1:numOfGenes,], kMeansModelAll$cluster,main="Centroid Plot(Combined Samples)")
plotcluster(dataTopDiabetes[1:numOfGenes,], kMeansModelDiabetes$cluster,main="Centroid Plot(Diabetes)")
plotcluster(dataTopLeukemia[1:numOfGenes,], kMeansModelLeukemia$cluster,main="Centroid Plot(Leukemia)")
plotcluster(dataTopParkinson[1:numOfGenes,], kMeansModelParkinson$cluster,main="Centroid Plot(Parkinson)")
plotcluster(dataTopAsthma[1:numOfGenes,], kMeansModelAsthma$cluster,main="Centroid Plot(Asthma)")
plotcluster(dataTopPancreatic[1:numOfGenes,], kMeansModelPancreatic$cluster,main="Centroid Plot(Pancreatic Cancer)")
plotcluster(dataTopLung[1:numOfGenes,], kMeansModelLung$cluster,main="Centroid Plot(Lung Cancer)")
plotcluster(dataTopBreast[1:numOfGenes,], kMeansModelBreast$cluster,main="Centroid Plot(Breast Cancer)")

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

mdsPlot=function(data,header){
  dist <- dist(data) # euclidean distances between the rows
  cmdModel <- cmdscale(dist,eig=TRUE, k=2) # k is the number of dim
  # plot solution 
  xc <- cmdModel$points[,1]
  yc <- cmdModel$points[,2]
  plot(xc, yc, xlab="Coordinate 1", ylab="Coordinate 2", 
     main=header,	type="n")
  text(xc, yc, labels = row.names(data), cex=.7)
}

mdsPlot(dataTopAll[1:numOfGenes,],"Metric	MDS(Combined Sample)")
mdsPlot(dataTopDiabetes[1:numOfGenes,],"Metric	MDS(Diabetes Sample)")
mdsPlot(dataTopLeukemia[1:numOfGenes,],"Metric	MDS(Leukemia Sample)")
mdsPlot(dataTopParkinson[1:numOfGenes,],"Metric	MDS(Parkinson Sample)")
mdsPlot(dataTopAsthma[1:numOfGenes,],"Metric	MDS(Asthma Sample)")
mdsPlot(dataTopPancreatic[1:numOfGenes,],"Metric	MDS(Pancreatic Cancer Sample)")
mdsPlot(dataTopLung[1:numOfGenes,],"Metric	MDS(Lung CancerSample)")
mdsPlot(dataTopBreast[1:numOfGenes,],"Metric	MDS(Breast Cancer Sample)")

hist(varGenesAll[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Combined Samples)",
     xlab="Gene Variance for Combined Samples")
hist(varGenesDiabetes[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Diabetes)",
     xlab="Gene Variance for Diabetes Samples")
hist(varGenesLeukemia[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Leukemia)",
     xlab="Gene Variance for Leukemia Samples")
hist(varGenesParkinson[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Parkinson)",
     xlab="Gene Variance for Parkinson Samples")
hist(varGenesAsthma[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Asthma)",
     xlab="Gene Variance for Asthma Samples")
hist(varGenesPancreatic[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Pancreatic Cancer)",
     xlab="Gene Variance for Pancreatic Samples")
hist(varGenesLung[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Lung Camcer)",
     xlab="Gene Variance for Lung Samples")
hist(varGenesBreast[[1]], breaks=100,
     main="Historgram of gene-wise distribution (Breast Cancer)",
     xlab="Gene Variance for Breast Samples")


## Plot the expression profiles of the two genes with highest variance
gn1 <- geneLstTopAll[[1]][1]
gn2 <- geneLstTopAll[[1]][2]
print(list(gn1,gn2))

g1 <- as.vector(as.matrix(allData[gn1,]))
g2 <- as.vector(as.matrix(allData[gn2,]))
plot(g1,g2,
     col=all.group.colors,
     type='n',
     panel.first=grid(col='black'), 
     main="2 genes with the highest variance", 
     xlab=paste('gene', gn1), 
     ylab=paste('gene', gn2))
text(g1, g2,labels=all.group.abbrev,col=all.group.colors,pch=0.5)
legend('topright',col=all.group.colors, 
       legend=names(all.group.colors),pch=0.5,cex=0.5,bg='white',inset=0.01,x.intersp=2,xjust=0,yjust=0)

library(knitr)
tbl <- NULL
ln <- c("ALL","Diabetes","Leukemia","Parkinson","Asthma","Pancreatic","Lung","Breast")
tbl <- rbind(ln)
for(i in 1:10){
  ln <- c(geneLstTopAll[[1]][i],geneLstTopDiabetes[[1]][i],
          geneLstTopLeukemia[[1]][i],geneLstTopParkinson[[1]][i],
          geneLstTopAsthma[[1]][i],geneLstTopPancreatic[[1]][i],
          geneLstTopLung[[1]][i],geneLstTopBreast[[1]][i])
  tbl <- rbind(tbl,ln)
}
rownames(tbl)[1] <- "Disease"
for(i in 2:11)
  rownames(tbl)[i] <- paste("Gene ",as.character(i-1))
print(tbl)

sink('tbl.txt')
print(tbl,right=F)
sink()

library(stats) 
numOfGenes = 500
## Perform the PCA transformation
allData.prcomp <- prcomp(t(dataTopAll[1:numOfGenes,]),cor=TRUE)
diabetesData.prcomp <- prcomp(t(dataTopDiabetes[1:numOfGenes,]),cor=TRUE)
leukemiaData.prcomp <- prcomp(t(dataTopLeukemia[1:numOfGenes,]),cor=TRUE)
asthmaData.prcomp <- prcomp(t(dataTopAsthma[1:numOfGenes,]),cor=TRUE)
pancreaticData.prcomp <- prcomp(t(dataTopPancreatic[1:numOfGenes,]),cor=TRUE)
parkinsonData.prcomp <- prcomp(t(dataTopParkinson[1:numOfGenes,]),cor=TRUE)
lungData.prcomp <- prcomp(t(dataTopLung[1:numOfGenes,]),cor=TRUE)
breastData.prcomp <- prcomp(t(dataTopBreast[1:numOfGenes,]),cor=TRUE)

## Analyze the content of the prcomp result: 
## the result of the method prcomp() is an object 
## belonging to the class "prcomp"

#class(allData.prcomp) 
#attributes(allData.prcomp)
#plot(allData.prcomp, main='Variance  per component', xlab='Component')
#all.sd.per.pc <- allData.prcomp$sdev
#all.var.per.pc <- all.sd.per.pc^2

## Display the percentage of total variance explained by each 
#all.sd.per.pc.percent <- all.sd.per.pc/sum(all.sd.per.pc)
#all.var.per.pc.percent <- all.var.per.pc/sum(all.var.per.pc)

#barplot(all.var.per.pc.percent[1:10], main='Percent of variance  per component'
#        , xlab='Component', ylab='Percent variance', col='#BBDDFF')


plotPCA=function(prcomp,data,sampleType){
  if(trimws(tolower(sampleType)) == "combined"){
    sample.subtypes <- as.vector(projetab$disease)
    sample.labels <- all.group.abbrev[sample.subtypes]
    names(sample.labels) <- names(data)
    sample.colors <- all.group.colors[as.vector(projetab$disease)]
    names(sample.colors) <- names(data)}
  else {
    sample.subtypes <- as.vector(projetab[c(projetab$disease == trimws(tolower(sampleType))),])
    sample.labels <- all.group.abbrev[sample.subtypes]
    names(sample.labels) <- names(data)
    sample.colors <- all.group.colors[as.vector(projetab[c(projetab$disease ==trimws(tolower(sampleType))),])]
    names(sample.colors) <- names(data)
  }
  
  biplot(prcomp,var.axes=FALSE,
       panel.first=grid(col='black'), 
       main=paste('PCA Plot ',sampleType,
                  ncol(data), 'samples *', nrow(data), 'genes', sep=' '), 
       xlab='First component', ylab='Second component')
  plot(prcomp$x[,1:2],
     col=sample.colors,
     type='n',
     panel.first=grid(col='black'), 
     main=paste('PCA Plot ',sampleType,
                ncol(data), 'samples *', nrow(data), 'genes', sep=' '), 
     xlab='PC1', ylab='PC2')
  text(prcomp$x[,1:2],labels=sample.labels,col=sample.colors,pch=0.5,font=4)
#legend('bottomleft',col=all.group.colors, 
#       legend=names(all.group.colors),pch=1,cex=0.7,bg='white',bty='o')
## Plot components PC2 and PC3
  plot(prcomp$x[,2:3],
     col=sample.colors,
     type='n',
     panel.first=grid(col='black'), 
     main=paste('PCA Plot ',sampleType,
                ncol(data), 'samples *', nrow(data), 'genes', sep=' '), 
     xlab='PC2', ylab='PC3')
  text(prcomp$x[,2:3],labels=sample.labels,col=sample.colors,pch=0.5,font=4)     
#  legend('bottomleft',col=all.group.colors, 
#       legend=names(all.group.colors),pch=1,cex=0.7,bg='white',bty='o')
}

plotPCA(allData.prcomp,dataTopAll[1:numOfGenes,],'Combined ')
plotPCA(diabetesData.prcomp,dataTopDiabetes[1:numOfGenes,],'Diabetes ')
plotPCA(leukemiaData.prcomp,dataTopLeukemia[1:numOfGenes,],'Leukemia ')
plotPCA(parkinsonData.prcomp,dataTopParkinson[1:numOfGenes,],'Parkinson ')
plotPCA(pancreaticData.prcomp,dataTopPancreatic[1:numOfGenes,],'Pancreatic')
plotPCA(asthmaData.prcomp,dataTopAsthma[1:numOfGenes,],'Asthma ')
plotPCA(lungData.prcomp,dataTopLung[1:numOfGenes,],'Lung')
plotPCA(breastData.prcomp,dataTopBreast[1:numOfGenes,],'Breast')


numOfGenes=500
pcaFitAll <- princomp(dataTopAll[1:numOfGenes,], cor=TRUE)
summary(pcaFitAll) # print variance accounted for 
loadings(pcaFitAll,cutoff = 0.1) # pc loadings 
plot(pcaFitAll,type="lines",main = "PCA-Combined") # scree plot 
#fit$scores # the principal components
biplot(pcaFitAll,main="PCA Biplot-Combined")
lstPcaAll <- pcaFitAll$scores
importanceOrder=order(-lstPcaAll[,1])
pcaNamesAll=rownames(lstPcaAll)[importanceOrder][1:50]
cat("List of First 50 features of PCA1")
pcaNamesAll

pcaFitDiabetes <- princomp(dataTopDiabetes[1:numOfGenes,], cor=TRUE)
summary(pcaFitDiabetes) # print variance accounted for 
loadings(pcaFitDiabetes,cutoff = 0.1) # pc loadings 
plot(pcaFitDiabetes,type="lines",main = "PCA-Diabetes") # screen plot 
#fit$scores # the principal components
biplot(pcaFitDiabetes,main="PCA Biplot-Diabetes")
lstPcaDiabetes <- pcaFitDiabetes$scores
importanceOrder=order(-lstPcaDiabetes[,1])
pcaNamesDiabetes=rownames(lstPcaDiabetes)[importanceOrder][1:50]
cat("List of First 50 features of PCA1")
pcaNamesDiabetes

pcaFitLung <- princomp(dataTopLung[1:numOfGenes,], cor=TRUE)
summary(pcaFitLung) # print variance accounted for 
loadings(pcaFitLung,cutoff = 0.1) # pc loadings 
plot(pcaFitLung,type="lines",main = "PCA-Lung") # screen plot 
#fit$scores # the principal components
biplot(pcaFitLung,main="PCA Biplot-Lung")
lstPcaLung <- pcaFitLung$scores
importanceOrder=order(-lstPcaLung[,1])
pcaNamesLung=rownames(lstPcaLung)[importanceOrder][1:50]
cat("List of First 50 features of PCA1")
pcaNamesLung



#x <- (rownames(allData) == "")
#count(x,value= TRUE)
#count(x, value=FALSE)

pca.out=prcomp(dataTopAll[1:numOfGenes,], scale=TRUE)
biplot(pca.out, scale=0)
pca_out_df <- as.data.frame(pca.out$x)
library(scales)
ramp <- colorRamp(c("yellow", "blue"))
colours_by_mean <- rgb( 
  ramp( as.vector(rescale(rowMeans(pca_out_df),c(0,1)))), 
  max = 255 )
plot(PC1~PC2, data=pca_out_df,
     main= "PCA Mean Dist",
     cex = .1, lty = "solid", col=colours_by_mean)
text(PC1~PC2, data=pca_out_df, 
     labels=rownames(allData),
     cex=.8, col=colours_by_mean)
lst <- pca.out$rotation
importanceOrder=order(-lst[,"PC1"])
names=rownames(lst)[importanceOrder][1:50]
cat("List of First 50 features of PCA1")
names
names=rownames(lst)[order(-lst[,"PC2"])][1:50]
cat("List of First 50 features of PCA2")
names


# Random Forests
library(randomForest)
#library("cluster")

to.dendrogram <- function(dfrep,rownum=1,height.increment=0.1){
  
  if(dfrep[rownum,'status'] == -1){
    rval <- list()
    
    attr(rval,"members") <- 1
    attr(rval,"height") <- 0.0
    attr(rval,"label") <- dfrep[rownum,'prediction']
    attr(rval,"leaf") <- TRUE
    
  }else{##note the change "to.dendrogram" and not "to.dendogram"
    left <- to.dendrogram(dfrep,dfrep[rownum,'left daughter'],height.increment)
    right <- to.dendrogram(dfrep,dfrep[rownum,'right daughter'],height.increment)
    rval <- list(left,right)
    
    attr(rval,"members") <- attr(left,"members") + attr(right,"members")
    attr(rval,"height") <- max(attr(left,"height"),attr(right,"height")) + height.increment
    attr(rval,"leaf") <- FALSE
    attr(rval,"edgetext") <- dfrep[rownum,'split var']
  }
  
  class(rval) <- "dendrogram"
  
  return(rval)
}

rfCalc=function(data,sample.type){
  rfModel=randomForest(t(data),importance=TRUE,keep.forest=TRUE)
  bestKmeans(rfModel$proximity,min(9,nrow(rfModel$proximity)-1),sample.type)
  lst <- rfModel$importance
  importanceOrder=order(-lst[,1])
  topGenesRF=rownames(lst)[importanceOrder][1:5000]
  varImpPlot(rfModel,type=2,sort=TRUE, main=paste("Gene Importance Order for ",sample.type))
  return(list(rfModel,topGenesRF)) 
}
  
rfCalcClust=function(data,rfModelT,nClust,sample.type){
  E_rf <- kmeans(rfModelT$proximity, 3, iter.max = 100, nstart = 10)
  rfModel <- randomForest(t(data),as.factor(E_rf$cluster),ntree = 500)
  lst <- rfModel$importance
  importanceOrder=order(-lst[,1])
  topGenesRF=rownames(lst)[importanceOrder][1:5000]
  #cat("List of First 5000 features of Random Forest for",sampleType)
  #topGenesRFAll
  #cat("List of common genes between Variance ordering and Random Forest importance order")
  #intersect(topGenesRFAll[1:500],geneLstTopAll[[1]][1:500])
  varImpPlot(rfModel,type=2,sort=TRUE, main=paste("Gene Importance Order for ",sample.type))
  tree <- getTree(rfModel,1,labelVar=TRUE)
  dnd <- to.dendrogram(tree)
  str(dnd)
  plot(dnd,center=TRUE,leaflab='none',edgePar=list(t.cex=1,p.col=NA,p.lty=0),
       main=paste("Dendogram of Gene Importance Order for ",sample.type))
  return(list(rfModel,topGenesRF)) 
}

ret=rfCalc(allData,"Combined Samples")
topGenesRFAll=ret[[2]]
rfModelAll=ret[[1]] 

ret=rfCalc(logDiabetes,"Diabetes Samples")
topGenesRFDiabetes=ret[[2]]
rfModelDiabetes=ret[[1]] 

ret=rfCalc(logLeukemia,"Leukemia Samples")
topGenesRFLeukemia=ret[[2]]
rfModelLeukemia=ret[[1]] 

ret=rfCalc(logAsthma,"Asthma Samples")
topGenesRFAsthma=ret[[2]]
rfModelAsthma=ret[[1]] 

ret=rfCalc(logParkinson,"Parkinson Samples")
topGenesRFParkinson=ret[[2]]
rfModelParkinson=ret[[1]] 

ret=rfCalc(logPancreatic,"Pancreatic Cancer Samples")
topGenesRFPancreatic=ret[[2]]
rfModelPancreatic=ret[[1]] 

ret=rfCalc(logLung,"Lung Cancer Samples")
topGenesRFLung=ret[[2]]
rfModelLung=ret[[1]] 

ret=rfCalc(logBreast,"Breast Cancer Samples")
topGenesRFBreast=ret[[2]]
rfModeBreast=ret[[1]] 









