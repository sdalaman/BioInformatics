# Introduction to Microarray data analysis
library('pheatmap')
library('golubEsets')
library('GEOquery')
library('affy')
library('matrixStats')
library('dendextend')
######################################################################################
#Download GPL file nad annotations file, put it in the current directory, and load it:
######################################################################################
gpl97 <- getGEO(filename='GPL97.annot')
Meta(gpl97)$title
colnames(Table(gpl97))
Table(gpl97)[1:10,1:4]
Table(gpl97)[1:10,c("ID","Gene title")]
Table(gpl97)[1:10,c("ID","Gene title","Gene symbol","Gene ID","GenBank Accession")]
IDs <- attr(dataTable(gpl97), "table")[,c("ID", "Gene symbol")] #create a table of ID and gene symbols

anno <- Table(gpl97)
projetab <- read.table('projectData.tab', sep='\t', head=TRUE, row=1)


all.group.abbrev <- c(
  'diabetes'='ds',
  'leukemia'='la',
  'parkinson'='pn',
  'asthma'='aa',
  'pancreatic'='pc',
  'lung'='lg',
  'breast'='bt'
)

all.group.colors <- c(
  'diabetes'='cyan',
  'leukemia'='black',
  'parkinson'='darkgray',
  'asthma'='green',
  'pancreatic'='orange',
  'lung'='violet',
  'breast'='blue'
)

## Unknown Gene Symbol assigned as Gene ID
lst <- which(IDs[,2] == "")
IDs[,1] <- as.character(IDs[,1] )
IDs[,2] <- as.character(IDs[,2] )
IDs[lst,2] <- IDs[lst,1]
IDs[,1] <- as.factor(IDs[,1] )
IDs[,2] <- as.factor(IDs[,2] )

#############################################
#Read the table data
#############################################
diabetes <- read.table("GSE9006-GPL97_series_matrix.txt",skip = 74, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
diabetes<-diabetes[-22646,]  
diabetes <- diabetes[,c("GSM254184","GSM254185", "GSM254187", "GSM254189", "GSM254190",
  "GSM254186", "GSM254188", "GSM254194", "GSM254195", "GSM254196")]
dim(diabetes)
names(diabetes)
round(apply(diabetes, 2, summary))

leukemia <- read.table("GSE22529-GPL97_series_matrix.txt",skip = 59, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
leukemia<-leukemia[-22646,]  
leukemia<-leukemia[,c("GSM559433","GSM559434","GSM559436","GSM559437","GSM559438",
                      "GSM559440","GSM559441","GSM559442","GSM559444","GSM559445")]
dim(leukemia)
names(leukemia)
round(apply(leukemia, 2, summary))

parkinson <- read.table("GSE8397-GPL97_series_matrix.txt",skip = 58, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
parkinson<-parkinson[-22646,]  
parkinson<-parkinson[,c("GSM208700","GSM208701","GSM208702","GSM208703","GSM208704",
                      "GSM208705","GSM208706","GSM208707","GSM208708","GSM208709")]
dim(parkinson)
names(parkinson)
round(apply(parkinson, 2, summary))

asthma <- read.table("GSE473-GPL97_series_matrix.txt",skip = 64, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
asthma<-asthma[-22646,]  
asthma<-asthma[,c("GSM3922","GSM3924","GSM3926","GSM3928","GSM3930",
                        "GSM3932","GSM3934","GSM3936","GSM3938","GSM3940")]
dim(asthma)
names(asthma)
round(apply(asthma, 2, summary))

pancreatic <- read.table("GSE43288-GPL97_series_matrix.txt",skip = 76, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
pancreatic<-pancreatic[-22646,]  
pancreatic<-pancreatic[,c("GSM1060008","GSM1060009","GSM1060010","GSM1060011","GSM1060019",
                  "GSM1060020","GSM1060021","GSM1060022","GSM1060023","GSM1060024")]
dim(pancreatic)
names(pancreatic)
round(apply(pancreatic, 2, summary))

lung <- read.table("GSE31908-GPL97_series_matrix.txt",skip = 86, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
lung<-lung[-22646,]  
lung<-lung[,c("GSM783029",	"GSM783030"	,"GSM783031",	"GSM783032","GSM783033",	"GSM783034"	,
              "GSM783035",	"GSM783036","GSM783037",	"GSM783038"  )]
dim(lung)
names(lung)
round(apply(lung, 2, summary))

breast <- read.table("GSE6569-GPL97_series_matrix.txt",skip = 57, header = TRUE, sep = "\t", row.names = 1,fill=TRUE)
#remove the last line from the matrix
breast<-breast[-22646,]  
breast<-breast[,c("GSM151910"	,"GSM151911",	"GSM151912",	"GSM151913",	"GSM151914",
                  "GSM151915",	"GSM151916"	,"GSM151917",	"GSM151918",	"GSM151919" )]
dim(breast)
names(breast)
round(apply(breast, 2, summary))

sum(is.na(diabetes))
sum(is.na(leukemia))
sum(is.na(parkinson))
sum(is.na(asthma))
sum(is.na(pancreatic))
sum(is.na(lung))
sum(is.na(breast))

logDiabetes <- log2(diabetes)
logLeukemia <- log2(leukemia)
logParkinson <- log2(parkinson)
logAsthma <- log2(asthma)
logPancreatic <- log2(pancreatic)
logLung <- log2(lung)
logBreast <- log2(breast)

allData = cbind(logDiabetes,logLeukemia,logParkinson,logAsthma,
                logPancreatic,logLung,logBreast,deparse.level = 2)

###########################################################
#Replace the rownames with the corresponding gene symbols.
##########################################################
allData<-as.matrix(allData)
rows <- rownames(allData)
rows[3200]
which(IDs[,1] == rows[3200])
IDs[which(IDs[,1] == rows[3200]), 2]
symb <- rep(0, length(rows))
for (i in 1:length(rows)) {
  #symb[i] <- as.character(IDs[which(IDs[,1] == rows[i]), 2])
  symb[i] <- paste(IDs[which(IDs[,1] == rows[i]), 2],"|",
                   IDs[which(IDs[,1] == rows[i]), 1])
}
rownames(allData) <- symb
rownames(logDiabetes) <- symb
rownames(logLeukemia) <- symb
rownames(logParkinson) <- symb
rownames(logAsthma) <- symb
rownames(logPancreatic) <- symb
rownames(logLung) <- symb
rownames(logBreast) <- symb

write.table(logDiabetes,'allData.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logDiabetes.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logLeukemia.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logParkinson.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logAsthma.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logPancreatic.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logLung.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logBreast.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)

nrm=function(x){ return( (x - mean(x))/sd(x))}
logDiabetes <- apply(logDiabetes,2,nrm)
logLeukemia <- apply(logLeukemia,2,nrm)
logParkinson <- apply(logParkinson,2,nrm)
logAsthma <- apply(logAsthma,2,nrm)
logPancreatic <- apply(logPancreatic,2,nrm)
logLung <- apply(logLung,2,nrm)
logBreast <- apply(logBreast,2,nrm)

allData = cbind(logDiabetes,logLeukemia,logParkinson,logAsthma,
                logPancreatic,logLung,logBreast,deparse.level = 2)


write.table(logDiabetes,'allDataN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logDiabetesN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logLeukemiaN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logParkinsonN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logAsthmaN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logPancreaticN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logLungN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)
write.table(logDiabetes,'logBreastN.txt', sep='\t', row.names = TRUE,
            col.names = TRUE)

all.sample.subtypes <- as.vector(projetab$disease)
all.sample.labels <- all.group.abbrev[all.sample.subtypes]
names(all.sample.labels) <- colnames(allData)

## Plot the expression profiles of the two selected genes 
all.sample.colors <- all.group.colors[as.vector(projetab$disease)]
names(all.sample.colors) <- names(allData)
sample(all.sample.colors,size=20)
gn1 <- 236
gn2 <- 1213
gn1x <- as.vector(as.matrix(allData[gn1,]))
gn2y <- as.vector(as.matrix(allData[gn2,]))
plot(gn1x,gn2y,
     col=all.sample.colors,
     type='n',
     panel.first=grid(col='black'), 
     main='two selected genes', 
     xlab=paste('gene', rownames(allData)[gn1]), ylab=paste('gene', rownames(allData)[gn2]))
text(gn1x, gn2y,labels=all.sample.labels,col=all.sample.colors,pch=0.5)
legend('topright',col=all.group.colors, 
       legend=names(all.group.colors),pch=0.5,cex=0.5,bg='white',inset=0.01,x.intersp=2,xjust=0,yjust=0)


## Compute gene-wise variance

geneLstTop=function(x,cnt){ 
  var.per.gene <- apply(x, 1, var)
  ## Inspect the distribution of gene-wise variance
  #hist(var.per.gene, breaks=100)
  ## Sort genes per decreasing variance
  genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)
  ### Print the 5 genes with highest variance
  head(genes.by.decr.var)
  geneLst <- names(genes.by.decr.var)[1:cnt]
  #colnames(geneLst)[1] = "Gene Name"
  return(list(geneLst,var.per.gene))
}

ret=geneLstTop(allData,5000)
geneLstTopAll=ret[1]
varGenesAll=ret[2]
ret=geneLstTopDiabetes=geneLstTop(logDiabetes,5000)
geneLstTopDiabetes=ret[1]
varGenesDiabetes=ret[2]
ret=geneLstTop(logLeukemia,5000)
geneLstTopLeukemia=ret[1]
varGenesLeukemia=ret[2]
ret=geneLstTop(logParkinson,5000)
geneLstTopParkinson=ret[1]
varGenesParkinson=ret[2]
ret=geneLstTop(logAsthma,5000)
geneLstTopAsthma=ret[1]
varGenesAsthma=ret[2]
ret=geneLstTop(logPancreatic,5000)
geneLstTopPancreatic=ret[1]
varGenesPancreatic=ret[2]
geneLstTopLung=geneLstTop(logLung,5000)
geneLstTopLung=ret[1]
varGenesLung=ret[2]
ret=geneLstTop(logBreast,5000)
geneLstTopBreast=ret[1]
varGenesBreast=ret[2]

dataTopAll = allData[geneLstTopAll[[1]],]
dataTopDiabetes = allData[geneLstTopDiabetes[[1]],1:10]
dataTopLeukemia = allData[geneLstTopLeukemia[[1]],11:20]
dataTopParkinson = allData[geneLstTopParkinson[[1]],21:30]
dataTopAsthma = allData[geneLstTopAsthma[[1]],31:40]
dataTopPancreatic = allData[geneLstTopPancreatic[[1]],41:50]
dataTopLung = allData[geneLstTopLung[[1]],51:60]
dataTopBreast = allData[geneLstTopBreast[[1]],61:70]

dim(dataTopAll)
dim(dataTopDiabetes)
dim(dataTopLeukemia)
dim(dataTopParkinson)
dim(dataTopAsthma)
dim(dataTopPancreatic)
dim(dataTopLung)
dim(dataTopBreast)


rownames(dataTopAll[1:10,])

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

sink('generanks.txt')
print(tbl,right=F)
sink()




