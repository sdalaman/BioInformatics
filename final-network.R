# Correlations with significance levels
library(Hmisc)
library(igraph)
library(NetIndices)
library(sand)

numOfGenes = 500
crGenesAll <- cor(t(dataTopAll[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesDiabetes <- cor(t(dataTopDiabetes[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesLeukemia <- cor(t(dataTopLeukemia[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesAsthma <- cor(t(dataTopAsthma[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesPancreatic <- cor(t(dataTopPancreatic[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesLung <- cor(t(dataTopLung[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesParkinson <- cor(t(dataTopParkinson[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesBreast <- cor(t(dataTopBreast[1:numOfGenes,]), use="complete.obs", method="pearson")


createNetworkFile=function(crData,fname){
  NRows <- dim(crData)[1]
  NCols <- dim(crData)[2]
  fn <- file(fname,"w")
  writeLines("GeneId1,GeneId2,Corr",fn)
  threshold <- 0.85
  rwcnt <- 5000
  for (rw in c(2:NRows)){
    for (cl in c(1:rw)){
        if ((rw > cl) && (abs(crData[rw,cl]) >= threshold)){
          rwcnt <- rwcnt + 1
          line <- paste(rownames(crData)[rw],",",colnames(crData)[cl],",",crData[rw,cl],sep="")
          writeLines(line,fn)
        }
    }
    cat("row ",rw,"\n")
  }
  cat(rwcnt," line written\n")
  flush(fn)
  close(fn)
  cat(fname," created")
}

createNetworkFile(crGenesAll,"gene_links_all.csv")
createNetworkFile(crGenesDiabetes,"gene_links_diabetes.csv")
createNetworkFile(crGenesLeukemia,"gene_links_leukemia.csv")
createNetworkFile(crGenesAsthma,"gene_links_asthma.csv")
createNetworkFile(crGenesParkinson,"gene_links_parkinson.csv")
createNetworkFile(crGenesPancreatic,"gene_links_pancreatic.csv")
createNetworkFile(crGenesLung,"gene_links_lung.csv")
createNetworkFile(crGenesBreast,"gene_links_breast.csv")

all.links.data<-read.csv("gene_links_all.csv")
diabetes.links.data<-read.csv("gene_links_diabetes.csv")
leukemia.links.data<-read.csv("gene_links_leukemia.csv")
parkinson.links.data<-read.csv("gene_links_parkinson.csv")
pancreatic.links.data<-read.csv("gene_links_pancreatic.csv")
asthma.links.data<-read.csv("gene_links_asthma.csv")
lung.links.data<-read.csv("gene_links_lung.csv")
breast.links.data<-read.csv("gene_links_breast.csv")

# Column names for data
#colnames(all.links.data)
# Convert the data into a graph object using the first 2 columns of the dataset as an edgelist
all.graph<-graph.edgelist(as.matrix(all.links.data[,1:2]),directed = FALSE)
diabetes.graph<-graph.edgelist(as.matrix(diabetes.links.data[,1:2]),directed = FALSE)
leukemia.graph<-graph.edgelist(as.matrix(leukemia.links.data[,1:2]),directed = FALSE)
parkinson.graph<-graph.edgelist(as.matrix(parkinson.links.data[,1:2]),directed = FALSE)
pancreatic.graph<-graph.edgelist(as.matrix(pancreatic.links.data[,1:2]),directed = FALSE)
asthma.graph<-graph.edgelist(as.matrix(asthma.links.data[,1:2]),directed = FALSE)
lung.graph<-graph.edgelist(as.matrix(lung.links.data[,1:2]),directed = FALSE)
breast.graph<-graph.edgelist(as.matrix(breast.links.data[,1:2]),directed = FALSE)


#write.graph(all.graph, file='all.graph.txt', format="dot")



all.adjmatrix<-get.adjacency(all.graph,sparse=F)
diabetes.adjmatrix<-get.adjacency(diabetes.graph,sparse=F)
leukemia.adjmatrix<-get.adjacency(leukemia.graph,sparse=F)
parkinson.adjmatrix<-get.adjacency(parkinson.graph,sparse=F)
pancreatic.adjmatrix<-get.adjacency(pancreatic.graph,sparse=F)
asthma.adjmatrix<-get.adjacency(asthma.graph,sparse=F)
lung.adjmatrix<-get.adjacency(lung.graph,sparse=F)
breast.adjmatrix<-get.adjacency(breast.graph,sparse=F)


# Get the basic network indices from the matrices with GenInd()
#ind.all<-GenInd(all.adjmatrix)
# Now to plot these two webs to get a feel for what we are dealing with

#par(mar=c(.1,.1,.1,.1))
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout.circle,main="Circle Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_as_tree,main="Tree Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout.fruchterman.reingold,main="fruchterman.reingold Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_as_star,main="Star Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_in_circle,main="Circle Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_nicely,main="Nicely Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout.davidson.harel,main="davidson.harel Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_dh,main="Dh Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_gem,main="Gem Layout" )
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_graphopt,main="graphopt Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_on_grid,main="Grid Layout") 
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_mds,main="Mds Layout") 
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_components,main="components Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_sugiyama,main="sugiyama Layout")  
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_on_sphere,main="Sphere Layout") 
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_randomly,main="Randomly Layout") 
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_fr,main="Fr Layout") 
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_kk,main="Kk Layout")
#plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
#            layout=layout_with_lgl,main="Lgl Layout") 


plot.igraph(all.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Combined Samples (Top 500 Genes in variance order)")
plot.igraph(diabetes.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Diabetes Samples (Top 500 Genes in variance order)")
plot.igraph(leukemia.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Leukemia Samples (Top 500 Genes in variance order)")
plot.igraph(parkinson.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Parkinson Samples (Top 500 Genes in variance order)")
plot.igraph(pancreatic.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Pancreatic Cancer Samples (Top 500 Genes in variance order)")
plot.igraph(asthma.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Asthma Samples (Top 500 Genes in variance order)")
plot.igraph(lung.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Lung Cancer Samples (Top 500 Genes in variance order)")
plot.igraph(breast.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,
            layout=layout.fruchterman.reingold,main="Breast Cancer Samples (Top 500 Genes in variance order)")


# We can use a community detection algorithm to determine the most densely connected nodes in a graph.

#all.graph.community<-walktrap.community(all.graph)
#sizes(all.graph.community)
#barplot(sizes(all.graph.community))
#plot(all.graph.community,all.graph,layout=layout.fruchterman.reingold)
# This algorithm uses random walks to find the most densely connected subgraphs.
#all.members<-membership(all.graph.community)
# The members() function picks out the membership vector (list of nodes in the most densely connected subgraph) 
# from the communtiy object (e.g., walktrap community).
#par(mar=c(.1,.1,.1,.1))    # sets the edges of the plotting area
#plot.igraph(all.graph,
#            layout=layout.fruchterman.reingold,
#            #layout=layout_with_dh,
#            vertex.size=10,
#            vertex.label.cex=.5,
#            edge.arrow.size=.5,
#            mark.groups=list(all.members),
#            mark.col="green"
#)

# We can use a community detection algorithm to determine the most densely connected nodes in a graph.
getCommunityAndMembers=function(disease.graph,sample.type){
    disease.graph.community<-walktrap.community(disease.graph)
    #sizes(disease.graph.community)
    barplot(sizes(disease.graph.community),main=paste(sample.type,"-Network Community Sizes Plot"))
    plot(disease.graph.community,disease.graph,layout=layout.fruchterman.reingold,vertex.label="",
         main=paste(sample.type,"-Network Community Size Plot"))
  # This algorithm uses random walks to find the most densely connected subgraphs.
    disease.members<-membership(disease.graph.community)
  # The members() function picks out the membership vector (list of nodes in the most densely connected subgraph) 
  # from the communtiy object (e.g., walktrap community).
#    par(mar=c(.1,.1,.1,.1))    # sets the edges of the plotting area
#    plot.igraph(all.graph,
#              layout=layout.fruchterman.reingold,
#              #layout=layout_with_dh,
#              vertex.size=10,
#              vertex.label.cex=.5,
#              vertex.label="",
#              edge.arrow.size=.5,
#              mark.groups=list(disease.members),
#              mark.col="green",
#              main=paste(sample.type,"-The Most Densely Connected Subgraph")
#    ) 
    return(list(disease.graph.community,disease.members))
}

ret=getCommunityAndMembers(all.graph,"Combined Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]] 
ret=getCommunityAndMembers(diabetes.graph,"Diabetes Samples")
mostDenseDiabetes=ret[[2]]
communityDiabetes=ret[[1]]
ret=getCommunityAndMembers(leukemia.graph,"Leukemia Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]
ret=getCommunityAndMembers(asthma.graph,"Asthma Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]
ret=getCommunityAndMembers(parkinson.graph,"Parkinson Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]
ret=getCommunityAndMembers(pancreatic.graph,"Pancreatic Cancer Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]
ret=getCommunityAndMembers(lung.graph,"Lung Cancer Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]
ret=getCommunityAndMembers(breast.graph,"Breast Cancer Samples")
mostDenseAll=ret[[2]]
communityAll=ret[[1]]


ceb <- cluster_edge_betweenness(all.graph) 
dendPlot(ceb, mode="dendrogram")
plot(ceb, all.graph,vertex.label="") 

# The "GenInd()" function requires an input of an adjacency matrix
all.graph.adj<-get.adjacency(all.graph,sparse=F)
# in older versions of igraph the default was sparse=F,
# but now you must specify, other wise you get a matrix of 1s and .s
all.graph.properties<-GenInd(all.graph.adj)
all.graph.properties$N            #number of nodes
all.graph.properties$Ltot        #number of links
all.graph.properties$LD        #link density (average # of links per node)
all.graph.properties$C            #the connectance of the graph
# This function measures connectance as L/(N*(N-1)) where L is links, and N is nodes
# Connectance can also be calculated as L/(N^2)

#in.deg.all.graph<-degree(all.graph,v=V(all.graph),mode="in")
#out.deg.all.graph<-degree(all.graph,v=V(all.graph),mode="out")
#all.deg.all.graph<-degree(all.graph,v=V(all.graph),mode="all")

#sort(in.deg.all.graph,decreasing=TRUE)[1:5]
#sort(out.deg.all.graph,decreasing=TRUE)[1:5]
#sort(all.deg.all.graph,decreasing=TRUE)[1:5]

# Degree distribution is the cumulative frequency of nodes with a given degree
# this, like degree() can be specified as "in", "out", or "all"
#all.deg.distr<-degree.distribution(all.graph,cumulative=T,mode="all")

# Using the power.law.fit() function I can fit a power law to the degree distribution
#power<-power.law.fit(all.deg.all.graph)
#cat(power$KS.stat)  # smaller better
#if(power$KS.p < 0.05){
#  cat("The data can not be from fitted power-law")
#} else {
#  cat("The data can be from fitted power-law")
#}
# The output of the power.law.fit() function tells me what the exponent of the power law is ($alpha)
# and the log-likelihood of the parameters used to fit the power law distribution ($logLik)
# Also, it performs a Kolmogov-Smirnov test to test whether the given degree distribution could have
# been drawn from the fitted power law distribution.
# The function thus gives me the test statistic ($KS.stat) and p-vaule ($KS.p) for that test

# Then I can plot the degree distribution
#plot(all.deg.distr,log="xy",
#     ylim=c(.01,10),
#     bg="black",pch=21,
#     xlab="Degree",
#     ylab="Cumulative Frequency")

# And the expected power law distribution
#lines(1:20,10*(1:20)^((-power$alpha)+1))

diseaseNetworkDistribution=function(disease.graph,sample.type){
  disease.deg.distr<-degree.distribution(disease.graph,cumulative=T,mode="all")
    # Using the power.law.fit() function I can fit a power law to the degree distribution
  all.deg.disease.graph<-degree(disease.graph,v=V(disease.graph),mode="all")
  power<-power.law.fit(all.deg.disease.graph)
  cat(power$KS.stat)  # smaller better
  if(power$KS.p < 0.05){
    cat("The data can not be from fitted power-law")
  } else {
    cat("The data can be from fitted power-law")
  }
  # The output of the power.law.fit() function tells me what the exponent of the power law is ($alpha)
  # and the log-likelihood of the parameters used to fit the power law distribution ($logLik)
  # Also, it performs a Kolmogov-Smirnov test to test whether the given degree distribution could have
  # been drawn from the fitted power law distribution.
  # The function thus gives me the test statistic ($KS.stat) and p-vaule ($KS.p) for that test
  
  # Then I can plot the degree distribution
  plot(disease.deg.distr,log="xy",
       ylim=c(.01,10),
       bg="black",pch=21,
       xlab="Degree",
       ylab="Cumulative Frequency",
       main=paste(sample.type,"-Network Distribution Graph"))
  
  # And the expected power law distribution
  lines(1:20,10*(1:20)^((-power$alpha)+1))
}

diseaseNetworkDistribution(all.graph,"Combined Samples")
diseaseNetworkDistribution(diabetes.graph,"Diabetes Samples")
diseaseNetworkDistribution(leukemia.graph,"Leukemia Samples")
diseaseNetworkDistribution(asthma.graph,"Asthma Samples")
diseaseNetworkDistribution(parkinson.graph,"Parkinson Samples")
diseaseNetworkDistribution(pancreatic.graph,"Pancreatic Cancer Samples")
diseaseNetworkDistribution(lung.graph,"Lung Cancer Samples")
diseaseNetworkDistribution(breast.graph,"Breast Cancer Samples")


# Diameter is essentially the longest path between two vertices
#diameter(all.graph)
# Gives me the length of the diameter while

#nodes.diameter<-get.diameter(all.graph)
# Gives me the labels for each node that participates in the diameter

# I can look at the diameter graphically also
# First I will define the node and edge attributes
#V(all.graph)$color<-"skyblue"
# I want all the nodes to be skyblue
#V(all.graph)$size<-7
# I want all the nodes to be size=7
#V(all.graph)[nodes.diameter]$color<-"darkgreen"
#V(all.graph)[nodes.diameter]$size<-10
#V(all.graph)[nodes.diameter]$label.color<-"white"
# but the nodes in the diameter should be darkgreen and larger than the rest
# with a white label instead of black
# this will make the diameter pop out of the larger network
#E(all.graph)$color<-"grey"
# all non-diameter edges will be grey
#E(all.graph,path=nodes.diameter)$color<-"darkgreen"
#E(all.graph,path=nodes.diameter)$width<-2
# Edges in the diameter will be darkgreen and a little extra wide

# If you do not set the attributes of all of the nodes and edges then it will
# default such that you only see what you have defined

# Now when I plot the diameter will be larger than everything else, and darkgreen instead
# of grey/blue
#par(mar=c(.1,.1,.1,.1),lwd=0.1)
#plot.igraph(all.graph,
#            layout=layout.fruchterman.reingold,
#            vertex.label.cex=.5,
#            edge.arrow.size=.5)


diseaseLongestPath=function(disease.graph,sample.type){
  # Diameter is essentially the longest path between two vertices
  diameter(disease.graph)
  # Gives me the length of the diameter while
  
  nodes.diameter<-get.diameter(disease.graph)
  # Gives me the labels for each node that participates in the diameter
  
  # I can look at the diameter graphically also
  # First I will define the node and edge attributes
  V(disease.graph)$color<-"skyblue"
  # I want all the nodes to be skyblue
  V(disease.graph)$size<-7
  # I want all the nodes to be size=7
  V(disease.graph)[nodes.diameter]$color<-"darkgreen"
  V(disease.graph)[nodes.diameter]$size<-10
  V(disease.graph)[nodes.diameter]$label.color<-"white"
  # but the nodes in the diameter should be darkgreen and larger than the rest
  # with a white label instead of black
  # this will make the diameter pop out of the larger network
  E(disease.graph)$color<-"grey"
  # all non-diameter edges will be grey
  E(disease.graph,path=nodes.diameter)$color<-"darkgreen"
  E(disease.graph,path=nodes.diameter)$width<-2
  # Edges in the diameter will be darkgreen and a little extra wide
  
  # If you do not set the attributes of all of the nodes and edges then it will
  # default such that you only see what you have defined
  
  # Now when I plot the diameter will be larger than everything else, and darkgreen instead
  # of grey/blue
  par(mar=c(.1,.1,.1,.1),lwd=0.1)
  plot.igraph(disease.graph,
              layout=layout.fruchterman.reingold,
              vertex.label.cex=.5,
              edge.arrow.size=.5,
              main=paste(sample.type,"-Network Longest Path (Diameter)"))
  
}

diseaseLongestPath(all.graph,"Combined Samples")
diseaseLongestPath(diabetes.graph,"Diabetes Samples")
diseaseLongestPath(leukemia.graph,"Leukemia Samples")
diseaseLongestPath(asthma.graph,"Asthma Samples")
diseaseLongestPath(parkinson.graph,"Parkinson Samples")
diseaseLongestPath(pancreatic.graph,"Pancreatic Cancer Samples")
diseaseLongestPath(lung.graph,"Lung Cancer Samples")
diseaseLongestPath(breast.graph,"Breast Cancer Samples")

# Clustering coefficient is the proportion of
# a nodes neighbors that can be reached by other neighbors
# in igraph this property is apparently called "transitivity"

clusteringCoff=function(disease.graph,sample.type){
  global.trans <- transitivity(disease.graph)
  # gives the clustering coefficient of the whole network
  local.trns <- transitivity(disease.graph,type="local")
  # gives the clustering coefficient of each node  : ratio of triangles to connected triples each vertex is part of.
  hist(local.trns*1000, breaks=100, main=paste(sample.type,"-clustering coefficient of each node"))
  return(global.trans)
}

all.global.trans=clusteringCoff(all.graph,"Combined Samples")
diabetes.global.trans=clusteringCoff(diabetes.graph,"Diabetes Samples")
leukemia.global.trans=clusteringCoff(leukemia.graph,"Leukemia Samples")
asthma.global.trans=clusteringCoff(asthma.graph,"Asthma Samples")
parkinson.global.trans=clusteringCoff(parkinson.graph,"Parkinson Samples")
pancreatic.global.trans=clusteringCoff(pancreatic.graph,"Pancreatic Cancer Samples")
lung.global.trans=clusteringCoff(lung.graph,"Lung Cancer Samples")
breast.global.trans=clusteringCoff(breast.graph,"Breast Cancer Samples")

# Betweenness is the number of shortest paths between two nodes that go through each node of interest

nodeBetwenness=function(disease.graph,sample.type,node.number){
  disease.graph.betweenness<-betweenness(disease.graph,v=V(disease.graph))
  sorted.bw <- sort(disease.graph.betweenness,decreasing=TRUE)[1:node.number]
  vlist <- rep(0, length(sorted.bw))
  for (i in 1:length(labels(sorted.bw))) {
    vlist[i] <- which(labels(disease.graph.betweenness) == labels(sorted.bw)[i])
  }
  col <- rep("cyan", vcount(disease.graph))
  col[vlist] <- "red"
  plot(disease.graph, vertex.label="",vertex.color=col,layout=layout.fruchterman.reingold,
       main=paste(sample.type,"-nodes with highest betweenness "))
  return(disease.graph.betweenness)
}

all.betweenness=nodeBetwenness(all.graph,"Combined Samples",10)
diabetes.betweenness=nodeBetwenness(diabetes.graph,"Diabetes Samples",10)
leukemia.betweenness=nodeBetwenness(leukemia.graph,"Leukemia Samples",10)
asthma.betweenness=nodeBetwenness(asthma.graph,"Asthma Samples",10)
parkinson.betweenness=nodeBetwenness(parkinson.graph,"Parkinson Samples",10)
pancreatic.betweenness=nodeBetwenness(pancreatic.graph,"Pancreatic Cancer Samples",10)
lung.betweenness=nodeBetwenness(lung.graph,"Lung Cancer Samples",10)
breast.betweenness=nodeBetwenness(breast.graph,"Breast Cancer Samples",10)

nodeDegrees=function(disease.graph,sample.type){
  dg <- degree(disease.graph,mode="all")
  plot(disease.graph, vertex.label="",vertex.size=(round(10*abs(dg - mean(dg))/sd(dg))),layout=layout.fruchterman.reingold
     , main=paste(sample.type,"Plot of node degree"))
  hist(dg, breaks=1:vcount(all.graph)-1, main=paste(sample.type,"Histogram of node degree"))
  return(dg)
}

all.degrees=nodeDegrees(all.graph,"Combined Samples")
diabetes.degrees=nodeDegrees(diabetes.graph,"Diabetes Samples")
leukemia.degrees=nodeDegrees(leukemia.graph,"Leukemia Samples")
asthma.degrees=nodeDegrees(asthma.graph,"Asthma Samples")
parkinson.degrees=nodeDegrees(parkinson.graph,"Parkinson Samples")
pancreatic.degrees=nodeDegrees(pancreatic.graph,"Pancreatic Cancer Samples")
lung.degrees=nodeDegrees(lung.graph,"Lung Cancer Samples")
breast.degress=nodeDegrees(breast.graph,"Breast Cancer Samples")


hubScores=function(disease.graph,sample.type){
  hs <- hub_score(disease.graph, weights=NA)$vector
  #as <- authority_score(all.graph, weights=NA)$vector
  #par(mfrow=c(1,2))
  plot(disease.graph, vertex.label="",vertex.size=hs*15, main=paste(sample.type,"-Hubs"))
  #plot(all.graph, vertex.label="",vertex.size=as*10, main="Authorities")
  sorted.hs <- sort(hs,decreasing=TRUE)[1:5]
  vlist <- rep(0, length(sorted.hs))
  for (i in 1:length(labels(sorted.hs))) {
    vlist[i] <- which(labels(hs) == labels(sorted.hs)[i])
  }
  col <- rep("cyan", vcount(all.graph))
  col[vlist] <- "red"
  plot(disease.graph, vertex.label="",vertex.color=col, main=paste(sample.type,"-nodes with highest hub score "))
  return(hs)
}

all.hub.scores=hubScores(all.graph,"Combined Samples")
diabetes.hub.scores=hubScores(diabetes.graph,"Diabetes Samples")
leukemia.hub.scores=hubScores(leukemia.graph,"Leukemia Samples")
asthma.hub.scores=hubScores(asthma.graph,"Asthma Samples")
parkinson.hub.scores=hubScores(parkinson.graph,"Parkinson Samples")
pancreatic.hub.scores=hubScores(pancreatic.graph,"Pancreatic Cancer Samples")
lung.hub.scores=hubScores(lung.graph,"Lung Cancer Samples")
breast.hub.scores=hubScores(breast.graph,"Breast Cancer Samples")



library(knitr)
tbl <- NULL
a1 <- sort(all.hub.scores,decreasing=TRUE)[1:5]
a2 <- sort(diabetes.hub.scores,decreasing=TRUE)[1:5]
a3 <- sort(leukemia.hub.scores,decreasing=TRUE)[1:5]
a4 <- sort(parkinson.hub.scores,decreasing=TRUE)[1:5]
a5 <- sort(asthma.hub.scores,decreasing=TRUE)[1:5]
a6 <- sort(pancreatic.hub.scores,decreasing=TRUE)[1:5]
a7 <- sort(lung.hub.scores,decreasing=TRUE)[1:5]
a8 <- sort(breast.hub.scores,decreasing=TRUE)[1:5]
ln <- c("ALL","Diabetes","Leukemia","Parkinson","Asthma","Pancreatic","Lung","Breast")
tbl <- rbind(ln)
for(i in 1:5){
  ln <- c(labels(a1[i]),labels(a2[i]),labels(a3[i]),labels(a4[i]),labels(a5[i])
          ,labels(a6[i]),labels(a7[i]),labels(a8[i]))
  tbl <- rbind(tbl,ln)
}
rownames(tbl)[1] <- "Disease"
for(i in 2:6)
  rownames(tbl)[i] <- paste("Gene ",as.character(i-1))
print(tbl)

sink('hubtbl.txt')
print(tbl,right=F)
sink()


############################################################################33

numOfGenes = 100
crGenesAll2 <- cor(t(dataTopAll[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesDiabetes2 <- cor(t(dataTopDiabetes[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesLeukemia2 <- cor(t(dataTopLeukemia[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesAsthma2 <- cor(t(dataTopAsthma[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesPancreatic2 <- cor(t(dataTopPancreatic[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesLung2 <- cor(t(dataTopLung[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesParkinson2 <- cor(t(dataTopParkinson[1:numOfGenes,]), use="complete.obs", method="pearson")
crGenesBreast2 <- cor(t(dataTopBreast[1:numOfGenes,]), use="complete.obs", method="pearson")

createNetworkFile(crGenesAll2,"gene_links_all2.csv")
createNetworkFile(crGenesDiabetes2,"gene_links_diabetes2.csv")
createNetworkFile(crGenesLeukemia2,"gene_links_leukemia2.csv")
createNetworkFile(crGenesAsthma2,"gene_links_asthma2.csv")
createNetworkFile(crGenesParkinson2,"gene_links_parkinson2.csv")
createNetworkFile(crGenesPancreatic2,"gene_links_pancreatic2.csv")
createNetworkFile(crGenesLung2,"gene_links_lung2.csv")
createNetworkFile(crGenesBreast2,"gene_links_breast2.csv")

all.links.data2<-read.csv("gene_links_all2.csv")
diabetes.links.data2<-read.csv("gene_links_diabetes2.csv")
leukemia.links.data2<-read.csv("gene_links_leukemia2.csv")
parkinson.links.data2<-read.csv("gene_links_parkinson2.csv")
pancreatic.links.data2<-read.csv("gene_links_pancreatic2.csv")
asthma.links.data2<-read.csv("gene_links_asthma2.csv")
lung.links.data2<-read.csv("gene_links_lung2.csv")
breast.links.data2<-read.csv("gene_links_breast2.csv")

all.graph2<-graph.edgelist(as.matrix(all.links.data2[,1:2]),directed = FALSE)
diabetes.graph2<-graph.edgelist(as.matrix(diabetes.links.data2[,1:2]),directed = FALSE)
leukemia.graph2<-graph.edgelist(as.matrix(leukemia.links.data2[,1:2]),directed = FALSE)
parkinson.graph2<-graph.edgelist(as.matrix(parkinson.links.data2[,1:2]),directed = FALSE)
pancreatic.graph2<-graph.edgelist(as.matrix(pancreatic.links.data2[,1:2]),directed = FALSE)
asthma.graph2<-graph.edgelist(as.matrix(asthma.links.data2[,1:2]),directed = FALSE)
lung.graph2<-graph.edgelist(as.matrix(lung.links.data2[,1:2]),directed = FALSE)
breast.graph2<-graph.edgelist(as.matrix(breast.links.data2[,1:2]),directed = FALSE)


all.clq <- cliques(all.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
largest_cliques(all.graph2) # cliques with max number of nodes
vcol <- rep("grey80", vcount(all.graph2))
vcol[unlist(largest_cliques(all.graph2))] <- "gold"
plot(as.undirected(all.graph2), vertex.label=V(all.graph2)$name, vertex.color=vcol,
     main=paste("Combined","-clique with highest score "))

a1 <- largest_cliques(all.graph2)

all.clq <- cliques(diabetes.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
largest_cliques(all.graph2) # cliques with max number of nodes
vcol <- rep("grey80", vcount(all.graph2))
vcol[unlist(largest_cliques(all.graph2))] <- "gold"
plot(as.undirected(all.graph2), vertex.label=V(all.graph2)$name, vertex.color=vcol,
     main=paste("Diabetes","-clique with highest score "))

a2 <- largest_cliques(diabetes.graph2)

all.clq <- cliques(leukemia.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a3 <- largest_cliques(leukemia.graph2) # cliques with max number of nodes

all.clq <- cliques(parkinson.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a4 <- largest_cliques(parkinson.graph2) # cliques with max number of nodes

all.clq <- cliques(asthma.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a5 <- largest_cliques(asthma.graph2) # cliques with max number of nodes


all.clq <- cliques(pancreatic.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a6 <- largest_cliques(pancreatic.graph2) # cliques with max number of nodes


all.clq <- cliques(lung.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a7 <- largest_cliques(lung.graph2) # cliques with max number of nodes


all.clq <- cliques(breast.graph2) # list of cliques       
sapply(all.clq, length) # clique sizes
a8 <- largest_cliques(breast.graph2) # cliques with max number of nodes



tbl <- NULL
ln <- c("ALL","Diabetes","Leukemia","Parkinson","Asthma","Pancreatic","Lung","Breast")
tbl <- rbind(ln)
for(i in 1:5){
  ln <- c(labels(a1[i]),labels(a2[i]),labels(a3[i]),labels(a4[i]),labels(a5[i])
          ,labels(a6[i]),labels(a7[i]),labels(a8[i]))
  tbl <- rbind(tbl,ln)
}
rownames(tbl)[1] <- "Disease"
for(i in 2:6)
  rownames(tbl)[i] <- paste("Gene ",as.character(i-1))
print(tbl)

sink('clique')
print(tbl,right=F)
sink()





all.graph.edge.betweenness<-edge.betweenness(all.graph,e=E(all.graph))

# Closeness refers to how connected a node is to its neighbors

all.graph.closeness<-closeness(all.graph,vids=V(all.graph))

# Clustering coefficient, betweenness, and closeness
# all describe the small world properties of the network.
# A network with small world properties is one in which
# it takes a relatively short path to get from one node to the next
# (e.g., six degrees of separation)

# Every graph can be decomposed into its component n-node subgraphs.
# In particular there are 13 unique ways to arrange 3 nodes in directed graphs.
# Here are the adjacency matrices for each of the 13 subgraphs
s1<-matrix(c(0,1,0,0,0,1,0,0,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,0,0,1,0,0,0),nrow=3,ncol=3)
s3<-matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3)
s4<-matrix(c(0,0,1,0,0,1,0,0,0),nrow=3,ncol=3)
s5<-matrix(c(0,1,1,0,0,0,0,0,0),nrow=3,ncol=3)
d2<-matrix(c(0,1,1,1,0,1,0,0,0),nrow=3,ncol=3)
d1<-matrix(c(0,1,1,0,0,1,0,1,0),nrow=3,ncol=3)
d3<-matrix(c(0,0,1,1,0,0,1,0,0),nrow=3,ncol=3)
d4<-matrix(c(0,0,0,1,0,1,0,1,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,0,0,1,1,0,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,0,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# I then make the 13 matrices into a list
subgraph3.mat<-list(s1,s2,s3,s4,s5,d1,d2,d3,d4,d5,d6,d7,d8)
# And convert the matrices into graph objects
subgraph3.graph<-lapply(subgraph3.mat,graph.adjacency,mode="undirected")

# Here I have created a simple for loop to go through the list of subgraphs
# and count how many times that subgraph appears in the larger test.graph
subgraph.count<-c()
for(i in 1:13){
  subgraph.count[i]<-
    graph.count.subisomorphisms.vf2(all.graph,subgraph3.graph[[i]])
}


