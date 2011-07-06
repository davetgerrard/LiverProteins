

if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified


#### Go through func clusters and pull in topGO results into a single table across BP,MF,CC


elimCutOff <- 0.01


GOdata.BP <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "BP" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.MF <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "MF" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.CC <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "CC" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )




seedList




topDiffGenes <- function(allScore) {
  return(allScore )
}



summaryResultsList <- list()
#PC1
i<-1
for( i in 1:5)  {
compHead <- paste("Comp.",i,sep="")

#combResults.pc1 <-data.frame()
summaryResultsList[[i]] <- data.frame()

for(thisGOgraph in goGraphs )  {

geneList <- ubi.pca.5.scores[,compHead]
names(geneList) <- ubi.pca.5.scores$spAccession
GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = thisGOgraph ,
              allGenes = geneList,
		  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
		GO2genes=go2prot 
		               )


#GOdata <- updateGenes(GOdata,geneList,topDiffGenes)

test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
resultWilcox <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcox <- getSigGroups(GOdata, test.stat)

#resultElimWilcoxAdj <- resultElimWilcox 
#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
resultElimKS <- getSigGroups(GOdata, test.stat)

#resultElimKSAdj <- resultElimKS
#score(resultElimKSAdj ) <- qvalue(score(resultElimKS))$qvalue

# change scores to absolute values for Greater than test
geneList <- abs(ubi.pca.5.scores[,compHead])
names(geneList) <- ubi.pca.5.scores$spAccession
GOdata <- updateGenes(GOdata,geneList,topDiffGenes)
test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
resultWilcoxAbs <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)

#resultElimWilcoxAbsAdj <- resultElimWilcoxAbs 
#score(resultElimWilcoxAbsAdj ) <- qvalue(score(resultElimWilcoxAbs ))$qvalue



Wilcox <- data.frame(Wilcox=score(resultWilcox ),goTerm=names(score(resultWilcox)))
elimWilcox <- data.frame(elimWilcox=score(resultElimWilcox),goTerm=names(score(resultElimWilcox)))
KS <- data.frame(KS=score(resultKS),goTerm=names(score(resultKS)))
elimKS <- data.frame(elimKS=score(resultElimKS),goTerm=names(score(resultElimKS)))
absWilcox <- data.frame(absWilcox=score(resultWilcoxAbs),goTerm=names(score(resultWilcoxAbs)))
elimAbsWilcox <- data.frame(elimAbsWilcox=score(resultElimWilcoxAbs),goTerm=names(score(resultElimWilcoxAbs)))

thisGoResults <- NULL
thisGoResults <- merge(Wilcox,elimWilcox,by="goTerm")
thisGoResults <- merge(thisGoResults ,KS ,by="goTerm")
thisGoResults <- merge(thisGoResults ,elimKS ,by="goTerm")
thisGoResults <- merge(thisGoResults ,absWilcox ,by="goTerm")
thisGoResults <- merge(thisGoResults ,elimAbsWilcox ,by="goTerm")
thisGoResults$ontology <- thisGOgraph
thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(thisGoResults$goTerm, ontology(GOdata)))

goGroup <- as.character(thisGoResults$goTerm)
thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdata,x)[[1]])))

summaryResultsList[[i]] <- rbind(summaryResultsList[[i]],thisGoResults)

} # end of GO loops
} # end of PC loops


#### Use seedList and go through output table. 

#combResults.pc1
seedList

#### list proteins for a GO term, accounting for ontology
listProtsInGo <- function(goTerm,ontology)  {
	if(ontology == "BP") {
		genesInTerm(GOdata.BP,goTerm)[[1]]
	} else if(ontology == "MF") {
		genesInTerm(GOdata.MF,goTerm)[[1]]
	} else if(ontology == "CC") {
		genesInTerm(GOdata.CC,goTerm)[[1]]
	}
}


for(j in 1:5) {

combResults <- summaryResultsList[[j]]

clusterSummaries <- data.frame()
clusterResults <- list()
for(i in 1:length(seedList))  {

	thisCluster <- i

	# get GO terms for these proteins 
	
	# intersect for tested GO terms.

	# how many are left?


	#how many proteins must fit the GO term. 
	propTermShare <- 0.5
	topClusterTerms <- names(which(sort(table(as.character(unlist(prot2go[seedList[[thisCluster]]]))),decreasing=T ) >= ((length(seedList[[thisCluster ]])*propTermShare) )))
	allTerms <- names(table(as.character(unlist(prot2go[seedList[[thisCluster ]]]))))
	
	#mainGoTermsIndex <- which(colSums(binaryGrid.detected[seedList[[thisCluster],]) > nodeSizeValue)
	#mainGoTerms <- colnames(binaryGrid.detected)[mainGoTermsIndex]
	
	# only want terms that are included in results.
	topClusterTerms <- intersect(topClusterTerms,combResults$goTerm)
	allTerms <- intersect(allTerms ,combResults$goTerm)  # redundant?
	#allTerms <- mainGoTerms

	#allWeights <- 
	#xx[topClusterTerms ]

	termIndex <-  match(allTerms ,combResults$goTerm)
	thisTermResults <- combResults[termIndex,]
	clusterSize <- length(seedList[[thisCluster]])
	thisTermResults$clusNumb <- apply(thisTermResults,1, FUN= function(x) length(intersect(seedList[[thisCluster]],listProtsInGo(goTerm=x["goTerm"],ontology=x["ontology"])) ))
	thisTermResults <- thisTermResults[thisTermResults$clusNumb >= 10,]
	thisTermResults <- thisTermResults[order(thisTermResults$elimWilcox),]

	elimWilcox.mean <- mean(thisTermResults$elimWilcox)
	elimWilcox.weigthedMean <- weighted.mean(thisTermResults$elimWilcox,( thisTermResults$clusNumb/thisTermResults$number))
	
	elimWilcox.geoMean <- prod(thisTermResults$elimWilcox)^(1/(length(thisTermResults$elimWilcox)))  # geometric mean
	cluster.Enrichment <- -log10(elimWilcox.geoMean)
	
	summary <- data.frame(elimWilcoxMean=elimWilcox.mean, elimWilcoxWeightedMean=elimWilcox.weigthedMean, elimWilcoxGeoMean=elimWilcox.geoMean,clusterEnrichment=cluster.Enrichment )
	clusterSummaries <- rbind(clusterSummaries,summary)
	clusterResults[[i]] <- thisTermResults
}		

if(output)  {
printOrder <- order(clusterSummaries$elimWilcoxGeoMean)
outputFile <- paste("SummaryFuncClust_PC",j,".pdf",sep="")

pdf(file=outputFile , width=15,height=9)
header <- paste("Principal component:", j)
textplot(header)
for(i in 1:length(printOrder )) {
	thisCluster <- printOrder[i]
	thisTermResults <- clusterResults[[thisCluster]]
	if(nrow(thisTermResults) < 1) {next;}
	clusterSize <- length(seedList[[thisCluster]])
	layout(matrix(c(1,2,3,3), ncol=2,byrow=T),heights=c(1,4)) 
	textplot(paste("PC",j,"\nCluster",thisCluster ,":",clusterSize, "proteins\n elimWilcoxMean:",
		format.pval(clusterSummaries$elimWilcoxMean[thisCluster ]),
		"\n elimWilcoxWeightedMean", format.pval(clusterSummaries$elimWilcoxWeightedMean[thisCluster]),
		"\n GeoMean:",format.pval(clusterSummaries$elimWilcoxGeoMean[thisCluster]), "\n Cluster Enrichment: ",clusterSummaries$clusterEnrichment[thisCluster] ,
		sep=" "),cex=0.8,halign="left")
	proteinIds <- proteinByLiverSample.ubiquitous$spId[match(seedList[[thisCluster]],proteinByLiverSample.ubiquitous$spAccession)]
	outText <- ""
	for(l in 1:ceiling(length(proteinIds)/5)) {
		loopIndex <- ((l-1) *5) + 1:5
		thisLine <- paste(na.omit(proteinIds[loopIndex]),collapse=",")
		outText <- paste(outText,thisLine ,"\n")
	}

	textplot(outText)	
	textplot(thisTermResults,mar=c(1,1,1,1))
}
dev.off()

outputFile <- paste("SummaryTable_PC",j,".pdf",sep="")
pdf(file=outputFile , width=15,height=15)
outputResults <- combResults[order(combResults$elimWilcox),]
par(cex=1)
layout(matrix(c(1,2), byrow=T),heights=c(1,5)) 
title<- paste("PC",j,"\n Top 50 results. \n Sorted by elimWilcox")
textplot(title)
textplot(outputResults[1:50,])
dev.off()

}  # end of output loop

} # end of j loop

# Bloody hell! That took some doing.
#combResults.pc1$goTerm 



## or could use a weight average p-value

## need number of proteins IN THIS FUNCTIONAL GROUP annotated to the GO term. 



#also need to output non clustered GO terms.






