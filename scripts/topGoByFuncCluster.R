
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################

########## FUNCTIONS

#### list proteins for a GO term in a precalculated goData object, accounting for ontology
listProtsInGo <- function(goTerm,ontology)  {
	switch(ontology,
				"BP" = genesInTerm(GOdata.BP,goTerm)[[1]],
				"MF" = genesInTerm(GOdata.MF,goTerm)[[1]],
				"CC" = genesInTerm(GOdata.CC,goTerm)[[1]]
				)
}

listProtsInGoFromList <- function(goTerm,ontology,goDataList)  {
	genesInTerm(goDataList[[ontology]],goTerm)[[1]]
}




# the topGOdata objects require a method for gene/protein selection. For enrichment analysis, we keep all proteins.
topDiffGenes <- function(allScore) {
	return(allScore )
}



listBestOverlappingCluster <- function(goProteins.valid, seedList, goContainProp)  { 
	#length(goProteins.valid)
	#lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid)))
	# which clusters contain greater than 'goContainProp' of the proteins annotated to this goTerm.
	maxValue <-  max(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))))
	clusHitList <- which(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))) == maxValue )
	if((length(clusHitList) > 0) & (maxValue > goContainProp)) {
		clusHitList 
	} else {NA}
}

listAllOverlappingClusters <- function(goProteins.valid, seedList, goContainProp)  { 
	#length(goProteins.valid)
	#lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid)))
	# which clusters contain greater than 'goContainProp' of the proteins annotated to this goTerm.
	clusHitList <- which(lapply(seedList,FUN = function(x) (length(intersect(x,goProteins.valid))/length(goProteins.valid))) > goContainProp)
	if(length(clusHitList) > 0) {
		clusHitList <- paste(clusHitList,collapse=",")
	} else {NA}
}


########### PROCESS


if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

###################INFO: Re-run GO analyses grouping GO terms by functional clusters.


###INFO:  Go through func clusters and pull in topGO results into a single table across BP,MF,CC


elimCutOff <- 0.01


###INFO: set up a topGOdata object for each ontology based on the list of interesting proteins -> goDataCollection.ubiq

geneList.ubiq <- ubi.pca.5.scores$Comp.1			#take an aribitrary set of scores.
names(geneList.ubiq) <- ubi.pca.5.scores$spAccession	# the names are important.

goDataCollection.ubiq <- list()
for(thisGOgraph in goGraphs )  {
	goDataCollection.ubiq[[thisGOgraph]] <- new("topGOdata",
		description =  "Liver proteins data set",
		ontology = thisGOgraph,
		allGenes = geneList.ubiq,
  		geneSelectionFun = topDiffGenes,
		nodeSize = nodeSizeValue ,
		annot = annFUN.GO2genes,
		GO2genes=go2prot 
		)
}


#GOdata.BP <- new("topGOdata",
#		description =  "Liver proteins data set",
#		ontology = "BP" ,
#		allGenes = geneList,
 # 		geneSelectionFun = topDiffGenes,
#		nodeSize = nodeSizeValue ,
#		annot = annFUN.GO2genes,
#		GO2genes=go2prot 
#		)
#
#GOdata.MF <- new("topGOdata",
#		description =  "Liver proteins data set",
#		ontology = "MF" ,
#		allGenes = geneList,
#		geneSelectionFun = topDiffGenes,
#		nodeSize = nodeSizeValue ,
#		annot = annFUN.GO2genes,
#		GO2genes=go2prot 
#		)
#
#GOdata.CC <- new("topGOdata",
#		description =  "Liver proteins data set",
#		ontology = "CC" ,
#		allGenes = geneList,
#		geneSelectionFun = topDiffGenes,
#		nodeSize = nodeSizeValue ,
#		annot = annFUN.GO2genes,
#		GO2genes=go2prot 
#		)
#


## following structure iterates through PCs and within each through GO ontologies (BP, MF, CC)
## runs multiple tests for GO enrichment of PC scores and stores p-values per GO term
## aim is to create one table for each PC with the results from the ontologies combined. 

summaryResultsList <- list()	# the PCs are treated separately. Their individual results are stored in elements of a list. 
#PC1
i<-1
for( i in 1:5)  {
	compHead <- paste("Comp.",i,sep="")		# this PC

	#combResults.pc1 <-data.frame()
	summaryResultsList[[i]] <- data.frame()	#initiate the data frame for all results from this PC

	for(thisGOgraph in goGraphs )  {		# loop through 3 GO ontologies

		geneList <- ubi.pca.5.scores[,compHead]
		names(geneList) <- ubi.pca.5.scores$spAccession
		
		# load the correct pre-computed GO data object (different trees based on ontology.
		GOdata <- switch(thisGOgraph,
				"BP" = GOdata.BP,
				"MF" = GOdata.MF,
				"CC" = GOdata.CC
				)

		GOdata <- updateGenes(GOdata,geneList,topDiffGenes)

		####INFO: define test statistics and apply over all GO terms 
		#INFO: resultWilcox = 2-sided wilcox test. Good for outliers in one direction
		test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
		resultWilcox <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcox = elim version of resultWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcox <- getSigGroups(GOdata, test.stat)

		#resultElimWilcoxAdj <- resultElimWilcox 
		#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue
		#INFO: resultKS = Kolmogorov smirnov test: may be good for general distribution changes.
		test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
		resultKS <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimKS = elim version of resultKS
		test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
		resultElimKS <- getSigGroups(GOdata, test.stat)

		#resultElimKSAdj <- resultElimKS
		#score(resultElimKSAdj ) <- qvalue(score(resultElimKS))$qvalue

		#INFO: change scores to absolute values for Greater than test 
		geneList <- abs(ubi.pca.5.scores[,compHead])
		names(geneList) <- ubi.pca.5.scores$spAccession
		GOdata <- updateGenes(GOdata,geneList,topDiffGenes)
		#INFO: resultWilcoxAbs = test for outliers in both directions simultaneously. Less power than wilcox if true difference is unidireectional.
		test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
		resultWilcoxAbs <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcoxAbs = elim version of resultElimWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)

		#resultElimWilcoxAbsAdj <- resultElimWilcoxAbs 
		#score(resultElimWilcoxAbsAdj ) <- qvalue(score(resultElimWilcoxAbs ))$qvalue


		#INFO: collect and bind all the results together. 
		# This could be done with a function available in topGO but I wanted extra columns and control.
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
		
		#INFO: bind the results from this GO ontology as rows to the table for the current PC in 'summaryResultsList'
		summaryResultsList[[i]] <- rbind(summaryResultsList[[i]],thisGoResults)

	} # end of GO loops
} # end of PC loops


#### Use seedList and go through output table.  

#combResults.pc1
#seedList




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



#######INFO:  try condensing a GO results table based on functional clusters. 
#INFO:  Join GO terms into a group if the functional cluster contains a high proportion of the  proteins linked to that GO term
goContainProp <- 0.8  # same result with 0.7

# protein func clusters: 
# seedList
# better check how long this is. 

# go detected results table:
# summaryDetectResultsList 

#INFO:  ?need to order the results (and start at the top?)
# perhaps filter by significance cut-off
orderBy <- "elimFisher"
# or could do overlap for all GO terms? Approx 5000


#listProtsInGo(goTerm,ontology)


goTerm <- summaryDetectResultsList$goTerm[1]

validProts <- proteinByLiverSample.ubiquitous$spAccession
#goProteins <- na.omit(listProtsInGo(summaryDetectResultsList$goTerm[1],summaryDetectResultsList$ontology[1]))
goProteins <- na.omit(listProtsInGo("GO:0006414","BP"))
goProteins <- na.omit(listProtsInGo("GO:0005759","CC"))	# mito matrix
goProteins <- na.omit(listProtsInGo("GO:0005743","CC")) 	# mito inner membrane
goProteins <- na.omit(listProtsInGo("GO:0005739","CC"))	# mitochondrion
goProteins.valid <- intersect(validProts,goProteins)


 
goProteins.valid <- intersect(validProts,goProteins)



####INFO: First call clusters for detected/undetected GO table
#INFO: cannot use topGOdata objects from PCs as they removed nodes (GO terms with less than xxx proteins).


#listOverlappingClusters(goProteins.valid, seedList, goContainProp)
#apply(summaryDetectResultsList[1:5,],1, FUN= function(x) intersect(listProtsInGo(goTerm=x["goTerm"],ontology=x["ontology"]),validProts))
summaryDetectResultsList$allClusters  <- apply(summaryDetectResultsList,1, FUN= function(x) listAllOverlappingClusters(
						goProteins.valid=intersect(listProtsInGo(goTerm=x["goTerm"],ontology=x["ontology"]),validProts),
						seedList	,goContainProp ))

#summaryDetectResultsList <- subset(summaryDetectResultsList, select=-bestCluster)

summaryDetectResultsList$bestCluster  <- apply(summaryDetectResultsList,1, FUN= function(x) listBestOverlappingCluster(
						goProteins.valid=intersect(listProtsInGo(goTerm=x["goTerm"],ontology=x["ontology"]),validProts),
						seedList	,goContainProp ))



write.table(summaryDetectResultsList[order(summaryDetectResultsList$elimFisher),],file="detectGOWithClusterNumberElimFisher.tab",sep="\t",quote=F,row.names=F)


# format(seedList)

# find a way to output the table (or a subsection) with (significant) cluster members grouped. 
# Only works with single (best) cluster per GO term
detectTableSigThreshold <- 1.0e-05
sigTable <- subset(summaryDetectResultsList,summaryDetectResultsList[,orderBy] < detectTableSigThreshold)
#summaryDetectResultsList[order(summaryDetectResultsList[,orderBy]),]
sigTable <- sigTable[order(sigTable[,orderBy]),]
sigTable$outputRank <- rank(sigTable[,orderBy])

#levels(as.factor(as.character(sigTable$bestCluster)))

sigTable$bestCluster <- as.factor(as.character(sigTable$bestCluster))
for(thisCluster in levels(sigTable$bestCluster))  {
	if(thisCluster == "NA") {break}
	rankToSet <- min(sigTable[sigTable$bestCluster == thisCluster,"outputRank"])
	sigTable[sigTable$bestCluster == thisCluster,"outputRank"] <- rankToSet
}	
sigTable <- sigTable[order(sigTable$outputRank),]
write.table(sigTable,file="sigDetectGoTableByCluster.tab",quote=F,row.names=F,sep="\t")



## once the top terms have been defined, replace GO terms with clusters of proteins where appropriate. 

####INFO: Second, call clusters for PCs GO table
############## FOR PC TABLES ONLY!  - BELOW IS USING DETECTED RESULTS ON PC TABLES!


#INFO: Assign best functional cluster to GO terms, if applicable.
###PC table...

#summmaryPcResultList[[i]]


for(i in 1:length(summmaryPcResultList)) {
	
	summmaryPcResultList[[i]]$bestCluster  <- apply(summmaryPcResultList[[i]],1, FUN= function(x) listBestOverlappingCluster(
						goProteins.valid=intersect(listProtsInGo(goTerm=x["goTerm"],ontology=x["ontology"]),validProts),
						seedList	,goContainProp ))
}



#############INFO: GRAPHING PC SCORES FROM TOP SIG TERMS.


#cluster annotations. Cluster 4 (seedList[[4]]) are ribosomal proteins.
clusterAnnotations <- character(length(seedList))
names(clusterAnnotations) <- paste("CLUSTER", 1:length(seedList))

#clusterAnnotations[["CLUSTER 1"]]


# find a way to output the table (or a subsection) with (significant) cluster members grouped. 
# Only works with single (best) cluster per GO term
orderBy.PC <- "elimWilcox" ;pcTableSigThreshold <- 1.0e-04
# orderBy.PC <- "elimWilcox" ;pcTableSigThreshold <- 1.0e-05	
#orderBy.PC <- "elimAbsWilcox" ; pcTableSigThreshold <- 1.0e-03			# gives nice set of results
#orderBy.PC <- "elimKS" ; pcTableSigThreshold <- 1.0e-05
pcTableSigThreshold <- 1.0e-04
## how many groups to plot? A fixed number or a cut-off?
numbGraphResults <- 15



addGroupScoresToGroupTable <- function(groupTable,groupIds,groupName,scoreTable,idColumn,scoreColum) {
	index <- na.omit(match(groupIds , scoreTable[,idColumn]))
	scores <- scoreTable[index ,scoreColum]
	groupTable <- rbind(groupTable,data.frame(score=scores,group=groupName))
	
}



#pc.i <- 1
#INFO: Create violin plots
library(lattice)

bwPlotsAsPdf <- function(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15)  {
	pdfName <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,"pdf",sep=".")

	pdf(pdfName, paper="a4")
	frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results",sep=" ")
	textplot(frontPageText)
	for(pc.i in 1:length(summmaryPcResultList)) {
		plotPCsFromSummary(orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
	}
	dev.off()
}

bwPlotsAsFigs <- function(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, figType="tiff")  {
	#frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results",sep=" ")
	#textplot(frontPageText)
	#op <- par(cex=1.5)
	for(pc.i in 1:length(summmaryPcResultList)) {
		plotDir <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep="_")
		if(!file.exists(plotDir))  {dir.create(plotDir)}
		fileName <- paste("vioPlotsByPC",pc.i,orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep=".")
		fileName <- paste(plotDir,fileName,sep="/")
		tiff(fileName, compression="lzw",width=180, height=480,units="mm",res=300)	
		plotPCsFromSummary(orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
		dev.off()	
	}
	#par(op)
}


plotPCsFromSummary <- function(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, pc.i=1)  {

	sigTable <- subset(summmaryPcResultList[[pc.i]],summmaryPcResultList[[pc.i]][,orderBy.PC] < pcTableSigThreshold)
	sigTable <- sigTable[order(sigTable[,orderBy.PC]),]
	sigTable$outputRank <- rank(sigTable[,orderBy.PC])

	sigClusters <- levels(as.factor(as.character(sigTable$bestCluster)))
	for(thisCluster in sigClusters)  {
		if(is.na(thisCluster)) {break}
		rankToSet <- min(na.omit(sigTable[sigTable$bestCluster == thisCluster,"outputRank"]))
		sigTable[which(sigTable$bestCluster == thisCluster),"outputRank"] <- rankToSet
	}	
	sigTable <- sigTable[order(sigTable$outputRank),]

	## need to count each cluster only once and each independent go term once. 
	
	plotTable <- sigTable[match(unique(sigTable$outputRank),sigTable$outputRank),]
	numbGraphResultsLocal <- min(nrow(plotTable),numbGraphResults)
	if(numbGraphResultsLocal < 1)  {
		warningMess <- paste("No sig results for PC", pc.i,"\n using", orderBy.PC, "below", pcTableSigThreshold)
		textplot(warningMess)
		return(NULL)
	}
	plotTable <- plotTable[1:numbGraphResultsLocal,]

	# could try strwrp or substr
	compHead <- paste("Comp.",pc.i,sep="")
	baseScore <- data.frame(score=ubi.pca.5.scores[,compHead],group="                            All Proteins")
	groupTable <- baseScore
	plotList <- list()
	for(i in 1:nrow(plotTable))  {
		if(is.na(plotTable$bestCluster[i]))  {
			plotList[[i]] <- intersect(listProtsInGo(goTerm=plotTable[i,"goTerm"],ontology=plotTable[i,"ontology"]),validProts)
			goID <- paste(paste(strwrap(plotTable[i,"description"],width=40),collapse="\n"),plotTable[i,"goTerm"],sep="\n")
		}	else {
			plotList[[i]] <- seedList[[plotTable$bestCluster[i]]]
			goID <- paste("CLUSTER", plotTable$bestCluster[i])
		}
		groupTable <- addGroupScoresToGroupTable(groupTable=groupTable,groupIds=plotList[[i]],
					groupName=goID,scoreTable=ubi.pca.5.scores,
					idColumn="spAccession",scoreColum=compHead)
	}


	bwTitle <- compHead
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	bwTitle <- compHead
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	plot(bwplot(bymedian ~ score, groupTable, main=bwTitle,  aspect="xy", cex.main=2, cex.lab=2,
		par.settings = list(layout.widths = list(axis.left = 0, ylab.axis.padding = 40)),
	       panel = function(..., box.ratio) {
		   panel.violin(..., col = "darkgrey",
				varwidth = FALSE, box.ratio = box.ratio)
		   panel.bwplot(..., fill = "lightgrey", box.ratio = .1)
	       } )
	)
}



bwPlotsAsPdf(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10)
bwPlotsAsPdf(orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10)
bwPlotsAsPdf(orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10)

bwPlotsAsPdf(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-05, numbGraphResults =11)


bwPlotsAsFigs(orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10)
bwPlotsAsFigs(orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10)
bwPlotsAsFigs(orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10)




stopifnot(FALSE)

######  old gubbins

