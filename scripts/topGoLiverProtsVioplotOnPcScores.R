
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################

########## FUNCTIONS

addGroupScoresToGroupTable <- function(groupTable,groupIds,groupName,scoreTable,idColumn,scoreColum) {	#INFO: Gets data for a list of proteins and binds as rows onto a table. Resulting stack can be used as basis for large vioplot or boxplot of several groups.
	index <- na.omit(match(groupIds , scoreTable[,idColumn]))
	scores <- scoreTable[index ,scoreColum]
	groupTable <- rbind(groupTable,data.frame(score=scores,group=groupName))
	
}


bwPlotsAsPdf <- function(resultTableList,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, goDataList=goDataCollection.ubiq.scored)  {	#INFO: Control plotPCsFromSummary() for output in a pdf
	pdfName <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,"pdf",sep=".")

	pdf(pdfName, paper="a4")
	frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results",sep=" ")
	textplot(frontPageText)
	for(pc.i in 1:length(resultTableList)) {
		plotPCsFromSummary(resultTableList[[pc.i]],orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
	}
	dev.off()
}

bwPlotsAsFigs <- function(resultTableList,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, figType="tiff", goDataList=goDataCollection.ubiq.scored)  { #INFO: Control plotPCsFromSummary() for output in a figure
	for(pc.i in 1:length(resultTableList)) {
		plotDir <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep="_")
		if(!file.exists(plotDir))  {dir.create(plotDir)}
		fileName <- paste("vioPlotsByPC",pc.i,orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep=".")
		fileName <- paste(plotDir,fileName,sep="/")
		tiff(fileName, compression="lzw",width=180, height=480,units="mm",res=300)	
		plotPCsFromSummary(resultTableList[[pc.i]],orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i, goDataList=goDataList)
		dev.off()	
	}
	#par(op)
}




## Currently very specific to this project and data objects. Could generalise by specifying scoreTable and score column as variables.
plotPCsFromSummary <- function(resultsTable,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, pc.i=1, goDataList=goDataCollection.ubiq.scored)  {	#INFO: gets data and plots a nice set of violin plots
	#INFO: need to make a cluster table out of the resultsTable so that redundant terms are shown as a cluster.
	sigTable <- clusterTable(resultsTable, orderBy=orderBy.PC,detectTableSigThreshold = pcTableSigThreshold) 

	## need to count each cluster only once and each independent go term once. 
	# take the first instance of each outputRank
	plotTable <- sigTable[match(unique(sigTable$outputRank),sigTable$outputRank),]
	numbGraphResultsLocal <- min(nrow(plotTable),numbGraphResults)
	if(numbGraphResultsLocal < 1)  {	# are there going to be enough results to plot?
		warningMess <- paste("No sig results for PC", pc.i,"\n using", orderBy.PC, "below", pcTableSigThreshold)
		textplot(warningMess)
		return(NULL)
	}
	plotTable <- plotTable[1:numbGraphResultsLocal,]	# only going to try to plot the top terms

	# could try strwrp or substr
	compHead <- paste("Comp.",pc.i,sep="")
	#INFO: Want the scores for ALL proteins as first on plot for visual comparison.
	baseScore <- data.frame(score=ubi.pca.5.scores[,compHead],group="                            All Proteins")		# long name to handle padding of output
	groupTable <- baseScore
	plotList <- list()
	for(i in 1:nrow(plotTable))  {
		if(is.na(plotTable$bestCluster[i]))  {
			plotList[[i]] <- intersect(listProtsInGoFromList(goTerm=as.character(plotTable[i,"goTerm"]),ontology=as.character(plotTable[i,"ontology"]),goDataList),validProts)
			# some goDescriptions are very long. Currently wrapping the terms rather than truncating.
			goID <- paste(paste(strwrap(plotTable[i,"description"],width=40),collapse="\n"),plotTable[i,"goTerm"],sep="\n")
		}	else {
			# this should be a cluster
			plotList[[i]] <- seedList[[plotTable$bestCluster[i]]]
			goID <- paste("CLUSTER", plotTable$bestCluster[i])
		}
		# BEWARE: lots of dependencies here
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




##########INFO: PREREQUISITES


library(lattice)

#seedList, goDataCollection.ubiq.scored, ubi.pca.5.scores


##########INFO: PROCESS



####INFO: GRAPHING PC SCORES FROM TOP SIG TERMS.


#cluster annotations. Cluster 4 (seedList[[4]]) are ribosomal proteins.
clusterAnnotations <- character(length(seedList))
names(clusterAnnotations) <- paste("CLUSTER", 1:length(seedList))

#clusterAnnotations[["CLUSTER 1"]]


pcTableSigThreshold <- 1.0e-04
## how many groups to plot? A fixed number or a cut-off?
numbGraphResults <- 15

#INFO: Create violin plots as pdf and as figures (tiffs)

bwPlotsAsPdf(summaryScoreResults,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)
bwPlotsAsPdf(summaryScoreResults,orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)
bwPlotsAsPdf(summaryScoreResults,orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)

bwPlotsAsPdf(summaryScoreResults,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-05, numbGraphResults =11, goDataList=goDataCollection.ubiq.scored)

bwPlotsAsFigs(summaryScoreResults,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)
bwPlotsAsFigs(summaryScoreResults,orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)
bwPlotsAsFigs(summaryScoreResults,orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10, goDataList=goDataCollection.ubiq.scored)


