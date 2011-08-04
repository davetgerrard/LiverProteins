
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################


########### PRE-REQUISITES


# for this, a table of detected proteins.
# source("C:/Users/dave/LiverProteins/scripts/loadPCA_noOutput.R")

# A map of protein to GO annotations. 

# for the funcitonal clustering, the proteins need to have been clustered into functional groups based on shared annotations.
# e.g. source("C:/Users/dave/LiverProteins/scripts/DAVID_funcClustering.R")
# load utilities for performing functional clustering of GO result tables.
#source("C:/Users/dave/LiverProteins/scripts/funcClusterGoTable.R")
# load utilities for running multiple GO analyses.
#source("C:/Users/dave/LiverProteins/scripts/topGoLiverProteinsUtilities.R")


########### PROCESS


if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

###################INFO: Re-run GO analyses grouping GO terms by functional clusters.

###INFO:  Go through func clusters and pull in topGO results into a single table across BP,MF,CC

# standard topGO settings
elimCutOff <- 0.01
pcs <- 1:4
goGraphs <- c("BP","MF","CC")
nodeSizeValue <- 10	
topNodesValue <- 20
topTerms <- 50

###INFO: set up a topGOdata object for each ontology based on the list of interesting proteins -> goDataCollection.ubiq
## Used in performing GO tests and to retrieve appropriate protein lists from GO terms. 


geneList.ubiq <- ubi.pca.5.scores$Comp.1			#take an aribitrary set of scores.
names(geneList.ubiq) <- ubi.pca.5.scores$spAccession	# the names are important.

goDataCollection.ubiq.scored <- list()
for(thisGOgraph in goGraphs )  {
	goDataCollection.ubiq.scored[[thisGOgraph]] <- new("topGOdata",
		description =  "Liver proteins data set",
		ontology = thisGOgraph,
		allGenes = geneList.ubiq,
  		geneSelectionFun = topDiffGenes,
		nodeSize = nodeSizeValue ,
		annot = annFUN.GO2genes,
		GO2genes=go2prot 
		)
}


## following structure iterates through PCs and within each through GO ontologies (BP, MF, CC)
## runs multiple tests for GO enrichment of PC scores and stores p-values per GO term
## aim is to create one table for each PC with the results from the ontologies combined. 

summaryScoreResults <- list()	# the PCs are treated separately. Their individual results are stored in elements of a list. 

for( i in 1:5)  {
	compHead <- paste("Comp.",i,sep="")		# this PC
	geneList <- ubi.pca.5.scores[,compHead]	# get PC scores for proteins
	names(geneList) <- ubi.pca.5.scores$spAccession		# make sure scores have names of proteins

	summaryScoreResults[[i]] <- data.frame()	#initiate the data frame for all results from this PC

	for(thisGOgraph in goGraphs )  {		# loop through 3 GO ontologies

		
		goDataCollection.ubiq.scored[[thisGOgraph]] <- updateGenes(goDataCollection.ubiq.scored[[thisGOgraph]],geneList,topDiffGenes)		
		
		summaryScoreResults[[i]] <- rbind(summaryScoreResults[[i]],runScoreGoTests(goDataCollection.ubiq.scored[[thisGOgraph]],geneList,topDiffGenes))

	}
}

## output the plain GO on PC results
for(i in 1:length(summaryScoreResults)) {
	outFileName <- paste("GoSummaryByPc",i,"tab", sep=".")
	write.table(summaryScoreResults[[i]],file=outFileName,sep="\t",quote=F,row.names=F) 
	rm(outFileName)
}


#INFO:  Join GO terms into a group if the functional cluster contains a high proportion of the  proteins linked to that GO term
goContainProp <- 0.8  # same result with 0.7
detectTableSigThreshold <- 1.0e-03

#INFO:  ?need to order the results (and start at the top?)
# perhaps filter by significance cut-off
orderBy <- "elimWilcox"
#INFO: Assign best functional cluster to GO terms, if applicable.
###PC table...

#summmaryPcResultList[[i]]


for(i in 1:length(summaryScoreResults)) {
	
	summaryScoreResults[[i]]$bestCluster  <- apply(summaryScoreResults[[i]],1, FUN= function(x) listBestOverlappingCluster(
						goProteins.valid=intersect(listProtsInGoFromList(goTerm=x["goTerm"],ontology=x["ontology"],goDataCollection.ubiq.scored),validProts),
						seedList	,goContainProp ))
	outFileName <- paste("scoreGoWithClusterByPc",i,"tab", sep=".")
	write.table(summaryScoreResults[[i]][order(summaryScoreResults[[i]][,orderBy]),],file=outFileName,sep="\t",quote=F,row.names=F) 	
	rm(outFileName)
}


# use same detectTableSigThreshold for each PC table?



summaryScoreResults.sigClustered <- list()
for(i in 1:length(summaryScoreResults)) {
	summaryScoreResults.sigClustered[[i]] <- clusterTable(summaryScoreResults[[i]],orderBy=orderBy,detectTableSigThreshold=detectTableSigThreshold)     # function(table,orderBy,detectTableSigThreshold = 1.0e-05)  {
	outFileName <- paste("GoSummaryTopClusteredByPc",i,"tab", sep=".")
	write.table(summaryScoreResults.sigClustered[[i]],file=outFileName,sep="\t",quote=F,row.names=F) 	
	rm(outFileName)
}




