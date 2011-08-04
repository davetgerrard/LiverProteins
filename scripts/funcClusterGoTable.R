
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################


####INFO: utility functions to combine a GO results table with functional clusters.


############GENERIC

# give the index of the functional cluster which best represents this go term IF the cluster contains at least goContainProp of the proteins assigned to the go term.
listBestOverlappingCluster <- function(goProteins.valid, seedList, goContainProp)  { 
	# which clusters contain greater than 'goContainProp' of the proteins annotated to this goTerm.
	maxValue <-  max(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))))
	clusHitList <- which(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))) == maxValue )
	if((length(clusHitList) > 0) & (maxValue > goContainProp)) {
		clusHitList[1] 
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

### filter a table based on a treshold and reorders by keeping members of the same cluster together. 
## the table is filtered, marked and re-ordered, grouping GO terms which represent the same functional cluster of proteins/genes
clusterTable <- function(goTable,orderBy,detectTableSigThreshold = 0.05)  {
	# find a way to output the table (or a subsection) with (significant) cluster members grouped. 
	# Only works with single (best) cluster per GO term. Limit to terms below significance threshold.
	sigTable <- subset(goTable,goTable[,orderBy] < detectTableSigThreshold)
	# order the table to cluster lower terms with top term of each cluster.
	sigTable <- sigTable[order(sigTable[,orderBy]),]
	sigTable$outputRank <- sigTable$pureRank<- rank(sigTable[,orderBy])
	#levels(as.factor(as.character(sigTable$bestCluster)))
	# assign same rank to all significant members of same cluster.
	sigTable$bestCluster <- as.factor(as.character(sigTable$bestCluster))
	for(thisCluster in levels(sigTable$bestCluster))  {
		if(is.na(thisCluster)) {break}
		rankToSet <- min(na.omit(sigTable$outputRank[sigTable$bestCluster == thisCluster]))
		sigTable$outputRank[sigTable$bestCluster == thisCluster] <- rankToSet
	}	
	sigTable$subClusterTerm <- ifelse((sigTable$outputRank == sigTable$pureRank),FALSE,TRUE)
	sigTable <- sigTable[order(sigTable$outputRank),]
}










