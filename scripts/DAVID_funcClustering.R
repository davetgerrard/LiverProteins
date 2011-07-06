
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

##DAVID Clustering: A Heuristic Multiple Linkage Fuzzy Clustering Procedure

makeNew <- FALSE	# set to TRUE if want to re-calc binary matrix (add 30 minutes to run time)


###temp: get list of protein/GO annotaitons

protTerms <- prot2go[c(proteinByLiverSample.ubiquitous$spAccession)]

fullProtList <- proteinByLiverSample.ubiquitous$spAccession
fullTermList <- unique(as.character(unlist(protTerms)))

if(makeNew) {		# only create binaryGrid and kappaMatrix if none available

	### require binary grid of genes vs terms 

	binaryGrid <- matrix(0,nrow=length(fullProtList),ncol=length(fullTermList),dimnames=list(fullProtList,fullTermList))

	##binaryGrid["P61604",as.character(unlist(prot2go["P61604"]))] <- 1

	##possible to do this without FOR loop?
	for(i in 1:length(fullProtList))  {
		binaryGrid[fullProtList[i],as.character(unlist(prot2go[fullProtList[i]]))] <- 1
	}


	### create kappa matrix for genes against genes based on co-occurrence of terms
	## can eliminate genes with <4 (minTermSize) terms. 

	#number of terms per protein
	termCounts <- rowSums(binaryGrid)

	#contTable <- table(binaryGrid[37,],binaryGrid[38,])[c(2:1),c(2:1)]
	#contTable <- table(binaryGrid[1,],binaryGrid[2,])[c(2:1),c(2:1)]
	#contTable <- table(binaryGrid[1,],binaryGrid[1,])[c(2:1),c(2:1)]

	davidKappaFromTable <- function(contTable)  {

		obs <- (contTable[1,1] + contTable[2,2]) / sum(contTable)
		exp <- (rowSums(contTable)[1] * colSums(contTable)[1] + rowSums(contTable)[2] * colSums(contTable)[2])/ sum(contTable)^2
		davidKappa <- (obs - exp) / (1-exp)
	}
	(as.numeric(davidKappaFromTable(contTable)))


	#### THIS IS SLOW! About 30 mins. have saved the table
	kappaMatrix <-matrix(0,nrow=nrow(binaryGrid),ncol=nrow(binaryGrid),dimnames=list( row.names(binaryGrid), row.names(binaryGrid)))
	for(i in 1:nrow(kappaMatrix))  {
		for(j in 1:nrow(kappaMatrix)) {
			contTable <- table(binaryGrid[i,],binaryGrid[j,])[c(2:1),c(2:1)]
			kappaMatrix[i,j] <- as.numeric(davidKappaFromTable(contTable))
		}
	}


	#write.table(kappaMatrix, file="ubiProtsKappaMatrix.tab",quote=F,sep="\t")
} else {
	# use pre-computed kappaMatrix
	kappaMatrix <- read.delim("C:/Users/dave/LiverProteins/data/ubiProtsKappaMatrixDetected.tab",sep="\t",header=T)

}

### creates seed from each gene of clusters of genes with greater than 0.35 (start_kappa)
### all relationships within cluster must be greater than start_kappa

startKappa <- 0.65	# 0.65 is stronger threshold for kappaMatrix from deteceted protiens
propThreshold <- 0.5
minClusterSize <- 5

seedList <- list()
listCounter <- 1
for(i in 1:nrow(kappaMatrix)) {
	clusterIndex <- as.numeric(which(kappaMatrix[,i] > startKappa))
	#kappaMatrix[clusterIndex,clusterIndex]
	# remove members not passing startKappa threshold for propThreshold of pairs
	if(length(clusterIndex) < minClusterSize) {next;}
	# fliter out proteins that don't show kappa values above startKappa for >50% others in group.	
	clusterIndex.prop <- clusterIndex[(colSums(kappaMatrix[clusterIndex,clusterIndex] > startKappa) - 1)  > floor(length(clusterIndex) * propThreshold )]
	if(length(clusterIndex.prop) >= minClusterSize  )  {
		#add this seedCluster to the list.
		seedList[[listCounter]] <- row.names(kappaMatrix)[clusterIndex]
		listCounter <- listCounter+1
	}
}

#testList <- seedList

### iteratively merge clusters sharing more than 50% (cluster_similarity) of members

# start at one. if another cluster overlaps by 50%, merges these to a new list and remove from list. 
# Should just one match proportion or both?    
# If just one, then smaller groups will get rolled up into larger groups.
# If both, then fewer groups will get made.
propToMerge <- 0.5


## returns and index of the first pair of clusters to merge or 1 if no suitable pairs.
getVectorToMerge <- function(seedList)  {
	for(i in 2:length(seedList)) {
		minLength <- min(length(seedList[[1]]),length(seedList[[i]]))	# use this for single group limit
		#minLength <- max(length(seedList[[1]]),length(seedList[[i]]))	# use this for double group limit		
		if(length(intersect(seedList[[1]],seedList[[i]])) > floor(minLength * propToMerge ))  {
			return(c(1,i))
		}
	}
	return(1) 
}

#seedList <- newList
#seedList <- testList

mergeList <- function(seedList)  {		# one iteration to merge components of a term list.
	newList <- list()
	j <- 1
	while(length(seedList) > 1)  {
		mergeVector <- getVectorToMerge(seedList)
		if(length(mergeVector) > 1)  {
			newList[[j]] <- union(seedList[[mergeVector[1]]],seedList[[mergeVector[2]]])
			seedList[[mergeVector[2]]] <- NULL  # must be in reverse order for index to work!
			seedList[[mergeVector[1]]] <- NULL

		} else  {		# no suitable group to merge, pass the vector unchanged.
			newList[[j]] <- seedList[[mergeVector[1]]]
			seedList[[mergeVector[1]]] <- NULL
		}
		j <- j+1
	}
	if(length(seedList) > 0)  {
		newList[[j]] <- seedList[[1]]
		seedList[[1]] <- NULL
	}

	return(newList)

}

startLength <- 1
endLength <- 0

while(startLength > endLength)  {
	startLength <- length(seedList)	
	seedList <- mergeList(seedList)
	endLength <- length(seedList)	
}


#length(seedList)
### how is the cluster turned into a cluster of terms?



stopifnot(FALSE)




# GO terms for a cluster of proteins
#go.terms <- unique(as.character(unlist(prot2go[seedList[[1]]])))

#check if any of the functional groups look like the ribosomal proteins annotated as pancreas.
#lapply(seedList,FUN = function(x) intersect(x,panc.genes))
#go.terms <- unique(as.character(unlist(prot2go[seedList[[5]]])))


#tabulate common GO terms within a func group.
#table(as.character(unlist(prot2go[seedList[[5]]])))
#sorted table (1 to 8)
#sort(table(as.character(unlist(prot2go[seedList[[5]]]))),decreasing=T )[1:8]
#all terms shared by half or more.

#propTermShare <- 0.5
#topTerms <- names(which(sort(table(as.character(unlist(prot2go[seedList[[5]]]))),decreasing=T ) >= ((length(seedList[[5]])*propTermShare) )))

#xx[topTerms]


#### could run KS/wilcox tests on protein clusters or combine p-values of 'functional categories'



