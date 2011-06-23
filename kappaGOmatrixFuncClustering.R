
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

## this is clustering based on GO terms not proteins. 
## not what DAVID grouping suggests but might make sense


####check which binaryGrid was used.

binaryGrid.GO.detected <- t(binaryGrid.detected)

kappaMatrix.GO.detected <-matrix(0,nrow=nrow(binaryGrid.GO.detected),ncol=nrow(binaryGrid.GO.detected),
		dimnames=list( row.names(binaryGrid.GO.detected), row.names(binaryGrid.GO.detected)))
for(i in 1:nrow(kappaMatrix.GO.detected))  {
	for(j in 1:nrow(kappaMatrix.GO.detected)) {
		contTable <- table(binaryGrid.GO.detected[i,],binaryGrid.GO.detected[j,])[c(2:1),c(2:1)]
		kappaMatrix.GO.detected[i,j] <- as.numeric(davidKappaFromTable(contTable))
	}
}
write.table(kappaMatrix.GO.detected, file="ubiProtsKappaMatrixGO.Detected.tab",quote=F,sep="\t")



kappaMatrix.GO.detected <- read.delim("C:/Users/dave/LiverProteins/data/ubiProtsKappaMatrixGO.Detected.tab",header=T)


### creates seed from each gene of clusters of genes with greater than 0.35 (start_kappa)
### all relationships within cluster must be greater than start_kappa

startKappa.GO <- 0.65	# 0.65 is stronger threshold for kappaMatrix from deteceted protiens
propThreshold.GO <- 0.5
minClusterSize.GO <- 5

seedList.GO <- list()
listCounter <- 1
for(i in 1:nrow(kappaMatrix.GO.detected)) {
	clusterIndex <- as.numeric(which(kappaMatrix.GO.detected[,i] > startKappa.GO))
	#kappaMatrix.GO.detected[clusterIndex,clusterIndex]
	# remove members not passing startKappa threshold for propThreshold of pairs
	if(length(clusterIndex) < minClusterSize) {next;}
	# fliter out terms that don't show kappa values above startKappa for >50% others in group.	
	clusterIndex.prop <- clusterIndex[(colSums(kappaMatrix.GO.detected[clusterIndex,clusterIndex] > startKappa.GO) - 1)  > floor(length(clusterIndex) * propThreshold )]
	if(length(clusterIndex.prop) >= minClusterSize  )  {
		#add this seedCluster to the list.
		seedList.GO[[listCounter]] <- row.names(kappaMatrix.GO.detected)[clusterIndex]
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
	startLength <- length(seedList.GO)	
	seedList <- mergeList(seedList.GO)
	endLength <- length(seedList.GO)	
}

### what was the point of all this?

