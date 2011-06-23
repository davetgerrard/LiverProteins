## funcClustering using only select GO terms shown by 10+ proteins



#protTerms <- prot2go[c(proteinByLiverSample.ubiquitous$spAccession)]

fullProtList <- proteinByLiverSample.ubiquitous$spAccession
#fullTermList <- unique(as.character(unlist(protTerms)))


##hack to get the full list of relevant GO terms (accounting for GO structures)

goTerms.BP <- usedGO(GOdata.BP)
goTerms.MF  <- usedGO(GOdata.MF)
goTerms.CC <- usedGO(GOdata.CC)

bg.BP <- matrix(0,nrow=length(fullProtList),ncol=length(goTerms.BP),dimnames=list(fullProtList,goTerms.BP))
for(i in 1:length(goTerms.BP)) {
	bg.BP[as.character(unlist(genesInTerm(GOdata.BP,goTerms.BP[i]))),goTerms.BP[i]] <- 1
}
bg.MF <- matrix(0,nrow=length(fullProtList),ncol=length(goTerms.MF),dimnames=list(fullProtList,goTerms.MF))
for(i in 1:length(goTerms.MF)) {
	bg.MF[as.character(unlist(genesInTerm(GOdata.MF,goTerms.MF[i]))),goTerms.MF[i]] <- 1
}
bg.CC <- matrix(0,nrow=length(fullProtList),ncol=length(goTerms.CC),dimnames=list(fullProtList,goTerms.CC))
for(i in 1:length(goTerms.CC)) {
	bg.CC[as.character(unlist(genesInTerm(GOdata.CC,goTerms.CC[i]))),goTerms.CC[i]] <- 1
}

binaryGrid.detected <- cbind(bg.BP,bg.MF,bg.CC)
allDetectGoTerms <- c(goTerms.BP,goTerms.MF,goTerms.CC)

##takes about 10 minutes
kappaMatrix.detected <-matrix(0,nrow=nrow(binaryGrid.detected),ncol=nrow(binaryGrid.detected),
		dimnames=list( row.names(binaryGrid.detected), row.names(binaryGrid.detected)))
for(i in 1:nrow(kappaMatrix.detected))  {
	for(j in 1:nrow(kappaMatrix.detected)) {
		contTable <- table(binaryGrid.detected[i,],binaryGrid.detected[j,])[c(2:1),c(2:1)]
		kappaMatrix.detected[i,j] <- as.numeric(davidKappaFromTable(contTable))
	}
}
write.table(kappaMatrix.detected, file="ubiProtsKappaMatrixDetected.tab",quote=F,sep="\t")


startKappa <- 0.65
propThreshold <- 0.5
minClusterSize <- 10

seedList <- list()
listCounter <- 1
for(i in 1:nrow(kappaMatrix.detected)) {
	clusterIndex <- as.numeric(which(kappaMatrix.detected[,i] > startKappa))
	#kappaMatrix.detected[clusterIndex,clusterIndex]
	# remove members not passing startKappa threshold for propThreshold of pairs
	if(length(clusterIndex) < minClusterSize) {next;}
	# fliter out proteins that don't show kappa values above startKappa for >50% others in group.	
	clusterIndex.prop <- clusterIndex[(colSums(kappaMatrix.detected[clusterIndex,clusterIndex] > startKappa) - 1)  > floor(length(clusterIndex) * propThreshold )]
	if(length(clusterIndex.prop) >= minClusterSize  )  {
		#add this seedCluster to the list.
		seedList[[listCounter]] <- row.names(kappaMatrix.detected)[clusterIndex]
		listCounter <- listCounter+1
	}
}


### iteratively merge clusters sharing more than 50% (cluster_similarity) of members

# start at one. if another cluster overlaps by 50%, merges these to a new list and remove from list. 
# Should just one match proportion or both?    
# If just one, then smaller groups will get rolled up into larger groups.
# If both, then fewer groups will get made.
propToMerge <- 0.5


## returns and index of the first pair of clusters to merge or 1 if no suitable pairs.
getVectorToMerge <- function(seedList)  {
	for(i in 2:length(seedList)) {
		#minLength <- min(length(seedList[[1]]),length(seedList[[i]]))	# use this for single group limit
		minLength <- max(length(seedList[[1]]),length(seedList[[i]]))	# use this for double group limit		
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








#BPterms <- ls(GOBPTerm)
#apply(head(BPterms),FUN=function(x) genesInTerm(GOdata.BP,x)[[1]])


## need to get ALL proteins from a GO term. 


## subset binaryGrid for GO terms with >= 10 proteins, should be around 469 terms