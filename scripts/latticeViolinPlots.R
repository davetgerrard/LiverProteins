
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

## custom violin/boxplot for showing distribution of scores.

library(lattice)



##dev.code
#archetype.
bwplot(voice.part ~ height, singer,
       panel = function(..., box.ratio) {
           panel.violin(..., col = "transparent",
                        varwidth = FALSE, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } )





#### to use this, create a datatable of all scores with groups labeled
baseScore <- data.frame(score=ubi.pca.5.scores[,"Comp.1"],group="all")

goID <- "GO:0006805"	# "xenobiotic metabolic..."
go.genes.xeno <- genesInTerm(GOdata.BP,goID)[[1]]
goID <- "GO:0005759"	# mitochondrial matrix
go.genes <- genesInTerm(GOdata.CC,goID)[[1]]

addGroupScoresToGroupTable <- function(groupTable,groupIds,groupName,scoreTable,idColumn,scoreColum) {
	index <- na.omit(match(groupIds , scoreTable[,idColumn]))
	scores <- scoreTable[index ,scoreColum]
	groupTable <- rbind(groupTable,data.frame(score=scores,group=groupName))
	
}

baseScore <- data.frame(score=ubi.pca.5.scores[,"Comp.1"],group="all")
groupTable <- baseScore
groupTable <- addGroupScoresToGroupTable(groupTable=groupTable,groupIds=go.genes,
				groupName=goID,scoreTable=ubi.pca.5.scores,
				idColumn="spAccession",scoreColum="Comp.1")

bwplot(group ~ score, groupTable,
       panel = function(..., box.ratio) {
           panel.violin(..., col = "transparent",
                        varwidth = FALSE, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } )