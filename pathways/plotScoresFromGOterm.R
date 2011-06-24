setwd("C:/Users/dave/LiverProteins/pathways/")



goID <- "GO:0005759"	# "mitochondrial matrix" top term for CC, PC4. (REst vs. Highspot)
goName <-  as.character(topGO:::.getTermsDefinition(goID, ontology(GOdata)))
go.genes <- genesInTerm(GOdata.CC,goID)[[1]]
go.genes.index <- match(go.genes , proteinByLiverSample.ubiquitous$spAccession)
non.go.genes.index <- seq_len(length(proteinByLiverSample.ubiquitous$spAccession))[-go.genes.index]

comb.go.index <- rep(goName,nrow(proteinByLiverSample.ubiquitous))
comb.go.index[non.go.genes.index ] <- "background"

#boxplot(proteinByLiverSample.ubiquitous[1,3:23] ~ comb.go.index, notch=T)

#perhaps boxplots would be better?
# want mean plus stdev for all these proteins.

freshFetalIndex
highSpotFetalIndex
freshAdultIndex
highSpotAdultIndex
#check the indices are correct for each tissue type.
names(proteinByLiverSample.ubiquitous)[freshFetalIndex]
names(proteinByLiverSample.ubiquitous)[highSpotFetalIndex]
names(proteinByLiverSample.ubiquitous)[freshAdultIndex]
names(proteinByLiverSample.ubiquitous)[highSpotAdultIndex]

boxplot(proteinByLiverSample.ubiquitous[go.genes.index,3:23], notch=T)




##### simple boxplot expression values. Using ALL samples
freshFetal.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,freshFetalIndex])
highSpotFetal.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,highSpotFetalIndex])
freshAdult.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,freshAdultIndex])
highSpotAdult.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,highSpotAdultIndex])

boxplot(freshFetal.scores ,notch=T, xlim=c(0,5), ylim=c(0,3))
boxplot(highSpotFetal.scores ,notch=T, add=T, at=2)
boxplot(freshAdult.scores ,notch=T, add=T, at=3)
boxplot(highSpotAdult.scores ,notch=T, add=T, at=4)

## ALSO WITH BACKGROUND
freshFetal.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,freshFetalIndex])
highSpotFetal.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,highSpotFetalIndex])
freshAdult.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,freshAdultIndex])
highSpotAdult.scores <- unlist(proteinByLiverSample.ubiquitous[go.genes.index,highSpotAdultIndex])

freshFetal.bg <- unlist(proteinByLiverSample.ubiquitous[non.go.genes.index,freshFetalIndex])
highSpotFetal.bg <- unlist(proteinByLiverSample.ubiquitous[non.go.genes.index,highSpotFetalIndex])
freshAdult.bg <- unlist(proteinByLiverSample.ubiquitous[non.go.genes.index,freshAdultIndex])
highSpotAdult.bg <- unlist(proteinByLiverSample.ubiquitous[non.go.genes.index,highSpotAdultIndex])

boxplot(freshFetal.bg ,notch=T, xlim=c(0,8), ylim=c(0,4))
boxplot(freshFetal.scores ,notch=T, add=T, at=2)
boxplot(highSpotFetal.bg  ,notch=T, add=T, at=3)
boxplot(highSpotFetal.scores ,notch=T, add=T, at=4)
boxplot(freshAdult.bg ,notch=T, add=T, at=5)
boxplot(freshAdult.scores ,notch=T, add=T, at=6)
boxplot(highSpotAdult.bg ,notch=T, add=T, at=7)
boxplot(highSpotAdult.scores ,notch=T, add=T, at=8)


########### show focal and non-focal for each tissue (to allow for sample differences).






###### Could do plot of fetal vs. adult with mean values. 


##### expression difference values.
hsPairs <- c("FLM165","FHM165", "FLM213","FHM213","FLM257","FHM257","FA1","AH1","FA3","AH3","FA4","AH4")


diffMatrix <- subset(proteinByLiverSample.ubiquitous, select=hsPairs)

diffMatrix$diff.F165 <- diffMatrix$FLM165 - diffMatrix$FHM165
diffMatrix$diff.F213 <- diffMatrix$FLM213 - diffMatrix$FHM213
diffMatrix$diff.F257 <- diffMatrix$FLM257 - diffMatrix$FHM257
diffMatrix$diff.A1 <- diffMatrix$FA1 - diffMatrix$AH1
diffMatrix$diff.A3 <- diffMatrix$FA3 - diffMatrix$AH3
diffMatrix$diff.A4 <- diffMatrix$FA4 - diffMatrix$AH4

plot.diffMatrix <- diffMatrix[,grep("diff.",names(diffMatrix))]
boxplot(plot.diffMatrix , notch=T)
# weird effect of difference of ratios. Not good.


##### simple boxplot expression values. Using only paired samples
proteinByLiverSample.ubi.hsPairs <- subset(proteinByLiverSample.ubiquitous, select=hsPairs)

freshFetal.scores.hs <- unlist(proteinByLiverSample.ubi.hsPairs[go.genes.index,grep("FLM",names(proteinByLiverSample.ubi.hsPairs ))])
highSpotFetal.scores.hs <- unlist(proteinByLiverSample.ubi.hsPairs[go.genes.index,grep("FHM",names(proteinByLiverSample.ubi.hsPairs ))])
freshAdult.scores.hs <- unlist(proteinByLiverSample.ubi.hsPairs[go.genes.index,grep("FA",names(proteinByLiverSample.ubi.hsPairs ))])
highSpotAdult.scores.hs <- unlist(proteinByLiverSample.ubi.hsPairs[go.genes.index,grep("AH",names(proteinByLiverSample.ubi.hsPairs ))])

boxplot(freshFetal.scores.hs,notch=T, xlim=c(0,5), ylim=c(0,3) ) ; mtext("fetal\nfresh",side=1,at=1, line=2)
boxplot(highSpotFetal.scores.hs,notch=T, add=T, at=2); mtext("fetal\nhigh-spot",side=1,at=2, line=2)
boxplot(freshAdult.scores.hs,notch=T, add=T, at=3); mtext("adult\nfresh",side=1,at=3, line=2)
boxplot(highSpotAdult.scores.hs,notch=T, add=T, at=4); mtext("adult\nhigh-spot",side=1,at=4, line=2)




##### show all columns ordered into paired/triple samples.
orderedColIndex <- c("FLM134","FLM154","FLM158","FLM165","FHM165","FLM169","FLM171",
				"FLM213","FHM213","FLM257","FHM257","FA1","AM1","AH1","FA2",
				"FA3","AM3","AH3","FA4","AM4","AH4")

plotMatrix <- subset(proteinByLiverSample.ubiquitous, select=orderedColIndex)
#boxplot(plotMatrix[go.genes.index,], notch=T)
plot.HS <- c(grep("FHM", names(plotMatrix)),grep("AH", names(plotMatrix)))
plot.M <- grep("AM", names(plotMatrix))
plot.colours <- rep("white", length(names(plotMatrix)))
plot.colours[plot.HS] <- "black"
plot.colours[plot.M] <- "grey"
boxplot(plotMatrix[go.genes.index,], notch=T,col=plot.colours)

# use subset with pairs only
plotMatrix.hsPairs <- subset(plotMatrix, select=c("FLM165","FHM165","FLM213","FHM213","FLM257","FHM257",
				"FA1","AH1","FA3","AH3","FA4","AH4"))
plot.colours.hsPairs <- rep(c("white","black"),length(names(plotMatrix.hsPairs))/2)
par(mar=c(7,4,3,2), cex=1.5)
boxplot(plotMatrix.hsPairs[go.genes.index,], notch=T,col=plot.colours.hsPairs ,las=3,main=paste("GO term: ",goName),ylab="relative protein intensity")
legend("topleft",c("Fresh","High-spot"),fill=c("white","black"))
mtext("Fetal", side=1, at=3.5, line=5,cex=1.5)
mtext("Adult", side=1, at=9.5, line=5,cex=1.5)
abline(h=1,lty=2)
####### 

