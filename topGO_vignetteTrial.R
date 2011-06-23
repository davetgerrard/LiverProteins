library(topGO)
library(ALL)
data(ALL)
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(geneList))
sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP",
allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.db,
affyLib = affyLib)
sampleGOdata
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS,
elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher",
topNodes = 10)
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated/max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
pch = 19, cex = gSize, col = gCol)
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
pch = 19, cex = gSize, col = gCol)


BPterms <- ls(GOBPTerm)
head(BPterms)
library(genefilter)
selProbes <- genefilter(ALL, filterfun(pOverA(0.2, log2(100)), function(x) (IQR(x) >
0.25)))
eset <- ALL[selProbes, ]
selProbes
eset
geneID2GO <- readMappings(file="C:/Users/dave/Documents/R/win-library/2.11/topGO/examples/ensembl2go.map") ## replaced the file name here
geneID2GO
str(head(geneID2GO))
GO2geneID <- inverseList(geneID2GO)
str(head(GO2geneID))
geneNames <- names(geneID2GO)
head(geneNames)
myInterestingGenes <- sample(geneNames, length(geneNames)/10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO,
gene2GO = geneID2GO)
GOdata
y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == "T")))
table(y)
y
geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")
geneList
topDiffGenes <- function(allScore) {
	return(allScore < 0.01)
}
x <- topDiffGenes(geneList)
sum(x)
GOdata <- new("topGOdata", description = "GO analysis of ALL data; B-cell vs T-cell",
	ontology = "BP", allGenes = geneList, geneSel = topDiffGenes, annot = annFUN.db,
	nodeSize = 5, affyLib = affyLib)


