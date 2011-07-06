
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

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



#############more detritus 

#######THIS TRIAL SECTION WORKS
library(topGO)
library(org.Hs.eg.db) 
go2entrez <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "Entrez")
allGenes <- unique(unlist(go2entrez ))



someGenes <- sample(allGenes,400)
geneList <- 1:length(someGenes)
names(geneList) <- someGenes

GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = "BP",
              allGenes = geneList,
		  geneSelectionFun = topDiffGenes,
              nodeSize = 5,
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db",
              ID = "Entrez")

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.stat)

allRes <- GenTable(GOdata, KS = resultKS , elim=resultElim, orderBy = "elim", ranksOf =  "KS", topNodes=20)

#################
#prot_list <- list()
prot_list <- unique(protGoMap$DB_Object_ID)
	
listGoPerProt <- function(dataFrame,thisProt) {
	dataFrame$GO_ID[which(dataFrame$DB_Object_ID == thisProt)]
}

prot2go <- lapply(prot_list,FUN = function (x) listGoPerProt(protGoMap,x))
names(prot2go) <- prot_list
go2prot <- inverseList(prot2go)
prot2go <- inverseList(go2prot)


#topGoToProt <- annFUN.GO2genes("BP", feasibleGenes = NULL, go2prot)
#topGoProtToGo.BP <- annFUN.gene2GO("BP", feasibleGenes = NULL, prot2go )

allGenes <- unique(names(prot2go ))
someGenes <- sample(allGenes,400)
geneList <- 1:length(someGenes)
names(geneList) <- someGenes

#geneList <- 1:length(prot_list)
#names(geneList) <- prot_list


GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = "BP",
              allGenes = geneList,
		  geneSelectionFun = topDiffGenes,
              nodeSize = 5,
              annot = annFUN.GO2genes,
		GO2genes=go2prot 
		               )


test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.stat)

allRes <- GenTable(GOdata, KS = resultKS , elim=resultElim, orderBy = "elim", ranksOf =  "KS", topNodes=20)


######DEVELOPMENT



sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP",
allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.db,
affyLib = affyLib)

library(org.Hs.eg.db)  # not sure if this is already loaded.
raw.ent2up <- org.Hs.egUNIPROT
    # Get the entrez gene IDs that are mapped to a Uniprot ID
    mapped_genes <- mappedkeys(raw.ent2up)
    # Convert to a list
    ent2up <- as.list(raw.ent2up[mapped_genes])
      #ent2up [1:5]

#inverseList(ent2up)   ## could use this to get entrez IDs from Uniprot?
up2entrez <- inverseList(ent2up)
up2entrez[[proteinByLiverSample$spAccession]]
up2entrez[proteinByLiverSample$spAccession]

## some proteins have multiple gene identifiers. 


go2entrez <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "Entrez")
head(go2entrez)
en2go <- inverseList(go2entrez)

##
BPterms <- ls(GOBPTerm)   # loaded with topGO



GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO,
gene2GO = geneID2GO)
y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == "T")))



allGenes <- unique(unlist(go2entrez ))
#myInterestedGenes <- sample(allGenes, 500)
#geneList <- factor(as.integer(allGenes 
#names(geneList) <- allGenes

################# From Using GO
GOTERM$"GO:0019083"
GOBPANCESTOR$"GO:0019083"
GOBPPARENTS$"GO:0019083"
GOMFCHILDREN$"GO:0019083"
GOTERM$"GO:0019084"
GOTERM$"GO:0019085"
GOTERM$"GO:0006350"
GOTERM$"GO:0019080"
GOTERM$"GO:0000003"
GOBPPARENTS$"GO:0019084"
GOTERM$"GO:0009299"
GOTERM$"GO:0019083"





#######THIS TRIAL SECTION WORKS
library(topGO)
library(org.Hs.eg.db) 
go2entrez <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "Entrez")
allGenes <- unique(unlist(go2entrez ))

topDiffGenes <- function(allScore) {
  return(allScore > 1)
}

someGenes <- sample(allGenes,400)
geneList <- 1:length(someGenes)
names(geneList) <- someGenes

GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = "BP",
              allGenes = geneList,
		  geneSelectionFun = topDiffGenes,
              nodeSize = 5,
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db",
              ID = "Entrez")

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
resultElim <- getSigGroups(GOdata, test.stat)

allRes <- GenTable(GOdata, KS = resultKS , elim=resultElim, orderBy = "elim", ranksOf =  "KS", topNodes=20)

#goID <- allRes[1, "GO.ID"]
#print(showGroupDensity(GOdata, goID, ranks = TRUE)) ## doesn't work unless Bioconductor annotated chip.

#showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = "all")
#showSigOfNodes(GOdata, score(resultWeight), firstSigNodes = 5, useInfo = "def")
#printGraph(GOdata, resultKS, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#printGraph(GOdata, resultElim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE)
#####END OF TRIAL SECTION


##############from the vignette

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


