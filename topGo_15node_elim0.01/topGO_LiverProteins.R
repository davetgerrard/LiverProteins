library(topGO)


source("C:/Users/dave/LiverProteins/loadPCA_noOutput.R")
setwd("C:/Users/dave/LiverProteins/topGo_15node_elim0.01/")

#multiple runs. cycle through PCs, ontology terms.
pcs <- 1:4
goGraphs <- c("BP","MF","CC")
nodeSizeValue <- 15
topNodesValue <- 20

#ubi.pca.5.scores

#prot_list <- list()
#prot_list <- unique(protGoMap$DB_Object_ID)
	
prot2go <-  readMappings("C:/Users/dave/LiverProteins/data/go2prot.map")
go2prot <- inverseList(prot2go)


####
# insert standard fisher test of over-represented terms.
require(gplots)
library(qvalue)
elimCutOff <- 0.01
topTerms <- 50
resultSummaryAsText <-function (x)    # Function to collect data from GOdata results.
{
    text <-     paste("\nDescription:", description(x), "\n",
    	"'", algorithm(x), "' algorithm with the '", testName(x), 
        "' test\n", 
    length(score(x)), "GO terms scored:", sum(score(x) <= 
        elimCutOff), "terms with p < ",elimCutOff ,"\n")
    xg <- geneData(x)
    if ("Annotated" %in% names(xg)) {
    	text <- paste(text,"    Annotated proteins:", xg["Annotated"], "\n")
    	}
    if ("Significant" %in% names(xg)) {
            text <- paste(text,"    Significant proteins:", xg["Significant"], "\n")
        }
    text <- paste(text,"    Min. no. of proteins annotated to a GO:", xg["NodeSize"],"\n")
    if ("SigTerms" %in% names(xg)) {
        text <- paste(text,"    Nontrivial nodes:", xg["SigTerms"], "\n")
        }
    return(text) ;
    #.printGeneData(geneData(x))
}

pdf(file="liverProteinsUbiqBasicGo.pdf", width=10,height=9)
for(thisGOgraph in goGraphs )  {
allProtList <- names(prot2go)
detectedProts <- ubi.pca.5.scores$spAccession
geneList <- factor(as.integer(allProtList %in% detectedProts))
names(geneList) <- allProtList 
#head(geneList)
GOdataBase <- new("topGOdata",
  		description =  "Liver proteins data set",
              ontology = thisGOgraph,
              allGenes = geneList,
  			nodeSize = nodeSizeValue,
              annot = annFUN.GO2genes,
		GO2genes=go2prot 
               )
test.stat <- new("classicCount", testStatistic = GOFisherTest, name="Fisher test")
resultFisher <- getSigGroups(GOdataBase,test.stat)
test.stat <- new("elimCount", testStatistic = GOFisherTest, name="elimFisher", cutOff = elimCutOff)
resultElimFisher <- getSigGroups(GOdataBase,test.stat)
resultElimFisherAdj <- resultElimFisher
score(resultElimFisherAdj) <- qvalue(score(resultElimFisher))$qvalue
allResFT <- GenTable(GOdataBase, classic=resultFisher,elim=resultElimFisher,qvalElim=resultElimFisherAdj,
				orderBy="qvalElim",topNodes=topTerms )
#allResFT <- subset(allResFT, select=-"Rank in elim")
#allResFT
headText <- resultSummaryAsText(resultElimFisher)
headText <- paste(headText, "\nTop ",topTerms,"GO terms shown")
par(mar=c(0,0,0,0))
layout(matrix(c(1,2), byrow=T),heights=c(1,3)) 
textplot(headText,halign="left",mar=c(0,0,0,0))
textplot(allResFT,mar=c(0,0,0,0))
}
dev.off()
######## END OF BASE GO ANALYSIS


GOWilcoxTestGreater <- function (object,alternativeType="greater") 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}



GOWilcoxTest2Sided <- function (object,alternativeType="two.sided") 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}

#resultElimFisherAdj <- resultElimFisher
#score(resultElimFisherAdj) <- qvalue(score(resultElimFisher))$qvalue
#allResFT <- GenTable(GOdataBase, classic=resultFisher,elim=resultElimFisher,qvalElim=resultElimFisherAdj,#
#				orderBy="qvalElim",topNodes=topTerms )
#allResFT

#sample
i<-1
thisGOgraph <- "BP"
thisGOgraph <- "MF"
thisGOgraph <- "CC"


## faster to specify the GOdata object once for each ontology and update scores from each PC.
pdf(file="liverProteinsUbiqPCA_GO_byWilcox.pdf", width=11,height=9)
for(thisGOgraph in goGraphs )  {

#topGoToProt <- annFUN.GO2genes("BP", feasibleGenes = NULL, go2prot)
#topGoProtToGo.BP <- annFUN.gene2GO("BP", feasibleGenes = NULL, prot2go )

geneList <- ubi.pca.5.scores$Comp.1
names(geneList) <- ubi.pca.5.scores$spAccession

topDiffGenes <- function(allScore) {
  return(allScore )
}

GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = thisGOgraph ,
              allGenes = geneList,
		  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
		GO2genes=go2prot 
		               )
for (i in 1:5)  {
compHead <- paste("Comp.",i,sep="")
geneList <- ubi.pca.5.scores[,compHead]
names(geneList) <- ubi.pca.5.scores$spAccession
GOdata <- updateGenes(GOdata,geneList,topDiffGenes)

test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
resultWilcox <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcox <- getSigGroups(GOdata, test.stat)

#resultElimWilcoxAdj <- resultElimWilcox 
#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
resultElimKS <- getSigGroups(GOdata, test.stat)

#resultElimKSAdj <- resultElimKS
#score(resultElimKSAdj ) <- qvalue(score(resultElimKS))$qvalue

# change scores to absolute values for Greater than test
geneList <- abs(ubi.pca.5.scores[,compHead])
names(geneList) <- ubi.pca.5.scores$spAccession
GOdata <- updateGenes(GOdata,geneList,topDiffGenes)
test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
resultWilcoxAbs <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)

#resultElimWilcoxAbsAdj <- resultElimWilcoxAbs 
#score(resultElimWilcoxAbsAdj ) <- qvalue(score(resultElimWilcoxAbs ))$qvalue


allResKS <- GenTable(GOdata, Wilcox = resultWilcox,  elimWilcox = resultElimWilcox ,
					absWilcox = resultWilcoxAbs, elimAbsWilcox=resultElimWilcoxAbs ,
					KS = resultKS , elimKS=resultElimKS, 
					 orderBy = "elimAbsWilcox", 
					topNodes=topTerms)
allResKS <- subset(allResKS, select=-c(Annotated,Significant,Expected))
headText <- paste("Principal Component: ",i,"\n",resultSummaryAsText(resultElimWilcox))
headText <- paste(headText, "\nTop ",topTerms,"GO terms shown")
par(mar=c(0,0,0,0))
layout(matrix(c(1,2), byrow=T),heights=c(1,3)) 
textplot(headText,halign="left",mar=c(0,0,0,0))
textplot(allResKS ,mar=c(0,0,0,0))
}	#end of PC
}	#end of GO tree
dev.off()

###############################################

#goID <- allRes[1, "GO.ID"]
#print(showGroupDensity(GOdata, goID, ranks = TRUE))

goID <- allResKS[1,"GO.ID"]
print(showGroupDensity(GOdata, goID, ranks = FALSE,rm.one=F))  #rm.one: natively expects scores between 0 and 1. Not true here.

# get the genes.
go.genes <- genesInTerm(GOdata,goID)[[1]]

## What does it mean?  Show genes for a GO term plotted back on some original data.
goID <- "GO:0005759"	# "mitochondrial matrix" top term for CC, PC4. (REst vs. Highspot)
thisGOgraph <- "CC"
GOdata <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = thisGOgraph ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )
### three separate GOdata structures for quickly pulling down genes with a GO term
GOdata.BP <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "BP" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.MF <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "MF" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.CC <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = "CC" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

goID <- "GO:0019083"	# viral transcription
goID <- "GO:0031018"	# endocrine pancreas development
goID <- "GO:0006414"	# translation elongation
goID <- "GO:0006415"	# translation termination

goID <- "GO:0055114"	# "oxidation reduction"
go.genes.oxRed <- genesInTerm(GOdata.BP,goID)[[1]]
goID <- "GO:0006094"	# "gluconeogenesis"
go.genes.gluco <- genesInTerm(GOdata.BP,goID)[[1]]
goID <- "GO:0006805"	# "xenobiotic metabolic..."
go.genes.xeno <- genesInTerm(GOdata.BP,goID)[[1]]
length(go.genes.oxRed);length(go.genes.gluco);length(go.genes.xeno)
length(intersect(go.genes.oxRed,go.genes.gluco))
length(intersect(go.genes.oxRed,go.genes.xeno))

goID <- "GO:0001889"  # won't load because less than 15 proteins in dataset annotated to "Liver development" (we detect 11/64)

goID <- "GO:0042470"  # CC: melanosome  ?

goID <- "GO:0020037"	# MF: Heme binding  (hemoglobins plus cytochromes)


goName <-  as.character(topGO:::.getTermsDefinition(goID, ontology(GOdata)))
go.genes <- genesInTerm(GOdata.BP,goID)[[1]]
go.genes <- genesInTerm(GOdata.MF,goID)[[1]]
go.genes <- genesInTerm(GOdata.CC,goID)[[1]]



#goName <-  as.character(topGO:::.getTermsDefinition(goID, ontology(GOdata)))
#go.genes <- genesInTerm(GOdata,goID)[[1]]
go.genes.index <- match(go.genes , proteinByLiverSample.ubiquitous$spAccession)
proteinByLiverSample.ubiquitous[go.genes.index,c("Accession", "Name")]
non.go.genes.index <- seq_len(length(proteinByLiverSample.ubiquitous$spAccession))[-go.genes.index]
par(mfrow=c(2,2),cex=1.0)
plot( proteinByLiverSample.ubiquitous$FA4[non.go.genes.index] ~ proteinByLiverSample.ubiquitous$FLM165[non.go.genes.index],
		xlab="FLM165, fresh fetal", ylab="FA4, fresh adult", main=goName )
abline(b=1,a=0)
points( proteinByLiverSample.ubiquitous$FA4[go.genes.index]~ proteinByLiverSample.ubiquitous$FLM165[go.genes.index], pch=17,col="blue")
plot( proteinByLiverSample.ubiquitous$AH4[non.go.genes.index] ~ proteinByLiverSample.ubiquitous$FHM165[non.go.genes.index],
		xlab="FHM165, high-spot fetal", ylab="AH4, high-spot adult", main=goName )
abline(b=1,a=0)
points( proteinByLiverSample.ubiquitous$AH4[go.genes.index]~ proteinByLiverSample.ubiquitous$FHM165[go.genes.index], pch=17,col="blue")
plot( proteinByLiverSample.ubiquitous$AH4[non.go.genes.index] ~ proteinByLiverSample.ubiquitous$FA4[non.go.genes.index],
		xlab="FA4, fresh adult", ylab="AH4, high-spot adult", main=goName )
abline(b=1,a=0)
points( proteinByLiverSample.ubiquitous$AH4[go.genes.index]~ proteinByLiverSample.ubiquitous$FA4[go.genes.index], pch=17,col="blue")
plot( proteinByLiverSample.ubiquitous$FHM165[non.go.genes.index] ~ proteinByLiverSample.ubiquitous$FLM165[non.go.genes.index],
		xlab="FLM165, fresh fetal", ylab="FHM165, high-spot fetal", main=goName)
abline(b=1,a=0)
points( proteinByLiverSample.ubiquitous$FHM165[go.genes.index] ~ proteinByLiverSample.ubiquitous$FLM165[go.genes.index], pch=17,col="blue")


######  Pancreas development is showing up because of large number of ribosomal proteins ascribed to this term.

goID <- "GO:0031018"
go.genes <- genesInTerm(GOdata,goID)[[1]]
thisGOgraph <- "BP"
GOdata <- new("topGOdata",
  description =  "Liver proteins data set",
              ontology = thisGOgraph ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenes,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )
go.genes <- genesInTerm(GOdata,goID)[[1]]
go.genes
go.genes.index <- match(go.genes , proteinByLiverSample.ubiquitous$spAccession)
proteinByLiverSample.ubiquitous[go.genes.index,]
proteinByLiverSample.ubiquitous[go.genes.index,c("Accession", "name")]
proteinByLiverSample.ubiquitous[go.genes.index,c("Accession", "Name")]










#geneList <- abs(ubi.pca.5.scores$Comp.4)
##names(geneList) <- ubi.pca.5.scores$spAccession
#GOdata <- updateGenes(GOdata,geneList,topDiffGenes)


#topGO:::.printGeneData
## how to put in wilcox.test?
#test.stat <- new("classicScore", testStatistic = function() wilcox.test(), name = "Wilcox tests")
#resultWT <- getSigGroups(GOdata, test.stat)







#####NEED to work out if these results are in anyway meaningful.....
#### the significant are all ZERO.
#### Also, whether using absolute value would be better.





stopifnot(FALSE)

#####################




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


