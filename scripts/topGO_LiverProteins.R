
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################


#use this source if running manually.
#source("C:/Users/dave/LiverProteins/loadPCA_noOutput.R")


if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

library(topGO)


#multiple runs. cycle through PCs, ontology terms.
pcs <- 1:4
goGraphs <- c("BP","MF","CC")
nodeSizeValue <- 10
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
elimCutOff <- 0.01   # 0.001 is too harsh
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


summaryDetectResultsList <- data.frame()

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
par(mar=c(1,1,1,1))
layout(matrix(c(1,2), byrow=T),heights=c(1,3)) 
textplot(headText,halign="left",mar=c(0,0,0,0))
textplot(allResFT,mar=c(0,0,0,0))

###  getMethods("GenTable")

### create summmary table manually and add to single table for all results across ontologies

manual.Fisher <- data.frame(Fisher=score(resultFisher),goTerm=names(score(resultFisher)))
manual.elimFisher <- data.frame(elimFisher=score(resultElimFisher),goTerm=names(score(resultElimFisher)))
manual.termStats <- termStat(GOdataBase,names(score(resultFisher)))
manual.termStats$goTerm <- row.names(manual.termStats)

thisGoResults <- NULL
thisGoResults <- merge(manual.termStats,manual.Fisher,by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimFisher ,by="goTerm")

thisGoResults$ontology <- thisGOgraph
thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(as.character(thisGoResults$goTerm), ontology(GOdataBase),numChar=200))
goGroup <- as.character(thisGoResults$goTerm)
#thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdataBase,x)[[1]])))

summaryDetectResultsList <- rbind(summaryDetectResultsList,thisGoResults)

}
dev.off()
######## END OF BASE GO ANALYSIS

write.table(summaryDetectResultsList[order(summaryDetectResultsList$elimFisher),],file="testDetectGOSummary.tab",quote=F,row.names=F,sep="\t")


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



summmaryPcResultList <- list()	# gonna have a separate data frame for each PC.
for (i in 1:5)  {
	summmaryPcResultList[[i]] <- data.frame()
}

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

### create summmary table manually and add to single table for all results across ontologies

manual.Wilcox <- data.frame(Wilcox=score(resultWilcox),goTerm=names(score(resultWilcox)))
manual.elimWilcox <- data.frame(elimWilcox=score(resultElimWilcox),goTerm=names(score(resultElimWilcox)))
manual.absWilcox <- data.frame(absWilcox =score(resultWilcoxAbs),goTerm=names(score(resultWilcoxAbs)))
manual.elimAbsWilcox <- data.frame(elimAbsWilcox=score(resultElimWilcoxAbs),goTerm=names(score(resultElimWilcoxAbs)))
manual.KS <- data.frame(KS=score(resultKS),goTerm=names(score(resultKS)))
manual.elimKS <- data.frame(elimKS=score(resultElimKS),goTerm=names(score(resultElimKS)))
manual.termStats <- termStat(GOdata,names(score(resultWilcox)))
manual.termStats <- subset(manual.termStats, select=-c(Significant,Expected))
manual.termStats$goTerm <- row.names(manual.termStats)

thisGoResults <- NULL
thisGoResults <- merge(manual.termStats, manual.Wilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.absWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimAbsWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.KS, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimKS, by="goTerm")

thisGoResults$ontology <- thisGOgraph
thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(as.character(thisGoResults$goTerm), ontology(GOdata),numChar=200))
goGroup <- as.character(thisGoResults$goTerm)
#thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdata,x)[[1]])))

summmaryPcResultList[[i]] <- rbind(summmaryPcResultList[[i]],thisGoResults)

}	#end of PC
}	#end of GO tree
dev.off()

for(i in 1:length(summmaryPcResultList)) {
	fileName <- paste("GoSummaryByPc",i,"tab",sep=".")
	write.table(summmaryPcResultList[[i]], file=fileName, sep="\t",row.names=F,quote=F)
}


stopifnot(FALSE)


#####################################STOP

####################################





#goID <- allRes[1, "GO.ID"]
#print(showGroupDensity(GOdata, goID, ranks = TRUE))

goID <- allResKS[1,"GO.ID"]
print(showGroupDensity(GOdata, goID, ranks = FALSE,rm.one=F))  #rm.one: natively expects scores between 0 and 1. Not true here.

# get the genes.
go.genes <- genesInTerm(GOdata,goID)[[1]]


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


## What does it mean?  Show genes for a GO term plotted back on some original data.
goID <- "GO:0005759"	# "mitochondrial matrix" top term for CC, PC4. (REst vs. Highspot)
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


#########
##show the PC1,4 biplot with only mitochondrial genes marked.
non.go.genes.index <- seq_len(length(proteinByLiverSample.ubiquitous$spAccession))[-go.genes.index]
mitoLabels <- proteinByLiverSample.ubiquitous$spAccession
mitoLabels[non.go.genes.index] <- "*" 
par(mfrow=c(1,2))
biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId,xlim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=mitoLabels,xlim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)



######  Pancreas development is showing up because of large number of ribosomal proteins ascribed to this term.

goID <- "GO:0031018"
go.genes <- genesInTerm(GOdata.BP,goID)[[1]]
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




