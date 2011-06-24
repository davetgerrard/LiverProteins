
#source()			
source("C:/Users/dave/LiverProteins/loadGoMaps.R")
nrow(proteinByLiverSample.ubi.pca)
nrow(allGo)
proteinByLiverSample.ubi.pca.go <- merge(proteinByLiverSample.ubi.pca,protGoMap,by.x="spAccession", by.y="DB_Object_ID") 
nrow(proteinByLiverSample.ubi.pca.go)

sampleGoTermCounts <- as.data.frame(table(as.character(proteinByLiverSample.ubi.pca.go$GO_ID)))
names(sampleGoTermCounts ) <- c("GO_ID","SampleGoWithProteinCount")
sampleGoTermCounts <- merge(sampleGoTermCounts,goTermCounts,by="GO_ID")
minSampleCounts <- 15
maxGlobalCounts <- 600

plot(sampleGoTermCounts$GoWithProteinCount,sampleGoTermCounts$SampleGoWithProteinCount)
abline(h=minSampleCounts)
abline(v=maxGlobalCounts)

sampleHighGo <- subset(sampleGoTermCounts, SampleGoWithProteinCount> minSampleCounts)
sampleHighGo <- subset(sampleHighGo, GoWithProteinCount < maxGlobalCounts)
sampleHighGo$GO_ID <- as.character(sampleHighGo$GO_ID) ## get rid of levels.
lengthHighGo <- nrow(sampleHighGo)
library("GO.db")
sampleHighGo$GO_NAME <- as.character(Term(as.character(sampleHighGo$GO_ID)))



pc.goTerms <- subset(merge(sampleHighGo,proteinByLiverSample.ubi.pca.go,by="GO_ID"),select=c(Comp.1,Comp.2,Comp.3,Comp.4,GO_ID,GO_NAME))
pc.goTerms$GO_ID <- as.character(pc.goTerms$GO_ID)
par(mar=c(5,20,5,2))
boxplot(pc.goTerms$Comp.1 ~ pc.goTerms$GO_NAME, horizontal=T,las=2, xlab="PC1 (Fetal : Adult)",xlim=c(0,lengthHighGo+1))
boxplot(pc.goTerms$Comp.2 ~ pc.goTerms$GO_NAME, horizontal=T,las=2, xlab="PC2 (HepG2 <-> Fetal)",xlim=c(0,lengthHighGo+1))
boxplot(pc.goTerms$Comp.3 ~ pc.goTerms$GO_NAME, horizontal=T,las=2, xlab="PC3 (Adult culturing)",xlim=c(0,lengthHighGo+1))
boxplot(pc.goTerms$Comp.4 ~ pc.goTerms$GO_NAME, horizontal=T,las=2, xlab="PC4 (Highspot <-> Other)",xlim=c(0,lengthHighGo+1))


##wilcox.tests on each group against the remainder. Make a table and adjust for multiple tests.
indexVersusRemainderWilcox <- function(dataVector,foreIndex)  {
	backIndex <- ifelse(foreIndex,FALSE,TRUE)
	pValue <- wilcox.test(dataVector[foreIndex],dataVector[backIndex ])[['p.value']]

}
indexGreaterThanRemainderWilcox <- function(dataVector,foreIndex)  {
	backIndex <- ifelse(foreIndex,FALSE,TRUE)
	pValue <- wilcox.test(dataVector[foreIndex],dataVector[backIndex ],alternative="g")[['p.value']]

}

indexVersusRemainderKS <- function(dataVector,foreIndex)  {
	backIndex <- ifelse(foreIndex,FALSE,TRUE)
	pValue <- ks.test(dataVector[foreIndex],dataVector[backIndex ])[['p.value']]
}
getPcaProteinIndexFromGoTerm <- function(thisGoTerm)  {
	thisGoTermProtList <- protGoMap$DB_Object_ID[protGoMap$GO_ID == thisGoTerm]
	thisGoTermProtIndexInPca <- as.numeric(na.omit(match(thisGoTermProtList, ubi.pca.5.scores$spAccession)))
	booleanIndex <- rep(FALSE,nrow(ubi.pca.5.scores))
	booleanIndex[thisGoTermProtIndexInPca] <- TRUE
	return(booleanIndex)
}


indexList <- lapply(levels(as.factor(sampleHighGo$GO_ID)),FUN=function(x) getPcaProteinIndexFromGoTerm(x))
countListVector <- as.numeric(lapply(indexList,sum))

pdf(file="liverProteinsPCvsGoTermColoured.pdf", width=10,height=10)
plot.new()
text(0.5,0.5,paste("Principal components scores for proteins grouped by GO terms.\n
		Only GO terms with >",
		minSampleCounts ,
		"proteins in the sample and <",
		maxGlobalCounts ,
		"proteins in the ontology are used.\n
		GO terms are redundant - the same protein may appear multiple times.\n
		Coloured boxes denote significantly higher absolute PC scores relative to remaining proteins in sample.\n
		(Wilcoxon test), no multiple correction. Yellow <0.05, Orange < 0.01, Red < 0.001",sep=" "))
par(mar=c(5,20,5,2))
pathNames <- sampleHighGo$GO_NAME[match(levels(as.factor(sampleHighGo$GO_ID)),sampleHighGo$GO_ID)]
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.1),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.1 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[1],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.1,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.1),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.2),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.2 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[2],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.2,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.2),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.3),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.3 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[3],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.3,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.3),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.4),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.4 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[4],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.4,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.4),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
dev.off()

pdf(file="liverProteinsPCvsGoTermColoured_KStest.pdf", width=10,height=10)
plot.new()
text(0.5,0.5,paste("Principal components scores for proteins grouped by GO terms.\n
		Only GO terms with >",
		minSampleCounts ,
		"proteins in the sample and <",
		maxGlobalCounts ,
		"proteins in the ontology are used.\n
		GO terms are redundant - the same protein may appear multiple times.\n
		Coloured boxes denote significantly higher absolute PC scores relative to remaining proteins in sample.\n
		(Wilcoxon test), no multiple correction. Yellow <0.05, Orange < 0.01, Red < 0.001",sep=" "))
par(mar=c(5,20,5,2))
pathNames <- sampleHighGo$GO_NAME[match(levels(as.factor(sampleHighGo$GO_ID)),sampleHighGo$GO_ID)]
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderKS(abs(ubi.pca.5.scores$Comp.1),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.1 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[1],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.1,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.1),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderKS(abs(ubi.pca.5.scores$Comp.2),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.2 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[2],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.2,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.2),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderKS(abs(ubi.pca.5.scores$Comp.3),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.3 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[3],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.3,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.3),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderKS(abs(ubi.pca.5.scores$Comp.4),x)))
colVector <- as.character(cut(pValList,c(1,0.05,0.01,0.001,0),labels=c("red","orange","yellow","white")))
boxplot(pc.goTerms$Comp.4 ~ pc.goTerms$GO_ID, horizontal=T,las=2, xlab=pc_defs[4],xlim=c(0,lengthHighGo+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.4,add=T,at=lengthHighGo+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.4),lty=2)
mtext(countListVector,side=4,at=1:lengthHighGo,las=2); mtext("ALL PROTEINS IN SAMPLE", side=2,at=lengthHighGo+1,las=2); mtext(side=3,"Protein count",adj=1)
dev.off()

