#############################################
#
# Dave Gerrard, University of Manchester
#
#############################################


setwd("C:/Users/dave/LiverProteins/pathways/")

output <- FALSE

dir()
baseProteinByLiverSample <- read.delim("C:/Users/dave/LiverProteins/Cliffs dataset.txt", header=T)

summary(baseProteinByLiverSample)
nrow(baseProteinByLiverSample)

# There are 2304 rows and many NAs (missing values). 
# Only a small number of proteins are present in all samples.
# Some samples have a maximum value of 9999. 

# make a copy of the dataset.
proteinByLiverSample <- baseProteinByLiverSample



############ MANIPULATIONS to clean and index data.
#!# Probably should split the Accession into three
proteinByLiverSample$spAccession <- matrix(unlist(strsplit(as.character(proteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]
proteinByLiverSample$spId <- matrix(unlist(strsplit(as.character(proteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,3]


# The second column has an odd name
names(proteinByLiverSample)[grep('Name.x.x',names(proteinByLiverSample))] <- "Name"
# One column is FL257 instead of FLM257
names(proteinByLiverSample)[grep('FL257',names(proteinByLiverSample))] <- "FLM257"
# The final column is not data. 
proteinByLiverSample <- subset(proteinByLiverSample, select= -X)

##
##THIS VERSION: REMOVE SAMPLES 290/292
##
 proteinByLiverSample <- subset(proteinByLiverSample, select=c(-FLM290,-FHM290,-FLM292,-FHM292))


##set up name indexes to access functional groups.
freshFetalIndex <- grep('FLM',names(proteinByLiverSample))
highSpotFetalIndex <- grep('FHM',names(proteinByLiverSample))
freshAdultIndex <- grep('FA',names(proteinByLiverSample))
highSpotAdultIndex <- grep('AH',names(proteinByLiverSample))
matrigelAdultIndex <- grep('AM',names(proteinByLiverSample))
fetalIndex <- c(freshFetalIndex,highSpotFetalIndex)	#excluded HepG2
adultIndex <- c(freshAdultIndex,highSpotAdultIndex ,matrigelAdultIndex ) 
hepG2Index <- grep('HepG2',names(proteinByLiverSample))
allTissueIndex <- c(fetalIndex ,adultIndex ,hepG2Index )
unCulturedIndex <- c(freshFetalIndex,freshAdultIndex,hepG2Index)

# One further manipulation: check which proteins have maxValue of 9999.
proteinByLiverSample$maxValue <- apply(proteinByLiverSample[,allTissueIndex] , 1,max,na.rm=T)
proteinByLiverSample[proteinByLiverSample$maxValue > 20,]
# there are only three none are ubiquitous. Remove from dataset.
proteinByLiverSample <- subset(proteinByLiverSample, maxValue < 21)

#!# Should probably load the other experimental factors like age and MS run.


############ QQ-plots to compare intensity value ranges across samples.
if(output) {
pdf(file="liverProteinsSampleBasicQQplots.pdf", width=10,height=10)
op <- par(mfrow=c(3,3))
numbSamples <- length(allTissueIndex)
for(i in 1:(numbSamples-1)) {
	for(j in (i+1):numbSamples) {
		qqplot(proteinByLiverSample[,allTissueIndex[i]],proteinByLiverSample[,allTissueIndex[j]],
			xlab=names(proteinByLiverSample)[allTissueIndex[i]],
			ylab=names(proteinByLiverSample)[allTissueIndex[j]])
		abline(a=0,b=1)
	}
}
par(op)
dev.off()
}

## count presence absence across all, across fetal/adult separately and across cultured?
proteinByLiverSample$detectCountAll <- length(allTissueIndex) -  rowSums(is.na(proteinByLiverSample[,allTissueIndex]))
proteinByLiverSample$detectCountFetal <- length(fetalIndex) -  rowSums(is.na(proteinByLiverSample[,fetalIndex]))
proteinByLiverSample$detectCountAdult <- length(adultIndex) -  rowSums(is.na(proteinByLiverSample[,adultIndex]))
proteinByLiverSample$detectCountHepG2 <- length(hepG2Index) -  rowSums(is.na(proteinByLiverSample[,hepG2Index]))
proteinByLiverSample$detectCountAM <- length(matrigelAdultIndex) -  rowSums(is.na(proteinByLiverSample[,matrigelAdultIndex]))
proteinByLiverSample$detectCountAH <- length(highSpotAdultIndex ) -  rowSums(is.na(proteinByLiverSample[,highSpotAdultIndex ]))
proteinByLiverSample$detectCountFHM <- length(highSpotFetalIndex ) -  rowSums(is.na(proteinByLiverSample[,highSpotFetalIndex ]))
proteinByLiverSample$detectCountFLM <- length(freshFetalIndex ) -  rowSums(is.na(proteinByLiverSample[,freshFetalIndex ]))
proteinByLiverSample$detectCountFA <- length(freshAdultIndex ) -  rowSums(is.na(proteinByLiverSample[,freshAdultIndex ]))


# what was I doing here? 
# I want to know if genes missing from some tissues show similarity within a tissue type. 
# Expect high count in fetal with low count in adults and vice versa.
# Don't reallly see it though. May as well show with heatmaps.
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountAdult,xlab="Total tissue count",ylab="Adult count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountFetal,xlab="Total tissue count",ylab="Fetal count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountHepG2,xlab="Total tissue count",ylab="HepG2 count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountAM,xlab="Total tissue count",ylab="Adult Matrigel Count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountAH,xlab="Total tissue count",ylab="Adult High-spot count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountFHM,xlab="Total tissue count",ylab="Fetal High-spot count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountFA,xlab="Total tissue count",ylab="Fresh Adult count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountFLM,xlab="Total tissue count",ylab="Fresh Fetal count") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountAdult, proteinByLiverSample$detectCountFetal,xlab="Adult tissue count",ylab="") ; abline(a=0,b=1)
#plot(proteinByLiverSample$detectCountHepG2, proteinByLiverSample$detectCountFetal,xlab="HepG2 tissue count",ylab="") ; abline(a=0,b=1)


## Make some subsets of data based on presence/absence of proteins across samples

# Subset of data with proteins detected in all samples
proteinByLiverSample.ubiquitous <- subset(proteinByLiverSample, detectCountAll == length(allTissueIndex))
numbUbiquitous  <- nrow(proteinByLiverSample.ubiquitous)

# Subset of data with proteins missing from one or more samples. Convert NA values to zero. 
proteinByLiverSample.selective <- subset(proteinByLiverSample, detectCountAll < length(allTissueIndex))  # <28
proteinByLiverSample.selective <- subset(proteinByLiverSample.selective, detectCountAll > 2)
proteinByLiverSample.selective[is.na(proteinByLiverSample.selective)] <- 0
numbSelective <- nrow(proteinByLiverSample.selective)

## Plot as HEATMAPS.
if(output) {
pdf(file="liverProteinsSampleBasicHeatmaps.pdf", width=15,height=10)
heatmap(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex]),labRow=proteinByLiverSample.ubiquitous$Name)
mtext(paste(numbUbiquitous,"Ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
heatmap(as.matrix(proteinByLiverSample.ubiquitous[unCulturedIndex]),labRow=proteinByLiverSample.ubiquitous$Name)
mtext(paste(numbUbiquitous,"Pre-cultured tissues"),side=3,line=2,adj=0,cex=1.5)
heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name)
mtext(paste(numbSelective,"Non-ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name,distfun=function(c) dist(c,method="binary"))
mtext(paste(numbSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=2,adj=0,cex=1.5)
# The tissue clustering is not replicated using the set of proteins shared by 10 to 27 tissues.
dev.off()  # end of "liverProteinsSampleBasicHeatmaps.pdf"
}

############ Pairwise plots of ubiquitous proteins to compare samples.
## Shows consistent positive correlations within fetal (+ HepG2) and within Adult but negative correlations between adult and fetal.
if(output) {
pdf(file="liverProteinsSampleUbiquitousPairplots.pdf", width=10,height=10)
op <- par(mfrow=c(3,3))
numbSamples <- length(allTissueIndex)
for(i in 1:(numbSamples-1)) {
	for(j in (i+1):numbSamples) {
		plot(proteinByLiverSample.ubiquitous[,allTissueIndex[i]],proteinByLiverSample.ubiquitous[,allTissueIndex[j]],
			xlab=names(proteinByLiverSample.ubiquitous)[allTissueIndex[i]],
			ylab=names(proteinByLiverSample.ubiquitous)[allTissueIndex[j]])
		abline(a=0,b=1)
	}
}
par(op)
dev.off()
}

########## PRINCIPAL COMPONENTS ANALYSIS

ubiquitous.pca <- princomp(proteinByLiverSample.ubiquitous[,allTissueIndex])
ubiquitous.pca
summary(ubiquitous.pca)

loadings(ubiquitous.pca)  
# PCA on with 431 ubiquitous proteins,
# PC1 (43%) separates adult from (fetal + HepG2) 
# pc2 (12%) separates fetal and HepG2
# pc3 (9%) separates fresh adult tissue from Matrigel
# pc4 (8%) partially separates high-spot + hepG2 from the rest
# pc5 (5%) partially separates high-spot from the remainder. (see 4 vs. 5)

if(output) {
pdf(file="liverProteinsSampleUbiquitousPCA.pdf", width=10,height=10)
plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(3,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(4,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
dev.off()	
}
## Make lists of genes correlated with each major principal component analysis.

ubi.pca.5.scores <- data.frame(ubiquitous.pca$scores[ , 1:5])
ubi.pca.5.scores$spId <- proteinByLiverSample.ubiquitous$spId
ubi.pca.5.scores$Name <- proteinByLiverSample.ubiquitous$Name
ubi.pca.5.scores$spAccession <- proteinByLiverSample.ubiquitous$spAccession
if(output) {
write.table(ubi.pca.5.scores, file="UbiquitousLiverProteinsFirst5PcaScores.tab",sep="\t",row.names=F,quote=F)
}
proteinByLiverSample.ubi.pca <- merge(proteinByLiverSample.ubiquitous,ubi.pca.5.scores,by="spAccession")

#!# NETWORK ANALYSIS
pc_defs <- c("PC1 (47%) Fetal+HepG2 <-> Adult", 
		"PC2 (13%) HepG2 <-> Fetal",
		"PC3 (11%) ADULT: Fresh <- High-spot -> Matrigel ",
		"PC4 (6%) Non-high-spot <-> High-spot")



source("C:/Users/dave/LiverProteins/loadKeggMaps.R")
proteinByLiverSample.ubi.pca.kegg <- merge(proteinByLiverSample.ubi.pca, protKeggMap, by.x="spAccession", by.y="uniprotId")

hist(proteinByLiverSample.ubi.pca.kegg$keggPathsWithProtCount,breaks=20)
# which pathways are represented and by how many proteins.
sampleKeggPathCounts <- as.data.frame(table(as.character(merge(proteinByLiverSample.ubi.pca.kegg,keggProtPathMap,by="keggProtId", all.x=T)$keggPathId)))
names(sampleKeggPathCounts) <- c("keggPathId","sampleCount")
sampleKeggPathCounts <- merge(sampleKeggPathCounts,pathCounts, by="keggPathId")
# some paths are highly represented in dataset already.
plot(sampleKeggPathCounts$keggProtsInPathCount, sampleKeggPathCounts$sampleCount)

sampleHighPaths <- subset(sampleKeggPathCounts, sampleCount > 10)
sampleHighPaths <- subset(sampleHighPaths , keggProtsInPathCount < 200)

#make table of PCs and all paths. (proteins duplicated).
pc.paths <- subset(merge(proteinByLiverSample.ubi.pca.kegg,keggProtPathMap,by="keggProtId"),select=c(Comp.1,Comp.2,Comp.3,Comp.4,keggPathId))

#filter on main paths of dataset

pc.paths.high <- merge(sampleHighPaths, pc.paths, by="keggPathId")
pc.paths.high$keggPathId <- as.factor(as.character(pc.paths.high$keggPathId))

if(output) {
pdf(file="liverProteinsPCvsKEGG.pdf", width=10,height=10)
par(mar=c(5,20,5,2))
boxplot(pc.paths.high$Comp.1 ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC1 (Fetal <-> Adult)") ; abline(v=0)
boxplot(pc.paths.high$Comp.2 ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC2 (HepG2 <-> Fetal)") ; abline(v=0)
boxplot(pc.paths.high$Comp.3 ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC3 (Adult culturing)") ; abline(v=0)
boxplot(pc.paths.high$Comp.4 ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC4 (Highspot <-> Other)") ; abline(v=0)
dev.off()
}



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


# get protein index for this kegg path from orignial prcomp list

getPcaProteinIndexFromKeggPath <- function(thisKeggPathId)  {
	thisPathProtList <- keggPathToUniprot$uniprotId[keggPathToUniprot$keggPathId == thisKeggPathId]
	thisPathProtIndexInPca <- as.numeric(na.omit(match(thisPathProtList, ubi.pca.5.scores$spAccession)))
	booleanIndex <- rep(FALSE,nrow(ubi.pca.5.scores))
	booleanIndex[thisPathProtIndexInPca] <- TRUE
	return(booleanIndex)
}


indexList <- lapply(levels(pc.paths.high$keggPathId),FUN=function(x) getPcaProteinIndexFromKeggPath(x))
countListVector <- as.numeric(lapply(indexList,sum))
#countListVector  <- sampleKeggPathCounts$sampleCount[match(levels(as.factor(pc.paths.high$keggPathName)),sampleKeggPathCounts$keggPathName)]
#indexList <- lapply(levels(pc.paths.high$keggPathId),FUN=function(x) (pc.paths.high$keggPathId == x))


lengthPathsHigh <- length(levels(pc.paths.high$keggPathId))

if(output) {
pdf(file="liverProteinsPCvsKeggColoured.pdf", width=10,height=10)
par(mar=c(5,20,5,2))
pathNames <- sampleKeggPathCounts$keggPathName[match(levels(pc.paths.high$keggPathId),sampleKeggPathCounts$keggPathId)]
# get p-values for each term against the remainder
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.1),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(pc.paths.high$Comp.1 ~ pc.paths.high$keggPathId, horizontal=T,las=2, xlab="PC1 (Fetal : Adult)",xlim=c(0,lengthPathsHigh+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.1,add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.1),lty=2)
mtext(countListVector,side=4,at=1:lengthPathsHigh,las=2)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.2),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(pc.paths.high$Comp.2 ~ pc.paths.high$keggPathId, horizontal=T,las=2, xlab="PC2 (HepG2 <-> Fetal)",xlim=c(0,lengthPathsHigh+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.2,add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.2),lty=2)
mtext(countListVector,side=4,at=1:lengthPathsHigh,las=2)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.3),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(pc.paths.high$Comp.3 ~ pc.paths.high$keggPathId, horizontal=T,las=2, xlab="PC3 (Adult culturing)",xlim=c(0,lengthPathsHigh+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.3,add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.3),lty=2)
mtext(countListVector,side=4,at=1:lengthPathsHigh,las=2)
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexGreaterThanRemainderWilcox(abs(ubi.pca.5.scores$Comp.4),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(pc.paths.high$Comp.4 ~ pc.paths.high$keggPathId, horizontal=T,las=2, xlab="PC4 (Highspot <-> Other)",xlim=c(0,lengthPathsHigh+1),col=colVector,names=pathNames) 
boxplot(ubi.pca.5.scores$Comp.4,add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2); abline(v=median(ubi.pca.5.scores$Comp.4),lty=2)
mtext(countListVector,side=4,at=1:lengthPathsHigh,las=2)
dev.off()
}


stopifnot(FALSE)

######NEW stuff

sum(proteinByLiverSample$detectCountAll == 28)
sum(proteinByLiverSample$detectCountFetal == 15)  # all proteins ubiquitous within Fetal are ubiquitious across all samples.


nonHepG2.pca <- princomp(proteinByLiverSample.ubiquitous[,c(fetalIndex,adultIndex)])
summary(nonHepG2.pca)
loadings(nonHepG2.pca)

pdf(file="liverProteinsSampleNonHepG2PCA.pdf", width=10,height=10)
plot(nonHepG2.pca, main="LiverProteinData\nVariance explained by Principal Components")			
biplot(nonHepG2.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(nonHepG2.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(nonHepG2.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(nonHepG2.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(nonHepG2.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(nonHepG2.pca,choices=c(3,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(nonHepG2.pca,choices=c(4,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
dev.off()


##########FETAL PCA
fetal.pca <- princomp(proteinByLiverSample.ubiquitous[,fetalIndex])
summary(fetal.pca)
loadings(fetal.pca)

pdf(file="liverProteinsSampleFetalPCA.pdf", width=10,height=10)
plot(fetal.pca, main="LiverProteinData\nVariance explained by Principal Components")			
biplot(fetal.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(fetal.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(fetal.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(fetal.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(fetal.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(fetal.pca,choices=c(3,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(fetal.pca,choices=c(4,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
dev.off()


##########ADULT PCA
adult.pca <- princomp(proteinByLiverSample.ubiquitous[,adultIndex])
summary(adult.pca)
loadings(adult.pca)

pdf(file="liverProteinsSampleAdultPCA.pdf", width=10,height=10)
plot(adult.pca, main="LiverProteinData\nVariance explained by Principal Components")			
biplot(adult.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(adult.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(adult.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(adult.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(adult.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)
biplot(adult.pca,choices=c(3,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
biplot(adult.pca,choices=c(4,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId)	
dev.off()


############ Development code. 

pdf(file="liverProteinsPCvsKeggColouredAbs.pdf", width=10,height=10)
par(mar=c(5,20,5,2))
# get p-values for each term against the remainder
pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderWilcox(abs(ubi.pca.5.scores$Comp.1),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(abs(pc.paths.high$Comp.1) ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC1 (Fetal : Adult)",xlim=c(0,lengthPathsHigh+1),col=colVector) 
boxplot(abs(ubi.pca.5.scores$Comp.1),add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2)

pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderWilcox(abs(ubi.pca.5.scores$Comp.2),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(abs(pc.paths.high$Comp.2) ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC2 (HepG2 <-> Fetal)",xlim=c(0,lengthPathsHigh+1),col=colVector) 
boxplot(abs(ubi.pca.5.scores$Comp.2),add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2)

pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderWilcox(abs(ubi.pca.5.scores$Comp.3),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(abs(pc.paths.high$Comp.3) ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC3 (Adult culturing)",xlim=c(0,lengthPathsHigh+1),col=colVector) 
boxplot(abs(ubi.pca.5.scores$Comp.3),add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2)

pValList <- as.numeric(lapply(indexList, FUN=function(x) indexVersusRemainderWilcox(abs(ubi.pca.5.scores$Comp.4),x)))
colVector <- ifelse(pValList < 0.05, "Red", "White")
boxplot(abs(pc.paths.high$Comp.4) ~ pc.paths.high$keggPathName, horizontal=T,las=2, xlab="PC4 (Highspot <-> Other)",xlim=c(0,lengthPathsHigh+1),col=colVector) 
boxplot(abs(ubi.pca.5.scores$Comp.4),add=T,at=lengthPathsHigh+1,horizontal=T,col="grey", names=c("All proteins"), las=2)
dev.off()


length(names(all))
head(all[,3:30])
pairs(all[,3:30]))
pairs(all[,3:30])
pairs(all[,3:10])
hist(all$FLM154)
pairs(all[,3:10])
pairs(all[,25:28])
all[,3:31] < 200

pairs(all[,20:28])


all[which(all$AH1 >300 ),]				# some have the value 9999
summary(all)
all[which(all$HepG2.2 >19 ),]

plot(all$FA1,all$AM1)
plot(all$FA1,all$AM1, xlim=c(0,5),ylim=c(0,5))
plot(all$FA1,all$AH1, xlim=c(0,5),ylim=c(0,5))
plot(all$AH1,all$AM1, xlim=c(0,5),ylim=c(0,5))
plot(all$FA1,all$AH1, xlim=c(0,5),ylim=c(0,5))
all[which(all$FLM154 >300 ),]
all[which(all$FLM169 >300 ),]


all$X <- apply(all[,3:30] , 1,max,na.rm=T) 	# find max protein ratio for each row.
hist(all$X)
all
head(all)
head(all)
not_all <- subset(all, X < 9999)			# subset without rows with a 9999  (3 rows)
nrow(not_all)
nrow(all)

not_all <- subset(not_all, X > 0)			# another ~800 rows have too few values to calc max()
nrow(not_all)

pairs(not_all[,3:30])



not_all_no_na <- na.omit(not_all)			# only 431 rows have no NAs.
nrow(not_all_no_na)


bare_matrix <- as.matrix(not_all_no_na[,3:30])

rownames(bare_matrix) <- not_all_no_na$Name.x.x
heatmap(bare_matrix)
## yellow is high, red is low.
# find values for protein with Cytochrome in title.
bare_matrix[grep("Cytochrome", row.names(bare_matrix)),]

pdf("Liver_not_all_no_na.pdf",height=20)
heatmap(bare_matrix)
dev.off()

heatmap(bare_matrix, cexRow=0.3)

png("Liver_not_all_no_na.png",height=2000)
heatmap(bare_matrix)
dev.off()


########## PRINCIPAL COMPONENTS ANALYSIS

pca_test <- princomp(bare_matrix)
pca_test
summary(pca_test)
biplot(pca_test)
loadings(pca_test)  
# in this iteration with 431 trouble-free genes,
# PC1 (43%) separates adult from fetal + HepG2 
# pc2 (12%) separates fetal and HepG2
# pc3 (9%) separates fresh tissue from Matrigel
# 
				
pca_test$scores[ , c(1, 2)]		#access the scores for each protein of the first two principal components.

firstThreeScores <- data.frame(pca_test$scores[ , 1:3])
firstThreeScores[order(firstThreeScores$Comp.1),]	# what's with all the mitochondrial proteins?
row.names(firstThreeScores[order(firstThreeScores$Comp.1),])	#extract row.names ordered

## TO DO: go back and re-work retaining swiss prot IDs with rows
### OR merge it in afterwards?
firstThreeScores$name <- row.names(firstThreeScores)
firstThreeScores <- merge(firstThreeScores,not_all_no_na[,1:2],by.x="name",by.y="Name.x.x")
write.table(firstThreeScores, file="firstThreePcScores.tab",sep="\t",quote=F,row.names=F)


##TO-DO: quantile or rank plots of pairs to look at sample to sample variation in values.
qqplot(all$FLM290,all$FHM290)
abline(a=0,b=1)
# but what does this show?

