#############################################
#
# Dave Gerrard, University of Manchester
#
#############################################


#setwd("C:/Users/dave/LiverProteins/pathways/")

if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

dir()
baseProteinByLiverSample <- read.delim("C:/Users/dave/LiverProteins/data/Cliffs dataset.txt", header=T)

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

# One further manipulation: check which proteins have maxValue of 9999.
proteinByLiverSample$maxValue <- apply(proteinByLiverSample[,allTissueIndex] , 1,max,na.rm=T)  # gives warnings, but not a problem
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
# PC1 (47%) separates adult from (fetal + HepG2) 
# pc2 (13%) separates fetal and HepG2
# pc3 (11%) separates fresh adult tissue from Matrigel
# pc4 (6%) partially separates high-spot from the rest


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



stopifnot(FALSE)

