
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################


#######INFO: LOAD & CLEAN PROTEIN DATA

#setwd("C:/Users/dave/LiverProteins/pathways/")

if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

dir()
##INFO: load the ITRAQ data
baseProteinByLiverSample <- read.delim("C:/Users/dave/LiverProteins/data/Cliffs dataset.txt", header=T)
##INFO: load the other experimental factors like age and MS run.
expDesignLiverProteins <- read.delim("C:/Users/dave/LiverProteins/data/CliffsExpDesign.txt", header=T)

summary(baseProteinByLiverSample)
nrow(baseProteinByLiverSample)

##INFO: There are 2304 rows and many NAs (missing values). 
##INFO: Only a small number of proteins are present in all samples.
##INFO: Some samples have a maximum value of 9999. 

##INFO: MANIPULATIONS to clean and index data.
#!# Probably should split the Accession into three
baseProteinByLiverSample$spAccession <- matrix(unlist(strsplit(as.character(baseProteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]
baseProteinByLiverSample$spId <- matrix(unlist(strsplit(as.character(baseProteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,3]


#INFO: The second column has an odd name
names(baseProteinByLiverSample)[grep('Name.x.x',names(baseProteinByLiverSample))] <- "Name"
#INFO: One column is FL257 instead of FLM257
names(baseProteinByLiverSample)[grep('FL257',names(baseProteinByLiverSample))] <- "FLM257"
#INFO: The final column is not data. 
baseProteinByLiverSample <- subset(baseProteinByLiverSample, select= -X)

#INFO: make a copy of the dataset BEFORE making counts across indices.
proteinByLiverSample <- baseProteinByLiverSample


#INFO: Set up indices to base data

freshFetalBaseIndex <- grep('FLM',names(baseProteinByLiverSample))
highSpotFetalBaseIndex <- grep('FHM',names(baseProteinByLiverSample))
freshAdultBaseIndex <- grep('FA',names(baseProteinByLiverSample))
highSpotAdultBaseIndex <- grep('AH',names(baseProteinByLiverSample))
matrigelAdultBaseIndex <- grep('AM',names(baseProteinByLiverSample))
fetalBaseIndex <- c(freshFetalBaseIndex,highSpotFetalBaseIndex)	#excluded HepG2
adultBaseIndex <- c(freshAdultBaseIndex,highSpotAdultBaseIndex ,matrigelAdultBaseIndex ) 
hepG2BaseIndex <- grep('HepG2',names(baseProteinByLiverSample))
allTissueBaseIndex <- c(fetalBaseIndex ,adultBaseIndex ,hepG2BaseIndex )

#INFO: count presence absence across all, across fetal/adult separately and across cultured?
baseProteinByLiverSample$baseDetectCountAll <- length(allTissueBaseIndex) -  rowSums(is.na(baseProteinByLiverSample[,allTissueBaseIndex]))
baseProteinByLiverSample$baseDetectCountFetal <- length(fetalBaseIndex) -  rowSums(is.na(baseProteinByLiverSample[,fetalBaseIndex]))
baseProteinByLiverSample$baseDetectCountAdult <- length(adultBaseIndex) -  rowSums(is.na(baseProteinByLiverSample[,adultBaseIndex]))
baseProteinByLiverSample$baseDetectCountHepG2 <- length(hepG2BaseIndex) -  rowSums(is.na(baseProteinByLiverSample[,hepG2BaseIndex]))
baseProteinByLiverSample$baseDetectCountAM <- length(matrigelAdultBaseIndex) -  rowSums(is.na(baseProteinByLiverSample[,matrigelAdultBaseIndex]))
baseProteinByLiverSample$baseDetectCountAH <- length(highSpotAdultBaseIndex ) -  rowSums(is.na(baseProteinByLiverSample[,highSpotAdultBaseIndex ]))
baseProteinByLiverSample$baseDetectCountFHM <- length(highSpotFetalBaseIndex ) -  rowSums(is.na(baseProteinByLiverSample[,highSpotFetalBaseIndex ]))
baseProteinByLiverSample$baseDetectCountFLM <- length(freshFetalBaseIndex ) -  rowSums(is.na(baseProteinByLiverSample[,freshFetalBaseIndex ]))
baseProteinByLiverSample$baseDetectCountFA <- length(freshAdultBaseIndex ) -  rowSums(is.na(baseProteinByLiverSample[,freshAdultBaseIndex ]))

expBaseIndices <- as.list(paste(levels(expDesignLiverProteins$expRun),"BaseIndex",sep=""))
names(expBaseIndices ) <- levels(expDesignLiverProteins$expRun)

#INFO: make detection counts. Do this separately for base data and foreground data.

for(thisExp in levels(expDesignLiverProteins$expRun))  {
	expSamples <- as.character(expDesignLiverProteins$sample[expDesignLiverProteins$expRun== thisExp])
	expBaseIndices[[thisExp]] <- match(expSamples,names(baseProteinByLiverSample))
	baseProteinByLiverSample[,paste(thisExp,"Count",sep="")] <- length(expBaseIndices[[thisExp]]) -  rowSums(is.na(baseProteinByLiverSample[,expBaseIndices[[thisExp]]]))
}





#INFO: start cleaning the copy of the data we are going to use

##
#INFO: THIS VERSION: REMOVE SAMPLES 290/292
##
proteinByLiverSample <- subset(proteinByLiverSample, select=c(-FLM290,-FHM290,-FLM292,-FHM292))


#INFO: set up name indexes to access functional groups.
freshFetalIndex <- grep('FLM',names(proteinByLiverSample))
highSpotFetalIndex <- grep('FHM',names(proteinByLiverSample))
freshAdultIndex <- grep('FA',names(proteinByLiverSample))
highSpotAdultIndex <- grep('AH',names(proteinByLiverSample))
matrigelAdultIndex <- grep('AM',names(proteinByLiverSample))
fetalIndex <- c(freshFetalIndex,highSpotFetalIndex)	#excluded HepG2
adultIndex <- c(freshAdultIndex,highSpotAdultIndex ,matrigelAdultIndex ) 
hepG2Index <- grep('HepG2',names(proteinByLiverSample))
allTissueIndex <- c(fetalIndex ,adultIndex ,hepG2Index )
freshTissueIndex <- c(freshFetalIndex, freshAdultIndex, hepG2Index)

expIndices <- as.list(paste(levels(expDesignLiverProteins$expRun),"Index",sep=""))
names(expIndices) <- levels(expDesignLiverProteins$expRun)

#INFO: make detection counts. Do this separately for base data and foreground data.

for(thisExp in levels(expDesignLiverProteins$expRun))  {
	expSamples <- as.character(expDesignLiverProteins$sample[expDesignLiverProteins$expRun== thisExp])
	expIndices[[thisExp]] <- as.integer(na.omit(match(expSamples,names(proteinByLiverSample))))
	proteinByLiverSample[,paste(thisExp,"Count",sep="")] <- length(expIndices[[thisExp]]) -  rowSums(is.na(proteinByLiverSample[,expIndices[[thisExp]]]))
}



#INFO: count presence absence across all, across fetal/adult separately and across cultured?
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


#INFO: One further manipulation: check which proteins have maxValue of 9999.
proteinByLiverSample$maxValue <- apply(proteinByLiverSample[,allTissueIndex] , 1,max,na.rm=T)  # gives warnings, but not a problem
proteinByLiverSample[proteinByLiverSample$maxValue > 20,]
#INFO: there are only three none are ubiquitous. Remove from dataset.
proteinByLiverSample <- subset(proteinByLiverSample, maxValue < 21)

############INFO: QQ-plots to compare intensity value ranges across samples.
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

#INFO:  Make some subsets of data based on presence/absence of proteins across samples

#INFO:  Subset of data with proteins detected in all samples
proteinByLiverSample.ubiquitous <- subset(proteinByLiverSample, detectCountAll == length(allTissueIndex))
numbUbiquitous  <- nrow(proteinByLiverSample.ubiquitous)

#INFO:  Subset of data with proteins missing from one or more samples. Convert NA values to zero. 
proteinByLiverSample.selective <- subset(proteinByLiverSample, detectCountAll < length(allTissueIndex) & detectCountAll > 2)  # <24
#proteinByLiverSample.selective <- subset(proteinByLiverSample.selective, detectCountAll > 2)
proteinByLiverSample.selective[is.na(proteinByLiverSample.selective)] <- 0
numbSelective <- nrow(proteinByLiverSample.selective)

#INFO:  Subset of BASE data (all samples) with proteins missing from some samples.
baseProteinByLiverSample.selective <- subset(baseProteinByLiverSample, baseDetectCountAll < length(allTissueBaseIndex) & baseDetectCountAll > 2)  # <28
baseProteinByLiverSample.selective[is.na(baseProteinByLiverSample.selective)] <- 0
numbBaseSelective <- nrow(baseProteinByLiverSample.selective)
niceNames <- expDesignLiverProteins$niceName[match(colnames(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex])),expDesignLiverProteins$sample)]
expLabels <- expDesignLiverProteins$expRun[match(colnames(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex])),expDesignLiverProteins$sample)]
combLabels <-  paste(niceNames,expLabels,sep=" : ")
par(mar=c(7,3,3,2))



#INFO:  Plot as Heatmaps.
if(output) {
	pdf(file="liverProteinsSampleBasicHeatmaps.pdf", width=15,height=10)
	heatmap(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex]),labRow=proteinByLiverSample.ubiquitous$Name)
	mtext(paste(numbUbiquitous,"Ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
	heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name)
	mtext(paste(numbSelective,"Non-ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
	heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name,distfun=function(c) dist(c,method="binary"))
	mtext(paste(numbSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=2,adj=0,cex=1.5)
	#INFO:  The tissue clustering is not replicated using the set of proteins shared by 10 to 27 tissues.
	#heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow="",labCol=combLabels ,distfun=function(c) dist(c,method="binary"),col=brewer.pal(3,"Greys"),mar=c(5,1)) # grey version
	par(mar=c(7,3,3,2))
	heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow="",labCol=combLabels ,distfun=function(c) dist(c,method="binary"),col=heat.colors(256),margins=c(10,1))
	mtext(paste(numbBaseSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=1,adj=0,cex=1.5)
	dev.off()  # end of "liverProteinsSampleBasicHeatmaps.pdf"

	##### heatmaps must be output separately. Can be joined in GIMP.
	thisFigure <- "3"   # ~ 15 samples
	tiff(paste("fig",thisFigure,"tif",sep="."),compression="lzw",width=180, height=120,units="mm",res=300)
	niceNames.fresh <- expDesignLiverProteins$niceName[match(colnames(as.matrix(proteinByLiverSample.ubiquitous[freshTissueIndex])),expDesignLiverProteins$sample)]
	expLabels.fresh <- expDesignLiverProteins$expRun[match(colnames(as.matrix(proteinByLiverSample.ubiquitous[freshTissueIndex])),expDesignLiverProteins$sample)]
	combLabels.fresh <-  paste(niceNames.fresh, expLabels.fresh,sep=" : ")
	test_heat  <- heatmap(t(as.matrix(proteinByLiverSample.ubiquitous[freshTissueIndex])),labCol="",labRow=niceNames.fresh,
				col=heat.colors(256),margins=c(1,8))
	mtext(thisFigure ,side=3,line=1,adj=0,cex=1.5)
	#mtext("1A: Ubiquitous proteins. Fresh adult, fetal and HepG2",side=3,line=1,cex=1.5)
	#plot(test_heat)
	dev.off()

	thisFigure <- "5"	# ~ 24 samples
	tiff(paste("fig",thisFigure,"tif",sep="."),compression="lzw",width=180, height=180,units="mm",res=300)
	niceNames.ubiq <- expDesignLiverProteins$niceName[match(colnames(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex])),expDesignLiverProteins$sample)]
	expLabels.ubiq <- expDesignLiverProteins$expRun[match(colnames(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex])),expDesignLiverProteins$sample)]
	combLabels.ubiq <-  paste(niceNames.ubiq, expLabels.ubiq,sep=" : ")
	test_heat  <- heatmap(t(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex])),labCol="",labRow=niceNames.ubiq,
				col=heat.colors(256),margins=c(1,8))
	mtext(thisFigure ,side=3,line=1,adj=0,cex=1.5)
	#mtext("dtg.1B: Ubiquitous proteins. All samples",side=3,line=1,cex=1.5)
	#plot(test_heat)
	dev.off()

	thisFigure <- "S1"		# some weirdness happens when I try and rotate this heatmap using t(). All goes red.
	tiff(paste("fig",thisFigure,"tif",sep="."),compression="lzw",width=180, height=180,units="mm",res=300)
	heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow="",labCol=combLabels ,distfun=function(c) dist(c,method="binary"),col=heat.colors(256),margins=c(10,1))
	mtext(thisFigure ,side=3,line=1,adj=0,cex=1.5)
	#mtext(paste(numbBaseSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=1,adj=0,cex=1.5)
	dev.off()


} #end of output if()




############INFO:  Pairwise plots of ubiquitous proteins to compare samples.
##INFO:  Shows consistent positive correlations within fetal (+ HepG2) and within Adult but negative correlations between adult and fetal.
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

##########INFO: PRINCIPAL COMPONENTS ANALYSIS

ubiquitous.pca <- princomp(proteinByLiverSample.ubiquitous[,allTissueIndex])
ubiquitous.pca
summary(ubiquitous.pca)

loadings(ubiquitous.pca)  
#INFO: PCA on with 431 ubiquitous proteins,
#INFO: PC1 (47%) separates adult from (fetal + HepG2) 
#INFO: pc2 (13%) separates fetal and HepG2
#INFO: pc3 (11%) separates fresh adult tissue from Matrigel
#INFO: pc4 (6%) partially separates high-spot from the rest


#INFO:  set up some symbols and numbers for biplots
sampleTypes <- c("Adult fresh", "Adult ALI-3D", "Adult Matrigel","Fetal fresh", "Fetal ALI-3D", "HepG2")
samplePch <- c(15,0,12,16,1,10,17)	# squares for adult, circles for fetal. Solid, empty, crossed for fresh, ALI-3D, Matrigel. 

expDesignLiverProteins$samplePch <- 0
for(i in 1:length(sampleTypes)) {
	thisSample <- sampleTypes[i]
	expDesignLiverProteins$samplePch[grep(thisSample, expDesignLiverProteins$niceName)] <- samplePch[i]
}

pchNames.ubiq <- expDesignLiverProteins$samplePch[match(colnames(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex])),expDesignLiverProteins$sample)]

#INFO: Graph  loadings and pc biplots -> "liverProteinsSampleUbiquitousPCA.pdf"
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

#par(mfrow=c(2,2))
#plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
#biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
#biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
#biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId,xlim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

#par(mfrow=c(2,2))
#plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
#biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
#biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
#biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

#par(mfrow=c(2,2))
#star.names <- rep("*",nrow(proteinByLiverSample.ubiquitous))
#plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
#biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
#biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
#biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

#par(mfrow=c(2,2))
#star.names <- rep("*",nrow(proteinByLiverSample.ubiquitous))
#plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
#biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
#biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
#biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=star.names, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

#tiff("fig.dtg.2.tiff",compression="lzw",width=180, height=180,units="mm",res=300)
#op <- par(mfrow=c(2,2),mar=c(4,4,4,2))
#star.names <- rep("*",nrow(proteinByLiverSample.ubiquitous))
#plot(ubiquitous.pca, main="") ; mtext("A" ,side=3,line=2,adj=0,cex=1.5)			
#biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,0.6),xlabs=star.names, arrow.len=0, ylabs=niceNames.ubiq, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) ; mtext("B" ,side=3,line=2,adj=0,cex=1.5)
#biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,0.6),xlabs=star.names, arrow.len=0, ylabs=niceNames.ubiq, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) ; mtext("C" ,side=3,line=2,adj=0,cex=1.5)	
#biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,0.6),xlabs=star.names, arrow.len=0, ylabs=niceNames.ubiq, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125)) ;abline(v=0,lty=2) ; abline(h=0,lty=2) ; mtext("D" ,side=3,line=2,adj=0,cex=1.5)
#par(op)
#dev.off()


##hacked into biplot to find out how to generate co-ordinates for arrows and use pch points.
#stats:::biplot.princomp
#

#ubiquitous.pca$loadings[,1:2]
#INFO: Graph selected loadings and pc biplots -> "fig.6.tiff"
tiff("fig.6.tiff",compression="lzw",width=180, height=180,units="mm",res=600)
op <- par(mfrow=c(2,2),mar=c(4,4,4,2))
star.names <- rep("*",nrow(proteinByLiverSample.ubiquitous))
blankSampleNames <- rep("",length(niceNames.ubiq))
# see stats:::print.summary.princomp for calculation of variance proportions.
prop.vars <- ubiquitous.pca$sdev^2/sum(ubiquitous.pca$sdev^2)
plot(ubiquitous.pca, main="") ; mtext("A" ,side=3,line=2,adj=0,cex=1.5)		
legend("topright", legend=sampleTypes, pch=samplePch) 
#thisChoice <- c(1,2)
#thisChoice <- c(1,4)
for(i in 2:4) {
	thisChoice <- c(1,i)
	xlabel <- paste(names(prop.vars)[thisChoice[1]]," (",round(prop.vars[thisChoice[1]]*100),"%)",sep="")
	ylabel <- paste(names(prop.vars)[thisChoice[2]]," (",round(prop.vars[thisChoice[2]]*100),"%)",sep="")
	biplot(ubiquitous.pca,choices=thisChoice, col=c("grey","black"),cex=c(0.5,0.6),xlab=xlabel, ylab=ylabel,xlabs=star.names, ylabs=blankSampleNames, xlim=c(-0.125,0.125), ylim=c(-0.125,0.125), var.axes=F) ;abline(v=0,lty=2) ; abline(h=0,lty=2) ; mtext(LETTERS[i] ,side=3,line=2,adj=0,cex=1.5)
	tempLam <- ubiquitous.pca$sdev[thisChoice] * sqrt(ubiquitous.pca$n.obs)
	points(t(t(ubiquitous.pca$loadings[,thisChoice]) * tempLam), pch=pchNames.ubiq, cex=1.5 )
	#legend("topleft", legend=sampleTypes, pch=samplePch, cex=0.8) 
}
par(op)
dev.off()

}
## Make lists of genes correlated with each major principal component analysis.

#INFO: Write protein pc scores as table -> "UbiquitousLiverProteinsFirst5PcaScores.tab"
ubi.pca.5.scores <- data.frame(ubiquitous.pca$scores[ , 1:5])
ubi.pca.5.scores$spId <- proteinByLiverSample.ubiquitous$spId
ubi.pca.5.scores$Name <- proteinByLiverSample.ubiquitous$Name
ubi.pca.5.scores$spAccession <- proteinByLiverSample.ubiquitous$spAccession
if(output) {
	write.table(ubi.pca.5.scores, file="UbiquitousLiverProteinsFirst5PcaScores.tab",sep="\t",row.names=F,quote=F)
}
proteinByLiverSample.ubi.pca <- merge(proteinByLiverSample.ubiquitous,ubi.pca.5.scores,by="spAccession")








stopifnot(FALSE)


