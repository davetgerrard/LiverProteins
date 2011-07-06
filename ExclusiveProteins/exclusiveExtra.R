
setwd("C:/Users/dave/LiverProteins/ExclusiveProteins/")


####################
####################

## load the ITRAQ data
baseProteinByLiverSample <- read.delim("C:/Users/dave/LiverProteins/data/Cliffs dataset.txt", header=T)
## load the other experimental factors like age and MS run.
expDesignLiverProteins <- read.delim("C:/Users/dave/LiverProteins/data/CliffsExpDesign.txt", header=T)

summary(baseProteinByLiverSample)
nrow(baseProteinByLiverSample)

# There are 2304 rows and many NAs (missing values). 
# Only a small number of proteins are present in all samples.
# Some samples have a maximum value of 9999. 

############ MANIPULATIONS to clean and index data.
#!# Probably should split the Accession into three
baseProteinByLiverSample$spAccession <- matrix(unlist(strsplit(as.character(baseProteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]
baseProteinByLiverSample$spId <- matrix(unlist(strsplit(as.character(baseProteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,3]


# The second column has an odd name
names(baseProteinByLiverSample)[grep('Name.x.x',names(baseProteinByLiverSample))] <- "Name"
# One column is FL257 instead of FLM257
names(baseProteinByLiverSample)[grep('FL257',names(baseProteinByLiverSample))] <- "FLM257"
# The final column is not data. 
baseProteinByLiverSample <- subset(baseProteinByLiverSample, select= -X)

# make a copy of the dataset BEFORE making counts across indices.
proteinByLiverSample <- baseProteinByLiverSample


## Set up indices to base data

freshFetalBaseIndex <- grep('FLM',names(baseProteinByLiverSample))
highSpotFetalBaseIndex <- grep('FHM',names(baseProteinByLiverSample))
freshAdultBaseIndex <- grep('FA',names(baseProteinByLiverSample))
highSpotAdultBaseIndex <- grep('AH',names(baseProteinByLiverSample))
matrigelAdultBaseIndex <- grep('AM',names(baseProteinByLiverSample))
fetalBaseIndex <- c(freshFetalBaseIndex,highSpotFetalBaseIndex)	#excluded HepG2
adultBaseIndex <- c(freshAdultBaseIndex,highSpotAdultBaseIndex ,matrigelAdultBaseIndex ) 
hepG2BaseIndex <- grep('HepG2',names(baseProteinByLiverSample))
allTissueBaseIndex <- c(fetalBaseIndex ,adultBaseIndex ,hepG2BaseIndex )

## count presence absence across all, across fetal/adult separately and across cultured?
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

## make detection counts. Do this separately for base data and foreground data.

for(thisExp in levels(expDesignLiverProteins$expRun))  {
	expSamples <- as.character(expDesignLiverProteins$sample[expDesignLiverProteins$expRun== thisExp])
	expBaseIndices[[thisExp]] <- match(expSamples,names(baseProteinByLiverSample))
	baseProteinByLiverSample[,paste(thisExp,"Count",sep="")] <- length(expBaseIndices[[thisExp]]) -  rowSums(is.na(baseProteinByLiverSample[,expBaseIndices[[thisExp]]]))
}





##### strart cleaning the copy of the data we are going to use

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

expIndices <- as.list(paste(levels(expDesignLiverProteins$expRun),"Index",sep=""))
names(expIndices) <- levels(expDesignLiverProteins$expRun)

## make detection counts. Do this separately for base data and foreground data.

for(thisExp in levels(expDesignLiverProteins$expRun))  {
	expSamples <- as.character(expDesignLiverProteins$sample[expDesignLiverProteins$expRun== thisExp])
	expIndices[[thisExp]] <- as.integer(na.omit(match(expSamples,names(proteinByLiverSample))))
	proteinByLiverSample[,paste(thisExp,"Count",sep="")] <- length(expIndices[[thisExp]]) -  rowSums(is.na(proteinByLiverSample[,expIndices[[thisExp]]]))
}



# One further manipulation: check which proteins have maxValue of 9999.
proteinByLiverSample$maxValue <- apply(proteinByLiverSample[,allTissueIndex] , 1,max,na.rm=T)  # gives warnings, but not a problem
proteinByLiverSample[proteinByLiverSample$maxValue > 20,]
# there are only three none are ubiquitous. Remove from dataset.
proteinByLiverSample <- subset(proteinByLiverSample, maxValue < 21)





####################
####################






## additional searches for 'exclusive' proteins.

### MAY DEPEND ON INCLUSION OF SAMPLES 290, 292



##### base data counts.
nrow(baseProteinByLiverSample[baseProteinByLiverSample$baseDetectCountFetal == 0,adultBaseIndex])


## look at counts of proteins across experimental runs. (esp. exp4)
# 493 proteins are detected in all 7 samples from a single run but only those samples.
nrow(subset(baseProteinByLiverSample, exp1Count ==7 & baseDetectCountAll  ==7))
nrow(subset(baseProteinByLiverSample, exp2Count ==7 & baseDetectCountAll  ==7))
nrow(subset(baseProteinByLiverSample, exp3Count ==7 & baseDetectCountAll  ==7))
nrow(subset(baseProteinByLiverSample, exp4Count ==7 & baseDetectCountAll  ==7))
#altogether
nrow(subset(baseProteinByLiverSample, (exp1Count ==7 |exp2Count ==7 | exp3Count ==7 | exp4Count ==7) & baseDetectCountAll  ==7))
##### foreground (non-base) data


baseProteins.ubi <- baseProteinByLiverSample$spAccession[baseProteinByLiverSample$baseDetectCountAll ==28]
baseProteins.ExpSpecific <- subset(baseProteinByLiverSample, (exp1Count ==7 |exp2Count ==7 | exp3Count ==7 | exp4Count ==7) & baseDetectCountAll  ==7)$spAccession

baseProteins.remainder <- setdiff(baseProteinByLiverSample$spAccession,c(baseProteins.ubi,baseProteins.ExpSpecific))

#intersect(baseProteins.remainder,baseProteins.ubi)	# test if remainder is correct

head(baseProteinByLiverSample[match(baseProteins.remainder,baseProteinByLiverSample$spAccession),])

proteinByLiverSample.remainder <- baseProteinByLiverSample[match(baseProteins.remainder,baseProteinByLiverSample$spAccession),]

table(proteinByLiverSample.remainder$baseDetectCountAll)  
# three common values, 0, 14, 21. 
# How can it be zero?
# 791 proteins have NA for ALL samples. (this is correct in original dataset).
head(proteinByLiverSample.remainder[proteinByLiverSample.remainder$baseDetectCountAll == 0,])
# 207 proteins are detected in all samples of two experiments and none of the samples from the other two.
# There are only 7s in the expCount values in the following tables.
table(proteinByLiverSample.remainder[proteinByLiverSample.remainder$baseDetectCountAll == 14,c("exp1Count","exp2Count","exp3Count","exp4Count")])
# Another 207 proteins are detected  in all samples of three experiments and none of the samples of the fourth experiment.
# The table has 209 rows, there are two exceptions
table(proteinByLiverSample.remainder[proteinByLiverSample.remainder$baseDetectCountAll == 21,c("exp1Count","exp2Count","exp3Count","exp4Count")])





432+493+791+207+207 # = 2130
nrow(baseProteinByLiverSample)
# leaves 174 proteins present in 'some' samples but not all, and not perfectly partitioned by experimental run.


# what of these 174 disparate proteins?
proteinByLiverSample.disparate <- subset(proteinByLiverSample.remainder, (baseDetectCountAll != 0 & baseDetectCountAll != 14 & baseDetectCountAll != 21))
hist(proteinByLiverSample.disparate$baseDetectCountAll,breaks=29)	# missing 0,7,14,21,28 as it should. Slight peak at 6 but otherwise nothing surprising.

##### what about proteins missing from all Fetal samples

baseProteinByLiverSample[baseProteinByLiverSample$baseDetectCountFetal == 0,adultIndex]   ## lots of NAs
proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountFetal == 0,adultBaseIndex]  ## Only obvious pattern is FA4,AH4, AM4 which were all exp4
proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountAdult == 0,fetalBaseIndex ] 

proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountFetal == 0,hepG2BaseIndex ]
proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountAdult == 0,hepG2BaseIndex ]

## for the disparate proteins with no detectinon in fetal, they mostly still cluster by experiment run. e.g. 3/3
proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountFetal == 0, ]
## same shit for disparate proteins with no detection in adult. Where they exist in more than one sample, all those samples are from the same exp run (max =3)
proteinByLiverSample.disparate[proteinByLiverSample.disparate$baseDetectCountAdult == 0, ]

proteinByLiverSample.disparate[,freshFetalBaseIndex ]
write.table(proteinByLiverSample.disparate,file="disparateProteins.tab",sep="\t",row.names=F,quote=F)

plot(proteinByLiverSample.disparate$baseDetectCountFetal,proteinByLiverSample.disparate$baseDetectCountAdult)
cor.test(proteinByLiverSample.disparate$baseDetectCountFetal,proteinByLiverSample.disparate$baseDetectCountAdult)
subset(proteinByLiverSample.disparate,baseDetectCountFetal >6 & baseDetectCountAdult < 4)
subset(proteinByLiverSample.disparate,baseDetectCountFetal < 4 & baseDetectCountAdult > 4)
subset(proteinByLiverSample.disparate,baseDetectCountFetal ==0 & baseDetectCountAdult > 1)
subset(proteinByLiverSample.disparate,baseDetectCountFetal > 1 & baseDetectCountAdult == 0 )

subset(proteinByLiverSample.disparate, baseDetectCountAdult == 0 )

baseProteinByLiverSample.selective <- subset(baseProteinByLiverSample, baseDetectCountAll < length(allTissueBaseIndex) & baseDetectCountAll > 2)  # <28
baseProteinByLiverSample.selective[is.na(baseProteinByLiverSample.selective)] <- 0
numbBaseSelective <- nrow(baseProteinByLiverSample.selective)

niceNames <- expDesignLiverProteins$niceName[match(colnames(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex])),expDesignLiverProteins$sample)]
expLabels <- expDesignLiverProteins$expRun[match(colnames(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex])),expDesignLiverProteins$sample)]
combLabels <-  paste(niceNames,expLabels,sep=" : ")
heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow=baseProteinByLiverSample.selective$Name)
mtext(paste(numbBaseSelective,"Non-ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
par(mar=c(7,3,3,2))
#heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow="",labCol=combLabels ,distfun=function(c) dist(c,method="binary"),col=brewer.pal(3,"Greys"),mar=c(5,1)) # grey version
heatmap(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex]),labRow="",labCol=combLabels ,distfun=function(c) dist(c,method="binary"),col=heat.colors(256),mar=c(5,1))
mtext(paste(numbBaseSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=1,adj=0,cex=1.5)
#mtext(expLabels,side=1,line=6,las=3,at=1:length(expLabels))


#################older stuff




#missing from one but in some of the opposite/ancestral

#not fresh fetal but in some adult
proteinByLiverSample[proteinByLiverSample$detectCountFetal == 0,]
proteinByLiverSample[proteinByLiverSample$detectCountFetal == 0,adultIndex]

# there are 959 proteins not detected in fetal tissue
nrow(proteinByLiverSample[proteinByLiverSample$detectCountFetal == 0,adultIndex])
# 123 of these are detected in 3+ other tissues
nrow(subset(proteinByLiverSample,detectCountFetal == 0  & detectCountAll >= 3))
# 121 of these have exactly 3 hits in adult
nrow(subset(proteinByLiverSample,detectCountFetal == 0  & detectCountAdult == 3))
# 119 of these appear in FA4, AH4 and AM4 samples.
# These 3 samples were all run together in i-TRAC with the 290 and 292 samples.

baseProteinByLiverSample$spAccession <- matrix(unlist(strsplit(as.character(baseProteinByLiverSample$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]
exp4prots <- subset(proteinByLiverSample,detectCountFetal == 0  & detectCountAdult == 3)$spAccession
exp4prot.baseIndex <- match(exp4prots, baseProteinByLiverSample$spAccession)
(baseProteinByLiverSample[exp4prot.baseIndex,])
# all but ~12 of these proteins are also detected in samples FLM290, FLM292, FHM290 & FHM292


#not fresh adult but in some fetal

#not fresh fetal but in some cultured fetal

#not fresh adult but in some cultured adult

#already tested for hepG2 exclusive, 
#what about shared with adult/fetal and missing from opposite?