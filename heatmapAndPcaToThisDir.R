
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################


setwd("C:/Users/dave/LiverProteins/")

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

##set up name indexes to access functional groups.
freshFetalIndex <- grep('FLM',names(proteinByLiverSample))
highSpotFetalIndex <- grep('FHM',names(proteinByLiverSample))
freshAdultIndex <- grep('FA',names(proteinByLiverSample))
highSpotAdultIndex <- grep('AH',names(proteinByLiverSample))
matrigelAdultIndex <- grep('AM',names(proteinByLiverSample))
fetalIndex <- c(freshFetalIndex,highSpotFetalIndex)
adultIndex <- c(freshAdultIndex,highSpotAdultIndex ,matrigelAdultIndex ) 
hepG2Index <- grep('HepG2',names(proteinByLiverSample))
allTissueIndex <- c(fetalIndex ,adultIndex ,hepG2Index )

# One further manipulation: check which proteins have maxValue of 9999.
proteinByLiverSample$maxValue <- apply(proteinByLiverSample[,allTissueIndex] , 1,max,na.rm=T)
proteinByLiverSample[proteinByLiverSample$maxValue > 20,]
# there are only three none are ubiquitous. Remove from dataset.
proteinByLiverSample <- subset(proteinByLiverSample, maxValue < 21)

#!# Should probably load the other experimental factors like age and MS run.


############ QQ-plots to compare intensity value ranges across samples.
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


## count presence absence across all, across fetal/adult separately and across cultured?
proteinByLiverSample$detectCountAll <- length(allTissueIndex) -  rowSums(is.na(proteinByLiverSample[,allTissueIndex]))
proteinByLiverSample$detectCountFetal <- length(fetalIndex) -  rowSums(is.na(proteinByLiverSample[,fetalIndex]))
proteinByLiverSample$detectCountAdult <- length(adultIndex) -  rowSums(is.na(proteinByLiverSample[,adultIndex]))

# what was I doing here? 
# I want to know if genes missing from some tissues show similarity within a tissue type. 
# Expect high count in fetal with low count in adults and vice versa.
# Don't reallly see it though. May as well show with heatmaps.
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountAdult)
#plot(proteinByLiverSample$detectCountAll, proteinByLiverSample$detectCountFetal)
#plot(proteinByLiverSample$detectCountAdult, proteinByLiverSample$detectCountFetal)


## Make some subsets of data based on presence/absence of proteins across samples

# Subset of data with proteins detected in all samples
proteinByLiverSample.ubiquitous <- subset(proteinByLiverSample, detectCountAll == length(allTissueIndex))
numbUbiquitous  <- nrow(proteinByLiverSample.ubiquitous)

# Subset of data with proteins missing from one or more samples. Convert NA values to zero. 
proteinByLiverSample.selective <- subset(proteinByLiverSample, detectCountAll < length(allTissueIndex))  # <28
proteinByLiverSample.selective <- subset(proteinByLiverSample.selective, detectCountAll > 9)
proteinByLiverSample.selective[is.na(proteinByLiverSample.selective)] <- 0
numbSelective <- nrow(proteinByLiverSample.selective)

## Plot as HEATMAPS.
pdf(file="liverProteinsSampleBasicHeatmaps.pdf", width=15,height=10)
heatmap(as.matrix(proteinByLiverSample.ubiquitous[allTissueIndex]),labRow=proteinByLiverSample.ubiquitous$Name)
mtext(paste(numbUbiquitous,"Ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name)
mtext(paste(numbSelective,"Non-ubiquitous proteins. Euclidean distance"),side=3,line=2,adj=0,cex=1.5)
heatmap(as.matrix(proteinByLiverSample.selective[allTissueIndex]),labRow=proteinByLiverSample.selective$Name,distfun=function(c) dist(c,method="binary"))
mtext(paste(numbSelective,"Non-ubiquitous proteins. Binary distance"),side=3,line=2,adj=0,cex=1.5)
# The tissue clustering is not replicated using the set of proteins shared by 10 to 27 tissues.
dev.off()  # end of "liverProteinsSampleBasicHeatmaps.pdf"


############ Pairwise plots of ubiquitous proteins to compare samples.
## Shows consistent positive correlations within fetal (+ HepG2) and within Adult but negative correlations between adult and fetal.
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

pdf(file="liverProteinsSampleUbiquitousPCA.pdf", width=10,height=10)
plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
plot(ubiquitous.pca, main="LiverProteinData\nVariance explained by Principal Components")			
biplot(ubiquitous.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(ubiquitous.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)	
biplot(ubiquitous.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(3,4), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
biplot(ubiquitous.pca,choices=c(4,5), col=c("grey","black"),cex=c(0.5,1),xlabs=proteinByLiverSample.ubiquitous$spId) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
dev.off()	

## Make lists of genes correlated with each major principal component analysis.

ubi.pca.5.scores <- data.frame(ubiquitous.pca$scores[ , 1:5])
ubi.pca.5.scores$spId <- proteinByLiverSample.ubiquitous$spId
ubi.pca.5.scores$Name <- proteinByLiverSample.ubiquitous$Name
write.table(ubi.pca.5.scores, file="UbiquitousLiverProteinsFirst5PcaScores.tab",sep="\t",row.names=F,quote=F)


#!# NETWORK ANALYSIS






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

