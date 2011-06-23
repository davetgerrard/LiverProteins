
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

#ubi.pca.5.scores[order(ubi.pca.5.scores$Comp.1),]
#ubi.pca.5.scores[order(ubi.pca.5.scores$Comp.2),]
#ubi.pca.5.scores[order(ubi.pca.5.scores$Comp.3),]
#ubi.pca.5.scores[order(ubi.pca.5.scores$Comp.4),]

if(output) {

## these individual proteins have been concatenated into 'intProts' below.
#plot one protein
thisProt <- "UD14_HUMAN"	# Supposedly Liver specific (UDP-glucuronosyltransferase 1-4). Archtypal PC1 +ve. High PC1 score. ~Zero other PC scores. High in adult, low other. HMGCL_HUMAN is similar
thisProt <- "MARCS_HUMAN"	# Archtypal PC1 -ve. Low PC1 score. ~Zero other PC scores. Low in adult, high other. SERPH_HUMAN is similar
thisProt <- "HBA_HUMAN"		# high PC2 score. Doesn't match fetal+HepG2 grouping. High in Fetal
thisProt <- "UGDH_HUMAN"	# low score on PC2. Doesn't match fetal+HepG2 grouping. High in HepG2
thisProt <- "G6PI_HUMAN"	# 0 for PC1, low PC2 score - High for HepG2, mixed for fetal. Also low PC3 score
thisProt <- "AK1BA_HUMAN"	# low in all fetal and fresh adult, high in cultured adult, med in HepG2. Low PC2 score. High PC4 score.
thisProt <- "AK1C2_HUMAN"	# Archetypal High PC4 score. Specifically high in adult High-spot.
thisProt <- "FIBB_HUMAN"	# fibrinogen. Specifically high in adult High-spot.
thisProt <- "ADH4_HUMAN" 	# high'ish PC1 score. Mostly splits adult vs. (fetal +HepG2). V.low PC3 (low in Matrigel)
thisProt <- "ADH1B_HUMAN"	# similar to ADH4_HUMAN
thisProt <- "ALDOA_HUMAN"	# Low'ish PC1 score. Very low PC2 score. High in HepG2, low in adult, fetal mixed.
thisProt <- "ALBU_HUMAN"	# 0 for PC1, high'ish PC2 score. Mixed in fetal, mod in adult, V. low in HepG2.
thisProt <- "K1C18_HUMAN"	# Keratin, high PC3. High in Matrigel specifically.
thisProt <- "K2C8_HUMAN"	# another Keratin, same profile as K1C18_HUMAN
thisProt <- "VTNC_HUMAN"	# Vitronectin. High PC3. High in Matrigel AND High-spot
thisProt <- "HBE_HUMAN"		# Low PC4 score, low'ish PC1 score. Only high in fresh fetal. 
thisProt <- "PCKGM_HUMAN"	# Low PC4 score. High in fresh adult and matrigel, low/med in high-spot.
thisProt <- "BLVRB_HUMAN"	# Sort of archetypal PC2. Reduced in HepG2 and Matrigel. Flavin reductase.
thisProt <- "ANXA2_HUMAN"	# High PC3.Up in all cultured. CLUS_HUMAN is similar.
thisProt <- "LAMP2_HUMAN"	# High PC3. Rest is a mixed bag.
thisProt <- "VIME_HUMAN"	# Vimentin. High PC4, low PC1. High in fetal high-spot. Relatively higher in adult High-spot.
thisProt <- "MOES_HUMAN"	# Moesin. High'ish PC4 low other. High in High-spot (BOTH ADULT AND FETAL). 

#thisProt <- ""

intProts <- list(UD14_HUMAN="Supposedly Liver specific (UDP-glucuronosyltransferase 1-4). Archtypal PC1 +ve. High PC1 score. ~Zero other PC scores. High in adult, low other. HMGCL_HUMAN is similar",MARCS_HUMAN="Archtypal PC1 -ve. Low PC1 score. ~Zero other PC scores. Low in adult, high other. SERPH_HUMAN is similar",HBA_HUMAN="high PC2 score. Doesn't match fetal+HepG2 grouping. High in Fetal",UGDH_HUMAN="low score on PC2. Doesn't match fetal+HepG2 grouping. High in HepG2",G6PI_HUMAN="0 for PC1, low PC2 score - High for HepG2, mixed for fetal. Also low PC3 score",AK1BA_HUMAN="low in all fetal and fresh adult, high in cultured adult, med in HepG2. Low PC2 score. High PC4 score.",AK1C2_HUMAN="Archetypal High PC4 score. Specifically high in adult High-spot.",FIBB_HUMAN="fibrinogen. Specifically high in adult High-spot.",ADH4_HUMAN="high'ish PC1 score. Mostly splits adult vs. (fetal +HepG2). V.low PC3 (low in Matrigel)",ADH1B_HUMAN="similar to ADH4_HUMAN",ALDOA_HUMAN="Low'ish PC1 score. Very low PC2 score. High in HepG2, low in adult, fetal mixed.",ALBU_HUMAN="0 for PC1, high'ish PC2 score. Mixed in fetal, mod in adult, V. low in HepG2.",K1C18_HUMAN="Keratin, high PC3. High in Matrigel specifically.",K2C8_HUMAN="another Keratin, same profile as K1C18_HUMAN",VTNC_HUMAN="Vitronectin. High PC3. High in Matrigel AND High-spot",HBE_HUMAN="Low PC4 score, low'ish PC1 score. Only high in fresh fetal. ",PCKGM_HUMAN="Low PC4 score. High in fresh adult and matrigel, low/med in high-spot.",BLVRB_HUMAN="Sort of archetypal PC2. Reduced in HepG2 and Matrigel. Flavin reductase.",ANXA2_HUMAN="High PC3.Up in all cultured. CLUS_HUMAN is similar.",LAMP2_HUMAN="High PC3. Rest is a mixed bag.",VIME_HUMAN="Vimentin. High PC4, low PC1. High in fetal high-spot. Relatively higher in adult High-spot.",MOES_HUMAN="Moesin. High'ish PC4 low other. High in High-spot (BOTH ADULT AND FETAL)" )

###################
###should replace with writeLines(strwrap(x, width = 60))  # can it be used with textplot?
###################
utility.insertNewLines <- function(thisText,length=60)  {
	#test if string contains whitespace
	#
	splitLine <- character()
	tokens <- unlist(strsplit(thisText," "))
	while(length(tokens) > 0)  {
		thisLineIndex <- which(cumsum(nchar(tokens)+1) < (length+2) )
		thisLine <- paste(tokens[thisLineIndex],collapse=" ")
		splitLine <- paste(splitLine,thisLine,"\n",sep="")
		tokens <- tokens[-thisLineIndex]

	}
	return(splitLine)
}
####################



proteinReport <- function(thisProt,thisComment="")  {
	#thisComment <- as.character(thisComment)
	layout(matrix(c(1,2,1,3,4,4),ncol=2, byrow=T),heights=c(1,1,4),widths=c(1,2))
	par(mar=c(4,4,4,2))
	barplot(as.matrix(ubi.pca.5.scores[ubi.pca.5.scores$spId == thisProt,1:5]),ylim=c(-3,+3),names.arg=1:5,ylab="PC Score")
	thisProtDesc <- as.character(baseProteinByLiverSample[baseProteinByLiverSample$spId==thisProt,"Name"])
	textplot(paste(thisProtDesc,thisProt,sep="\n"))
	if(nchar(thisComment) > 0) {
		thisComment <- utility.insertNewLines(thisComment,40)
		textplot(thisComment ,halign="left")
	} else {plot.new()}
	par(mar=c(4,9,2,2))
	barplot(log10(as.matrix(baseProteinByLiverSample[baseProteinByLiverSample$spId==thisProt,allTissueBaseIndex])), 
		names.arg=expDesignLiverProteins$niceName[match(colnames(as.matrix(baseProteinByLiverSample.selective[allTissueBaseIndex])),expDesignLiverProteins$sample)],
		horiz=T,las=2)
}

#proteinReport("PCKGM_HUMAN","what?")

pdf("proteinReportsInteresting.pdf")
for(i  in 1:length(intProts))  {
	thisProt <- intProts[i]
	proteinReport(names(thisProt),thisProt[[1]])
}
dev.off()

pdf("proteinReportsUbiq.pdf")
for(thisProt  in ubi.pca.5.scores$spId)  {
	proteinReport(thisProt)
}
dev.off()


##
goTerm <- ""
goTerm <- ""
goTerm <- "GO:0055114"	# oxidation reduction
goTerm <- "GO:0006805"	# xenobiotic metabolic process
goTerm <- "GO:0001889"	# liver development


Ontology(goTerm)[[1]]
Term(goTerm)[[1]]
protIds <- proteinByLiverSample[match(genesInTerm(GOdata,goTerm)[[1]],proteinByLiverSample.ubiquitous$spAccession),"spId"]
fileOutName <- paste("proteinReports",gsub(":","_",goTerm),gsub(" ","_",Term(goTerm)[[1]]),"pdf",sep=".")
pdf(fileOutName)
for(thisProt  in protIds )  {
	proteinReport(thisProt)
}
dev.off()

} # end of if(output)
