
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################

##INFO: Plot venn diagrams comparing Protein detection for Adult, Fetal and HepG2 and, seperately for the 4 Experimental groups.

require(gplots)


## useful columns in proteinByLiverSample:- 
# "exp1Count"        "exp2Count"        "exp3Count"        "exp4Count"        "detectCountAll"   "detectCountFetal" "detectCountAdult"
# "detectCountHepG2" "detectCountAM"    "detectCountAH"    "detectCountFHM"   "detectCountFLM"   "detectCountFA"


#### N.B. In the following, ALL groups means ALL in one group and not ALL in other groups.
#### What I actually want is ALL in one group and not ANY in other groups.


if(output)  {

	pdf(file="liverProteinVennDetected.pdf")
	textWidth <- 30
	# this section is outputting extra blank pages. I think it's something to do with venn but can't work out what. 
	section.text <- "Simple venn diagrams of sample type for proteins present in ALL or ANY samples of a group"
	textplot(strwrap(section.text , width =textWidth ))
	
	grouping.typeAll <- "HepG2, Fetal, Adult. Detected in ALL of each."
	list.types <- c("Adult", "Fetal", "HepG2")
	input.typeAll <- list()
	for(thisType in list.types)  {
		colName <- paste("detectCount",thisType,sep="")
		input.typeAll[[thisType]] <- proteinByLiverSample$spAccession[proteinByLiverSample[,colName] == max(proteinByLiverSample[,colName])]
	}
	#plot.new()  # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.

	venn(input.typeAll)
	mtext(grouping.typeAll, side=1, cex=1.2)

	grouping.typeAny <- "HepG2, Fetal, Adult. Detected in ANY ONE of each."
	list.types <- c("Adult", "Fetal", "HepG2")
	input.typeAny <- list()
	for(thisType in list.types)  {
		colName <- paste("detectCount",thisType,sep="")
		input.typeAny[[thisType]] <- proteinByLiverSample$spAccession[proteinByLiverSample[,colName] > 0]
	}
	if(names(dev.cur()) =="null device") {plot.new() } # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.typeAny)
	mtext(grouping.typeAny, side=1, cex=1.2)

	section.text <- "Complex venn diagrams of sample type for proteins present in ALL of one sample group and ANY samples of the other groups"
	textplot(strwrap(section.text , width =textWidth ))


	grouping.HepG2.Exclusive <- "ALL HepG2 vs ANY in Fetal, Adult."
	list.types <- c("Adult", "Fetal", "HepG2")
	input.HepG2.Exclusive <- list(Adult=input.typeAny[["Adult"]], Fetal=input.typeAny[["Fetal"]], HepG2=input.typeAll[["HepG2"]])
	plot.new()  # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.HepG2.Exclusive)
	mtext(grouping.HepG2.Exclusive, side=1, cex=1.2)
	
	grouping.Fetal.Exclusive <- "ALL Fetal vs ANY in HepG2, Adult."
	list.types <- c("Adult", "Fetal", "HepG2")
	input.Fetal.Exclusive <- list(Adult=input.typeAny[["Adult"]], Fetal=input.typeAll[["Fetal"]], HepG2=input.typeAny[["HepG2"]])
	plot.new()  # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.Fetal.Exclusive)
	mtext(grouping.Fetal.Exclusive, side=1, cex=1.2)

	grouping.Adult.Exclusive <- "ALL Adult vs ANY in Fetal, HepG2."
	list.types <- c("Adult", "Fetal", "HepG2")
	input.Adult.Exclusive <- list(Adult=input.typeAll[["Adult"]], Fetal=input.typeAny[["Fetal"]], HepG2=input.typeAny[["HepG2"]])
	plot.new()  # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.Adult.Exclusive)
	mtext(grouping.Adult.Exclusive, side=1, cex=1.2)


	#### Venn on Experimental runs
	section.text <- "Venn diagrams of proteins present in ANY and ALL samples within an experiment group"
	textplot(strwrap(section.text , width =textWidth ))

	grouping.expAll <- "The Four experiment groups. Detected in ALL of each."
	input.expAll <- list()
	for(thisExp in levels(expDesignLiverProteins$expRun))  {
		colName <- paste(thisExp,"Count",sep="")
		input.expAll[[thisExp]] <- proteinByLiverSample$spAccession[proteinByLiverSample[,colName] == max(proteinByLiverSample[,colName])]
	}
	if(names(dev.cur()) =="null device") {plot.new()  } # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.expAll)
	mtext(grouping.expAll, side=1, cex=1.2)

	grouping.expAny <- "The Four experiment groups. Detected in ANY ONE of each."
	input.expAny <- list()
	for(thisExp in levels(expDesignLiverProteins$expRun))  {
		colName <- paste(thisExp,"Count",sep="")
		input.expAny[[thisExp]] <- proteinByLiverSample$spAccession[proteinByLiverSample[,colName] > 0 ]
	}
	if(names(dev.cur()) =="null device") {plot.new()}  # required because venn does not seem to be setting plot.new(), which is a requirement of mtext.
	venn(input.expAny)
	mtext(grouping.expAny, side=1, cex=1.2)
	
	par(mfrow=c(2,1))
	section.text <- "A histogram shows the scale of experiment-dependent protein detection (plus the number of proteins listed but with no data)"
	textplot(strwrap(section.text , width =textWidth ))
	hist(baseProteinByLiverSample$baseDetectCountAll, xlab="Numb. samples detected", main="Histogram of protein detection across samples", sub="Note zero class and multiples of 7", breaks=seq(-1,30))


	dev.off()
}


## there are 4 proteins in ALL adult samples and not in any other samples



## there are 14 protein in ALL 3 HepG2 samples and  not in any other samples.
#spec.HepG2 <- setdiff(setdiff(input.typeAll[[3]],input.typeAll[[1]]),input.typeAll[[2]])

#proteinByLiverSample[match(spec.HepG2, proteinByLiverSample$spAccession),]

