
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################

##INFO: Utility functions to run multiple GO analyses over goData objects and store results in tables.

########## FUNCTIONS


listProtsInGoFromList <- function(goTerm,ontology,goDataList)  {	#INFO: list proteins for a GO term in a precalculated list of goData objects, accounting for ontology
	genesInTerm(goDataList[[ontology]],goTerm)[[1]]
}


topDiffGenes <- function(allScore) {	#INFO: the topGOdata objects require a method for gene/protein selection. For enrichment analysis, we keep all proteins.
	return(allScore )
}


resultSummaryAsText <-function (x)    #INFO: this was a hack to get results out of the GOdata object and used for text output. May be out-of-date.
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


runDetectGoTests <- function(GOdataBase )  {	#INFO: runs a suite of GO detected/undetected tests over a topGOdata object and returns a result table

	thisGOgraph <- GOdataBase@ontology
	#GOdataBase <- goDataCollection.ubiq.binary[[thisGOgraph]]
	# not sure if updateGenes() required here. 
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
	#goGroup <- as.character(thisGoResults$goTerm)
	#thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdataBase,x)[[1]])))

	#summaryDetectResults.ubiq <- rbind(summaryDetectResults.ubiq,thisGoResults)
	#rm(manual.Fisher, manual.elimFisher, manual.termStats,thisGoResults)
	#rm(resultElimFisherAdj, resultElimFisher,  resultFisher,allResFT,headText)
	return(thisGoResults)

}


runScoreGoTests <- function(GOdata,geneList,geneSelectionFun=topDiffGenes)  { #INFO: runs a suite of GO score tests over a topGOdata object and returns a result table
		thisGOgraph <- GOdata@ontology
	
		####INFO: define test statistics and apply over all GO terms 
		#INFO: resultWilcox = 2-sided wilcox test. Good for outliers in one direction
		test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
		resultWilcox <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcox = elim version of resultWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcox <- getSigGroups(GOdata, test.stat)

		#resultElimWilcoxAdj <- resultElimWilcox 
		#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue
		#INFO: resultKS = Kolmogorov smirnov test: may be good for general distribution changes.
		test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
		resultKS <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimKS = elim version of resultKS
		test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
		resultElimKS <- getSigGroups(GOdata, test.stat)

		#resultElimKSAdj <- resultElimKS
		#score(resultElimKSAdj ) <- qvalue(score(resultElimKS))$qvalue

		#INFO: change scores to absolute values for Greater than test 
		geneList <- abs(ubi.pca.5.scores[,compHead])
		names(geneList) <- ubi.pca.5.scores$spAccession
		GOdata <- updateGenes(GOdata,geneList,topDiffGenes)
		#INFO: resultWilcoxAbs = test for outliers in both directions simultaneously. Less power than wilcox if true difference is unidireectional.
		test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
		resultWilcoxAbs <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcoxAbs = elim version of resultElimWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)

		#resultElimWilcoxAbsAdj <- resultElimWilcoxAbs 
		#score(resultElimWilcoxAbsAdj ) <- qvalue(score(resultElimWilcoxAbs ))$qvalue


		#INFO: collect and bind all the results together. 
		# This could be done with a function available in topGO but I wanted extra columns and control.
		Wilcox <- data.frame(Wilcox=score(resultWilcox ),goTerm=names(score(resultWilcox)))
		elimWilcox <- data.frame(elimWilcox=score(resultElimWilcox),goTerm=names(score(resultElimWilcox)))
		KS <- data.frame(KS=score(resultKS),goTerm=names(score(resultKS)))
		elimKS <- data.frame(elimKS=score(resultElimKS),goTerm=names(score(resultElimKS)))
		absWilcox <- data.frame(absWilcox=score(resultWilcoxAbs),goTerm=names(score(resultWilcoxAbs)))
		elimAbsWilcox <- data.frame(elimAbsWilcox=score(resultElimWilcoxAbs),goTerm=names(score(resultElimWilcoxAbs)))

		thisGoResults <- NULL
		thisGoResults <- merge(Wilcox,elimWilcox,by="goTerm")
		thisGoResults <- merge(thisGoResults ,KS ,by="goTerm")
		thisGoResults <- merge(thisGoResults ,elimKS ,by="goTerm")
		thisGoResults <- merge(thisGoResults ,absWilcox ,by="goTerm")
		thisGoResults <- merge(thisGoResults ,elimAbsWilcox ,by="goTerm")
		thisGoResults$ontology <- thisGOgraph
		thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(thisGoResults$goTerm, ontology(GOdata)))

		goGroup <- as.character(thisGoResults$goTerm)
		thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdata,x)[[1]])))
		
		#INFO: bind the results from this GO ontology as rows to the table for the current PC in 'summaryResultsList'
		#summaryResultsList[[i]] <- rbind(summaryResultsList[[i]],thisGoResults)
		return(thisGoResults)
}


GOWilcoxTestGreater <- function (object,alternativeType="greater")	#INFO: Wilcox 1-sided (greater) test for use in topGO 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}



GOWilcoxTest2Sided <- function (object,alternativeType="two.sided")	#INFO: Wilcox 2-sided test for use in topGO  
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}


