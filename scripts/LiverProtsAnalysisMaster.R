
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#	2011				
#					
#########################################


test_headDir <- "C:/Users/dave/LiverProteins/testDir"
todayDir <- paste(test_headDir, Sys.Date(), sep="/")
if(!file.exists(todayDir)) {dir.create(todayDir) }
setwd(todayDir)  
#setwd("C:/Users/dave/LiverProteins/testDir")  # change this to wherever you want the output to come out.

#set this to TRUE if you want all the graphs to be output (mostly pdfs and tiffs)
output <- FALSE
#output <- TRUE


# change these absolute paths. Relative paths could be used if they were relative to the working directory above.
source("C:/Users/dave/LiverProteins/scripts/loadPCA_noOutput.R")


# plot venn diagrams of protein inclusion in sample groups.
source("C:/Users/dave/LiverProteins/scripts/LiverProteinsVennDiagrams.R")

# calculate functional clusters
# this version of DAVID_funcClustering.R loads a pre-computed matrix
source("C:/Users/dave/LiverProteins/scripts/DAVID_funcClustering.R")
# load utilities for performing functional clustering of GO result tables.
source("C:/Users/dave/LiverProteins/scripts/funcClusterGoTable.R")
# load utilities for running multiple GO analyses.
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProteinsUtilities.R")


# run detected/undetected GO analysis on all proteins
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProtsDetected.all.R")

# run detected/undetected GO analysis on ubiq proteins only
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProtsDetected.ubiq.R")

# run score based GO tests on PCs
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProtsPcScores.R")

# draw Vioplots of PC scores for functionally grouped significant GO terms
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProtsVioplotOnPcScores.R")

#source("C:/Users/dave/LiverProteins/scripts/topGO_LiverProteins.R")
#source("C:/Users/dave/LiverProteins/scripts/topGoByFuncCluster.R")






