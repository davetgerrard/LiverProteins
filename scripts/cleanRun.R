
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################


setwd("C:/Users/dave/LiverProteins/testDir")  # change this to wherever you want the output to come out.

#set this to TRUE if you want all the graphs to be output (mostly pdfs and tiffs)
output <- FALSE
#output <- TRUE


# change these absolute paths. Relative paths could be used if they were relative to the working directory above.
source("C:/Users/dave/LiverProteins/scripts/loadPCA_noOutput.R")

source("C:/Users/dave/LiverProteins/scripts/topGO_LiverProteins.R")

# this version of DAVID_funcClustering.R loads a pre-computed matrix
source("C:/Users/dave/LiverProteins/scripts/DAVID_funcClustering.R")

source("C:/Users/dave/LiverProteins/scripts/topGoByFuncCluster.R")






