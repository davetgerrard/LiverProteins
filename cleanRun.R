
setwd("C:/Users/dave/LiverProteins/testDir")

#output <- TRUE

source("C:/Users/dave/LiverProteins/loadPCA_noOutput.R")

source("C:/Users/dave/LiverProteins/topGO_LiverProteins.R")

# this version of DAVID_funcClustering.R loads a pre-computed matrix
source("C:/Users/dave/LiverProteins/DAVID_funcClustering.R")

source("C:/Users/dave/LiverProteins/topGoByFuncCluster.R")





###### NICER FORMAT


# parameters

runName <- paste(Sys.Date(),"_output",sep="")
setwd("C:/Users/dave/LiverProteins/")
workDir <- getwd()
dataDir <- paste(workDir,"/data/",sep="")
outputDir <- paste(workDir,runName,sep="/")
fileOutput <- TRUE
screenOutput <- TRUE

# load the data


# clean the data


# select ubiquitous proteins


# perform PCA


# perform GO on all detected proteins


# perform GO on PC scores


# filter by functional groups