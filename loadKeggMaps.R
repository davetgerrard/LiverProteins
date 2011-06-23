library(KEGG.db) # required for names.

rm(xx,protKeggMap,keggProtPathMap,pathCounts,proteinCounts,pathList)

xx <- as.list(KEGGPATHID2NAME)


protKeggMap <- read.delim("C:/Users/dave/LiverProteins/hsa_uniprot.list",header=F)
names(protKeggMap) <- c("keggProtId","uniprotId")
protKeggMap$uniprotId <- sub("up:","",protKeggMap$uniprotId,fixed=T)

keggProtPathMap <- read.delim("C:/Users/dave/LiverProteins/hsa_pathway.list",header=F)
names(keggProtPathMap) <- c("keggProtId","keggPathId")

proteinCounts <- as.data.frame(table(keggProtPathMap$keggProtId))
names(proteinCounts) <- c("keggProtId","keggPathsWithProtCount")
protKeggMap <- merge(protKeggMap,proteinCounts,by="keggProtId")

pathCounts <- as.data.frame(table(keggProtPathMap$keggPathId))
names(pathCounts) <- c("keggPathId","keggProtsInPathCount")
pathList <- sub('path:hsa','',fixed=T,as.character(pathCounts $keggPathId))
pathCounts$keggPathName <- as.character(lapply(pathList, FUN=function(x) xx[[x]]))

keggPathToUniprot <- merge(protKeggMap,keggProtPathMap,by="keggProtId")

rm(proteinCounts,pathList)



