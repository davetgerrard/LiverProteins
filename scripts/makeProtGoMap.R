
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

source("C:/Users/dave/LiverProteins/loadGoMaps.R")

prot_list <- unique(protGoMap$DB_Object_ID)
	
listGoPerProt <- function(dataFrame,thisProt) {
	dataFrame$GO_ID[which(dataFrame$DB_Object_ID == thisProt)]
}

prot2go <- lapply(prot_list,FUN = function (x) listGoPerProt(protGoMap,x))
names(prot2go) <- prot_list
go2prot <- inverseList(prot2go)
prot2go <- inverseList(go2prot)
####need to save this, saw something similar in a script:-
a <- lapply(prot2go, function(x) paste(x, collapse = ", "))
mat <- cbind(names(a), unlist(a, use.names = FALSE))
write.table(mat, file = "go2prot.map", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

### can be read in and made a mapping like this.
#library(topGO)
#prot2go <-  readMappings("go2prot.map")
