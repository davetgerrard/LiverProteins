## sharing of proteins between terms


# good example to illustrate that Go terms from different ontologies can be represented by the same set of proteins



goProteins.1 <- na.omit(listProtsInGo("GO:0005840","CC"))	# ribosome
goProteins.2 <- na.omit(listProtsInGo("GO:0003735","MF"))	# structural constituent

length(goProteins.1) ; length(goProteins.2); length(intersect(goProteins.1,goProteins.2))
# all 159 'struct const ribosome' proteins are included within the 200 'ribosome' proteins

goProteins.1 <- na.omit(listProtsInGo("GO:0055114","BP"))	# oxidation reduction
goProteins.2 <- na.omit(listProtsInGo("GO:0016491","MF"))	# oxidoreductase activity
 length(goProteins.1) ; length(goProteins.2); length(intersect(goProteins.1,goProteins.2))
#90% overlap.

goProteins.1 <- na.omit(listProtsInGo("GO:0003723","MF"))	# RNA binding
goProteins.2 <- na.omit(listProtsInGo("GO:0030529","CC"))	# ribonucleoprotein complex
length(goProteins.1) ; length(goProteins.2); length(intersect(goProteins.1,goProteins.2))
#763, 520, only 269 overlaps. 

goProteins.1 <- na.omit(listProtsInGo("GO:0030168","BP"))	# platelet activation
goProteins.2 <- na.omit(listProtsInGo("GO:0002576","BP"))	# platelet degranulation
length(goProteins.1) ; length(goProteins.2); length(intersect(goProteins.1,goProteins.2))
#236, 80, 80. Degranulation is entirely contained withing activation.

goProteins.1 <- na.omit(listProtsInGo("GO:0019825","MF"))	# oxygen binding
goProteins.2 <- na.omit(listProtsInGo("GO:0020037","MF"))	# heme binding
length(goProteins.1) ; length(goProteins.2); length(intersect(goProteins.1,goProteins.2))
# look at PC scores for heme binding
proteinByLiverSample.ubi.pca[match(goProteins.2,proteinByLiverSample.ubi.pca$spAccession),]


goProteins <- na.omit(listProtsInGo("GO:0005759","CC"))	# mito matrix
goProteins <- na.omit(listProtsInGo("GO:0005743","CC")) 	# mito inner membrane
goProteins <- na.omit(listProtsInGo("GO:0005739","CC"))	# mitochondrion
goProteins.valid <- intersect(validProts,goProteins)

goProteins.1 <- na.omit(listProtsInGo("GO:0055114","BP"))	# 
goProteins.2 <- na.omit(listProtsInGo("GO:0016491","MF"))	#