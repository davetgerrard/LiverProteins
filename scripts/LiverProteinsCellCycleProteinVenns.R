#### cell cycle GOs, which have much higher rank than 

go.list <- c("GO:0006977","GO:0051437","GO:0051436","GO:0000216")

set.list <- list()
for(thisGo in go.list) {
	set.list[[thisGo]] <- listProtsInGoFromList(thisGo,"BP",goDataCollection.all.binary)
}

venn(set.list)
mtext("All annotated proteins", side=1, cex=1.2)

set.list.detected <- list()
for(thisGo in go.list) {
	set.list.detected[[thisGo]] <- intersect(set.list[[thisGo]], detectedProts.all)
}


venn(set.list.detected)
mtext("Detected proteins", side=1, cex=1.2)


pdf("cellCycleGoTermProteins.pdf")
venn(set.list)
mtext("All annotated proteins", side=1, cex=1.2)
venn(set.list.detected)
mtext("Detected proteins", side=1, cex=1.2)
dev.off()