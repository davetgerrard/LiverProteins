allLoadings <- loadings(ubiquitous.pca)[,]
as.data.frame(allLoadings)
allLoadings <- as.data.frame(allLoadings)
allLoadings$age <- expDesignLiverProteins$sampleAge.weeks[match(row.names(allLoadings),expDesignLiverProteins$sample)]

plot(allLoadings[,1],allLoadings$age)
cor.test(allLoadings[,1],allLoadings$age)

corPvalues <- numeric()

for (i in 1:20) {
	corPvalues[i] <- cor.test(allLoadings[,i],allLoadings$age)$p.value
}

#perhaps #4?