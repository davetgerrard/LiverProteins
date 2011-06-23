



allGo <- read.delim("C:/Users/dave/LiverProteins/gene_association.goa_human",comment.char="!",header=F,colClasses = "character")
names(allGo)<- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID","DB_Reference","Evidence","With","Aspect","DB_Object_Name","Synonym","DB_Object_Type","Taxon_ID","Date","Assigned_By","Annotation_Extension ","Gene_Product_Form_ID")
nrow(allGo)

allGo <- allGo[grep('NOT',allGo$Qualifier,invert=T),]
nrow(allGo)

protGoMap <- unique(subset(allGo,select=c(DB_Object_ID,GO_ID)))
goTermCounts <- as.data.frame(table(protGoMap$GO_ID))
names(goTermCounts) <- c("GO_ID","GoWithProteinCount")
protGoTermCounts <- as.data.frame(table(protGoMap$DB_Object_ID))
names(protGoTermCounts ) <- c("DB_Object_ID","GoTermsForThisProtein")

rm(allGo)




#####

#goIds <- 		# put some GO ids here.

#library("GO.db")
#head(as.data.frame(Term(goIds)))
