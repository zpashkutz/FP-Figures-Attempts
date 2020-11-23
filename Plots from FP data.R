# Try to make a Beta Diversity plot using NMDS with the FP data by following the Alex Chase Tutorial

setwd("E:/School/ICU/Bio/199/Family Project/Figs Practice/Attempting FP figs from ABC")

library(vegan)
library(ggplot2)
library(EcolUtils)
library(biomformat)
library(stringr)
library(plyr)
library(gridExtra)

#import metadata
metadataFP <- read.delim("Microbiome_mapping_14Oct2019.txt")
str(metadataFP)

#get taxonomy - so functionally this will be equivalent to OTU_taxalevel from the Alex Chase tutorial 
OTU.taxonomyFP <- read.delim("taxonomy.txt", row.names = 1)

#here we can check to see how many unique phyla we can get - it should be 12
length(as.vector(unique(OTU.taxonomyFP$Phylum)))

#From here we can also look at the percentage of reads that are in each taxa level 
totalassignedtaxaFP <- length(which(is.na(OTU.taxonomyFP$Genus))) + 
  length(which(OTU.taxonomyFP$Genus==""))
  length(which(OTU.taxonomyFP$Kingdom == "Unassigned"))
(nrow(OTU.taxonomyFP)-totalassignedtaxaFP) / nrow(OTU.taxonomyFP) * 100
table(OTU.taxonomyFP$Genus)

#import OTU table and specify that the first column contains the row names 
OTU.tableFP <- read.delim("OTU_table.txt", row.names = 1)

#merge OTU taxonomy with OTU taxa level to filter? pls werk
OTU.table.with.taxonomyFP <- as.data.frame(merge(OTU.taxonomyFP, 
                                               OTU.tableFP, by.x = "row.names", by.y = "row.names"))
# hooray it works you goddamn monkey. Was that so hard, no? then why did it take you 3 hours. 

#Filtering Taxa, probably omit this step since now we're down to 87 samples and the original study had 95 where only one was omitted
#you should try this again but later
OTU.table.filteredtaxaFP <- OTU.table.with.taxonomyFP[!grepl("Unassigned|chloroplast", 
                                                         OTU.table.with.taxonomyFP$Kingdom),]
#subset only OTU table with the filtered taxa
FPfilteredtaxa <- OTU.table.filteredtaxaFP$Row.names

FPfilterOTU <- subset(OTU.tableFP, rownames(OTU.tableFP) %in% FPfilteredtaxa)

FPfinalOTU <- FPfilterOTU[,!names(FPfilterOTU) %in% c("Mock")]

# get quartile ranges for rarefaction
FPtransOTU <- rowSums(t(FPfinalOTU)) 
Q10 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.10)
Q15 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.15)

#look at ranges
barplot(sort(FPtransOTU), ylim = c(0, max(FPtransOTU)), 
        xlim = c(0, NROW(FPtransOTU)), col = "Blue", ylab = "Read Depth", xlab = "Sample") 
abline(h = c(Q10, Q15), col = c("red", "pink"))
plot.new()

#lines here did not plateau meaning that we did not sequence the entirety of the community 
#because of this I'm going to keep everything in the 85th percentile, so to the left of the pink line 
rarecurve(t(FPfinalOTU), step = 100, cex = 0.5)
abline(v = c(Q10, Q15), col = c("red", "pink"))

#deciding to keep everything to the left of the pink line in the rarefied OTU table 
FPrared.OTU <- as.data.frame((rrarefy.perm(t(FPfinalOTU), sample = Q15, n = 100, round.out = T)))
#only keep what meets the rarefaction cutoff
FPrared.OTU <- as.data.frame(FPrared.OTU[rowSums(FPrared.OTU) >= Q15 - (Q15 * 0.15), colSums(FPrared.OTU) >= 1]) 

#Lets look at this OTU to see if we messed it up 
View(FPrared.OTU)

#Beta-diversity NMDS plots 
# make that bray-curtis dissimilarity matrix 

NMDS.fam <- metaMDS(FPrared.OTU, distance = "bray", k = 2, trymax = 500)

#extract the two coordinates of the matrix 

coordinates.fam <- data.frame(NMDS.fam$points[,1:2]) 

#look at the plot 
plot(x = coordinates.fam$MDS1, y = coordinates.fam$MDS2) 

#add X in front of sample in metadataFP
metadataFP$X.NAME <- interaction( "X", metadataFP$X.NAME, sep = "")
rownames(metadataFP) <- metadataFP[,1]
metadataFP[,1] <- NULL 
View(metadataFP)

#merge with metadataFP
nmds.fam.metadata <- merge(coordinates.fam, metadataFP, by = 0)

Factor.F <- as.factor(nmds.fam.metadata$Family)

# Plot
ggplot(data = nmds.fam.metadata) +
  aes(x = MDS1, y = MDS2, color = Family)+ #Creates and colors legend to match, modify after $ here.
  geom_point(size = 3) +
  labs(col = "Family") + #Renames legend, modifiable. 
  theme_bw()

# for statistical tests, we need to get the data tidied up a bit
# so let's subset the metadata and keeps only the samples that passed filtering and rarefaction
filtersamplesFP <- rownames(FPrared.OTU)
filtermetaFP <- subset(metadataFP, rownames(metadataFP) %in% filtersamplesFP)

adonis(data = filtermetaFP, formula = FPrared.OTU ~ Family/ Individual/
         Individual + Family,
       permutations = 999, method = "bray")


ggplot(data = nmds.fam.metadata) +
  aes(x = MDS1, y = MDS2, color = Factor.F) + #Creates and colors legend to match, modify after $.
  geom_point() +
  labs(col = "Family") + #Renames legend, modifiable within quotes.
  ggtitle("NMDS of Families' Oral Microbiomes", subtitle = bquote(~R^2~ '= 0.199, p = 0.001')) +#Adds tittle and subtittle. Can modify p and r-squared values + title.
  theme_classic(base_size = 14, base_line_size = .5)

