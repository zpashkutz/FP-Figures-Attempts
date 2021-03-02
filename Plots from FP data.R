# Try to make a Beta Diversity plot using NMDS with the FP data by following the Alex Chase Tutorial

setwd("E:/School/ICU/Bio/199/Family project - Quarter 1/Family Project/Figs Practice/Attempting FP figs from ABC")

#libraries used 
library(vegan)
library(ggplot2)
library(EcolUtils)
library(biomformat)
library(stringr)
library(plyr)
library(gridExtra)
library(readxl) #maybe not this one

#import metadata
metadataFP <- read.delim("Microbiome_mapping_14Oct2019.txt", sep = "")
str(metadataFP)
#add X in front of sample in metadataFP
metadataFP$X.NAME <- interaction( "X", metadataFP$X.NAME, sep = "")
rownames(metadataFP) <- metadataFP[,1]
metadataFP[,1] <- NULL 


##md <- metadataFP
##names <- rownames(md)
##rownames(md) <- NULL
##newmd <- cbind(names, md)


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

#merge OTU taxonomy with OTU taxa level to filter?
OTU.table.with.taxonomyFP <- as.data.frame(merge(OTU.taxonomyFP, 
                                               OTU.tableFP, by.x = "row.names", by.y = "row.names"))

#OMITTING FILTERING BEFORE RAREFACTION undo "##" to enable filtering 
##################################################################
#Filtering Taxa - lets try this by omitting filtering and going straight to rarefaction

OTU.table.filteredtaxaFP <- OTU.table.with.taxonomyFP[!grepl("Unassigned", 
                                                         OTU.table.with.taxonomyFP$Kingdom),]
#subset only OTU table with the filtered taxa
FPfilteredtaxa <- OTU.table.filteredtaxaFP$Row.names
FPfilterOTU <- subset(OTU.tableFP, rownames(OTU.tableFP) %in% FPfilteredtaxa)

#For us the below line doesn't do anything but, if there were any mock community in qiita and qiime2 we should have filtered them out - to my understanding
#FPfinalOTU <- FPfilterOTU[,!names(FPfilterOTU) %in% c("Mock")]

# get quartile ranges for rarefaction
FPtransOTU <- rowSums(t(OTU.tableFP)) 
Q0 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.0)
Q10 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.10)
Q15 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.15)
Q25 <- quantile(FPtransOTU[order(FPtransOTU, decreasing = TRUE)], 0.25)
#look at ranges
barplot(sort(FPtransOTU), ylim = c(0, max(FPtransOTU)), 
        xlim = c(0, NROW(FPtransOTU)), col = "Blue", ylab = "Read Depth", xlab = "Sample") 
abline(h = c(Q0, Q10, Q15,Q25), col = c("green","red", "pink", "cyan"))
plot.new()

######################### Right here is where I think I messed up, the samples got cut in size from 95 to 89
#lines here did not plateau meaning that we did not sequence the entirety of the community 
#because of this I'm going to keep everything in the 90th percentile, so to the left of the red line 
rarecurve(t(OTU.tableFP), step = 100, cex = 0.5)
abline(v = c(Q0, Q10, Q15, Q25), col = c("green","red", "pink", "cyan"))

#deciding to keep everything to the left of the red line in the rarefied OTU table 
FPrared.OTU <- as.data.frame(rrarefy.perm(t(OTU.tableFP)), sample = Q0, n = 100, round.out = T)

raredreads <- as.data.frame(rowSums(FPrared.OTU))
colnames(raredreads) [1] <- "read count"
#only keep what meets the rarefaction cutoff -  so don't do this then?   This function throws out all samples 
#to the left of the chosen line
#FPrared.OTU <- as.data.frame(FPrared.OTU[rowSums(FPrared.OTU) >= Q0 - (Q0 * 0.1), colSums(FPrared.OTU) >= 1]) 

#Lets look at this OTU to see if we messed it up 
View(FPrared.OTU)


# in the following exercises: 
#######################################################################
#NO FILTERING AND NO RAREFATION USE "OTU.tableFP" 
#RUN WITHOUT ANY LINES DISABLED FOR FILTERING AND RAREFACTION, USE,"FPrared.OTU"
#RUN WITH ONLY FILTERING DISABLED FOR NO FILTERING BUT RAREFACTION, USE "FPrared.OTU"
#RUN WITH ONLY FILTERING ENALBED AND DISABLE RAREFACTION FOR ONLY FILTERING, USE "FPfinalOTU"
#######################################################################



#Beta-diversity NMDS plots 
# make that bray-curtis dissimilarity matrix 

NMDS.fam <- metaMDS(FPrared.OTU, distance = "bray", k = 2, trymax = 500, wascores = T)

#extract the two coordinates of the matrix 

coordinates.fam <- data.frame(NMDS.fam$points[,1:2]) 

#look at the plot 
plot(x = coordinates.fam$MDS1, y = coordinates.fam$MDS2) 

#merge with metadataFP
nmds.fam.metadata <- merge(coordinates.fam, metadataFP, by = 0, group=metadataFP$Family)
rownames(nmds.fam.metadata) <- nmds.fam.metadata[,1]
nmds.fam.metadata[,1] <- NULL 



Factor.F <- as.factor(nmds.fam.metadata$Family)

# Plot
ggplot(data = nmds.fam.metadata) +
  aes(x = MDS1, y = MDS2, color = Family)+ #Creates and colors legend to match, modify after $ here.
  geom_point(size = 3) +
  labs(col = "Family") + #Renames legend, modifiable. 
  theme_bw()

# for statistical tests, we need to get the data tidied up a bit
# so let's subset the metadata and keeps only the samples that passed filtering and rarefaction

# If you want to try this without rarefying first sub FPrared.OTU for OTU.tableFP
filtersamplesFP <- rownames(FPrared.OTU)
filtermetaFP <- subset(metadataFP, rownames(metadataFP) %in% filtersamplesFP)

adonis(data = filtermetaFP, formula = FPrared.OTU ~ Family/ Individual/
         Individual + Family,
       permutations = 999, method = "bray")

#plot the actual NMDS 
ggplot(data = nmds.fam.metadata) +
  aes(x = MDS1, y = MDS2, color = Factor.F) + #Creates and colors legend to match, modify after $.
  geom_point() + 
  labs(col = "Family") + #Renames legend, modifiable within quotes.
  ggtitle("NMDS of Families' Oral Microbiomes", subtitle = bquote(~R^2~ '= 0.057, p = 0.006')) + #Adds tittle and subtitle. Can modify p and r-squared values + title.
  theme_classic(base_size = 14, base_line_size = .5) +
  stat_ellipse(size=1)


###############
###############
###############

#plot the actual NMDS - This one should be exactly the same as the one above
#Here I'm going to try to use ordiellipse to make the ellipses - idk if it's gonna work tho

#Okay so this catastrophically failed, I'll keep this in here for now to revisit -  try different methods for now
# all the lines with "###" in front designate code for the failed ordiellipse attempt

###############
###############
###############

###ord <-ordiellipse(NMDS.fam, nmds.fam.metadata$Family, display = "sites", 
             ###   kind = "sd", conf = 0.95, label = F)

###df_ell <- data.frame()
###  for(g in levels(nmds.fam.metadata$Family)){
  ###  df_ell <- rbind(df_ell, cbind(as.data.frame(with(nmds.fam.metadata[nmds.fam.metadata$Family==g,],
   ###                                                 veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                  ###               ,Factor.F=g))
###}


###ggplot(data = nmds.fam.metadata) +
 ### aes(x = MDS1, y = MDS2, color = Factor.F) + #Creates and colors legend to match, modify after $.
###  geom_point(aes(color = Factor.F)) + geom_path(data = df_ell, aes(x=MDS1, y=MDS2, color=Factor.F), size =1, linetype=1) +
 ### labs(col = "Family") + #Renames legend, modifiable within quotes.
 ### ggtitle("NMDS of Families' Oral Microbiomes", subtitle = bquote(~R^2~ '= 0.199, p = 0.001')) + #Adds tittle and subtitle. Can modify p and r-squared values + title.
  ###theme_classic(base_size = 14, base_line_size = .5)



#Make some alpha diversity plots and see if they line up with what I can generate from MicrobiomeAnalyst 
#Originally had FPrared.OTU (2 instances below) in place of OTU.tableFP because rarefied values were used, note you must transpose the "OTU.tableFP to use" 
FPrichness <- as.data.frame(specnumber(FPrared.OTU))
colnames(FPrichness) <- c("speciesrich")
#Merge with metadata for plot
FPmerged.rich <- merge(FPrichness, metadataFP, by = 0)
rownames(FPmerged.rich) <- FPmerged.rich$Row.names
FPmerged.rich$Row.names <- NULL

#using Shannon diversity as a measure of alpha diversity 
FPshannon <- as.data.frame(diversity(FPrared.OTU, index = "shannon"))
colnames(FPshannon) <- c('alpha_shannon')

#get all data in one data.frame 
#make the sample identifiers "X##" the row names again instead of first column 

FPmerged.fam.alpha <- merge(FPmerged.rich, FPshannon, by = 0)
rownames(FPmerged.fam.alpha) <- FPmerged.fam.alpha$Row.names
FPmerged.fam.alpha[[1]] <- NULL
View(FPmerged.fam.alpha)

#Plot it Can do this with either species richness on the y axis or the shannon index it's the 
#same code just pick one and change the parameters to fit what you want to plot 
#I'm doing shannon index for ease of comparison with MA

Family.alpha <- as.factor(FPmerged.fam.alpha$Family)

FPp1 <- ggplot(data = FPmerged.rich)+
  aes(x = Family.alpha, y = FPmerged.fam.alpha$alpha_shannon,
      fill = Family.alpha) +
  geom_boxplot(outlier.shape = NA, lwd = 1) + 
  labs(title = "Alpha Diverisity Within Families", x = "Family", y = "Shannon Index", 
       fill = " Family") +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_jitter(width = .2) +
  theme(legend.position = "none")

FPp1

#well I guess i did it with the richness too 
FPp2 <- ggplot(data = FPmerged.rich)+
  aes(x = Family.alpha, y = FPmerged.fam.alpha$speciesrich,
      fill = Family.alpha) +
  geom_boxplot(outlier.shape = NA, lwd = 1) + 
  labs(title = "Species Richness Between Families", x = "Family", y = "Species Richness", 
       fill = " Family") +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_jitter(width = .2) +
  theme(legend.position = "none")

FPp2

