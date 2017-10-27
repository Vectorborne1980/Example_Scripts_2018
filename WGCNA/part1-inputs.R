################################################################################################
### WGCNA PART 1
# This script will create a datExpr data frame based on expression data from an eset.
# For WGCNA the columns of datExpr will be probes and the rows will be arrays. The script
# will also create a datTraits frame where columns are clinical traits and rows are arrays.
# The row names of datExpr must match the row names of datTraits.
# CURRENTLY, as of February 2017, AMPEL is using ONLY GCRMA normalized esets annotated with
# manufacturer CDFs (i.e., Affy, Illumina, etc.), and NOT Brainarray.
################################################################################################

### LOAD AND CONFIGURE WGCNA LIBRARY
library(Biobase)
library(limma)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

active.base <- "GSE10325 CD19 Standardized_CTL_SLEDAI_NA"
title.base <- "GSE10325 CD19 Standardized_CTL_SLEDAI_NA"
file.base <- "GSE10325_CD19_Standardized_CTL_SLEDAI_NA"
analysisDir <- "/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/Standardized_CTL_SLEDAI_NA/"
esetDir <- "/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/Standardized_CTL_SLEDAI_NA/"
globalDir <- "/mnt/disks/rserver-hdd/projects/"
projectDir <- "/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/"
scriptDir <- "/mnt/disks/rserver-hdd/projects/master_scripts/"
setwd(analysisDir)
save(globalDir,projectDir,scriptDir,esetDir,analysisDir,file.base,title.base,active.base,
     file=paste0(file.base,"_analysis_parameters.RData"))

###LOAD AFFY GCRMA NORMALIZED EXPRESSION DATA
setwd(analysisDir)
load("./eset_affy_gcrma_ampel.RData")
eset <- eset_affy_gcrma_ampel #Rename loaded eset
rm(eset_affy_gcrma_ampel) #Remove original eset

###Filter out any cell type or cohort that is not relevant to the present comparisons
table(eset$cell_type) #Check which cell types are there
index <- which(eset$cell_type == "CD19 B")
eset <- eset[,index]
eset$cell_type <- droplevels(eset$cell_type) #Only keep the data that have cell type associations
eset$cohort <- as.character(eset$cohort)
table(eset$cohort) #Check that SLE and CTL (controls) are present

###Load patient data
pdata <- read.table("./GSE10325_CD19-Bcells_Active_Patient-Data (revised_NSG_23Oct2017).txt", header=TRUE, sep="\t")
rownames(pdata) <- pdata[,1] #Make the first column be the row names
pdata <- pdata[,-1] #Remove the first column
sampleNames(eset) #Look at the sample names
rownames(pdata) #Look at the patient data row names
identical(rownames(pdata),sampleNames(eset)) #If FALSE, rename the rows in pData to match sample names in eset

###THIS SECION WILL BE UNIQUE FOR EACH DATASET!!!
table(pdata$cell_type) #Get rid of everything but B cells
index <- which(pdata$cell_type == "CD19 B")
pdata <- pdata[index,]
table(pdata$SLEDAI) #Only retain active SLE (SLEDAI >= 6)
index <- which(pdata$SLEDAI >= 6)
actives <- pdata[index,]
table(pdata$cohort) #Keep the healthy controls
index <- which(pdata$cohort == "healthy")
controls <- pdata[index,]
pdata <- rbind(actives,controls) #Put the controls and active SLE together
pdata$cohort <- gsub("healthy","CTL",pdata$cohort) #Rename the cohort variable to a standard type
pdata$SLEDAI[which(pdata$cohort=="CTL")] <- gsub(0,NA,pdata$SLEDAI[which(pdata$cohort=="CTL")])
rownames(pdata) <- paste(pdata$cohort, rownames(pdata),pdata$SLEDAI,sep="." ) #Add cohort to rownames
pdata <- pdata[order(rownames(pdata)),]
eset$SLEDAI[which(eset$cohort=="CTL")] <- gsub(0,NA,eset$SLEDAI[which(eset$cohort=="CTL")])
sampleNames(eset) <- paste(eset$cohort,sampleNames(eset),eset$SLEDAI,sep=".")
eset <- eset[,order(sampleNames(eset))]
sampleNames(eset) #Look at the sample names
rownames(pdata) #Look at the patient data row names
index <- which(sampleNames(eset)%in%rownames(pdata)) #Remove the inactives from the eset
eset <- eset[,index]
index <- which(rownames(pdata) != "CTL.GSM260892.NA") #One outlier was removed for DE analysis in the eset
pdata <- pdata[index,]
index <- which(sampleNames(eset)!="CTL.GSM260892.NA")
eset <- eset[,index]
identical(rownames(pdata),sampleNames(eset)) #Double check
rm(actives,controls,index)

pData(eset) <- pdata #Substitute patient data for the phenoData in the eset

###Create SYMBOL column required by some a4 functions
fData(eset)[,"SYMBOL"] <- as.character(fData(eset)[,"geneSymbol"]) #Take information from the featureData (fData) in the eset
fData(eset)[,"ENTREZ"] <- as.character(fData(eset)[,"geneEntrezID"]) #Take information from the featureData (fData) in the eset

###REMOVE AFFY PROBES
index <- which( substr(rownames(eset),1,4)!="AFFX" )
eset <- eset[index, ]
rm(index)

#######################################################################################
# Plot PCA

setwd(analysisDir)
library(affycoretools)
options(rgl.useNULL=TRUE)
library(rgl)

eset$cohort <- factor(eset$cohort)
pdf("PCA_Pre-filter.pdf",width = 16, height = 9)
plotPCA(eset, groups=as.numeric(eset$cohort), groupnames=levels(eset$cohort), main="Pre-filter",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), addtext=sampleNames(eset))
dev.off()
# plotPCA(eset, groups=as.numeric(eset$cohort), groupnames=levels(eset$cohort),
#         pcs=c(1,2,3), plot3d=TRUE, col=c("red","blue"), pch=c(16,17,17))

#######################################################################################
# Draw a quick dendrogram to look for outliers in the unfiltered set
library(flashClust)

sampleTree <- flashClust(dist(t(exprs(eset))), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
pdf("Dendrogram_Pre-filter.pdf",width = 16, height = 9)
par(cex = 0.9);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Pre-filter", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1.5)
dev.off()
### STOP: IF OUTLIERS ARE DETECTED IN THE DENDROGRAM, USE THIS
# Find clusters cut by the line and plot a line to show the cut
# abline(h = 115, col = "red"); #Determine cluster under the line
# clust <- cutreeStatic(sampleTree, cutHeight = 115, minSize = 10)
# table(clust)
# keepSamples <- (clust==1) #Cluster 1 contains the samples we want to keep.
# eset <- eset[,keepSamples] #SLE-GSM260919 was removed as an outlier
# nrow(eset)

# ### ALTERNATIVELY: IF OUTLIERS ARE DETECTED IN THE DENDROGRAM, USE THIS
# sampleNames(eset)
# eset <- eset[,sampleNames(eset)!="SLE-GSM260919"] #Remove the outlier sample(s) by name

#######################################################################################
###Filter out any unannotated probes (Currently NOT performing this)
# index <- which(is.na(fData(eset)[,"ENTREZ"])|fData(eset)[,"ENTREZ"]=="")
# eset <- eset[-index,]
# nrow(eset)
# 
# #######################################################################################
# # PLOT HISTOGRAM OF INTENSITIES AND FILTER LOW PASS PROBES IF NO FILTER HAS YET BEEN APPLIED
# 
pdf("Low_intensity_probe_filtration.pdf",width = 16,height = 9)
par(mfrow=c(1,2)) #Setting graphics parameters A vector of length 2, where the first argument specifies the number of rows and the second the number of columns of plots

nrow <- format(as.numeric(nrow(eset)),big.mark=",",scientific=F) #Make an object out of the number of rows that originally are in the eset

hist(rowMeans(exprs(eset)), breaks=100, xlab="Average log2 Probe Intensity",
     main="Active SLE All Intensities",
     sub=paste(nrow, "total probes"), xlim=c(2,12), col="yellow")

index <- which((rowMeans(exprs(eset))>2.301)|(rowMeans(exprs(eset))<2.199)) #Get the probes' row numbers whose expression averages are NOT between 2.2 and 2.3
nrow(eset) #How many probes are in the eset
length(index) #How many probes did NOT have average expression between 2.2 and 2.3
nrow(eset)-length(index) #How many probes do have average expression between 2.2 and 2.3

# Basically, zoom in to the data
hist(rowMeans(exprs(eset)), breaks=4000, xlab="Average log2 Probe Intensity",
     main="Intensities 2.2-2.3",
     sub=paste((nrow(eset)-length(index))," with 2.2 < mean < 2.3, out of ", nrow), xlim=c(2.2,2.3), col="yellow") #Write in a calculation to find how many probes are represented in this range

abline(v=2.22, lty=2, lwd=2, col="red")
dev.off()
# APPLY THE LOW PASS FILTER
index <- which(rowMeans(exprs(eset))<2.22) #Get row numbers for probes with average expression less than 2.27
eset <- eset[-index,] #Remove the probes with average expression less than 2.27
nrow(eset) #Check the number of probes remaining

#######################################################################################
###Create the datExpr0 frame for WGCNA
datExpr0 <- exprs(eset) #Take out the expression data ONLY from the eset dataframe
nrow(datExpr0) #Make sure the number of rows/probes matches the filtered eset

# IMPORTANT! TRANSPOSE THE ORIGINAL EXPRESSION MATRIX
# WGCNA requires columns be probes and rows be arrays
datExpr0 <- t(datExpr0) 

### CHECK FOR MISSING EXCESSIVE VALUES AND IDENTIFY OUTLIER MICROARRAY SAMPLES

# Check for gene entries with NAs and very low counts
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# Display and remove these bad entries
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    badGenes <- colnames(datExpr0)[!gsg$goodGenes];
  printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
} else {
  badGenes<-NULL
}

# We must remove probe columns containing only one value. Failure to do so breaks downstream TOM creation,
# and it must be performed now before module colors are assigned.
datExpr1 <- datExpr0[,apply(datExpr0,2,function(x) any(c(FALSE,x[-length(x)]!=x[-1])))]
ncol(datExpr0) - ncol(datExpr1)
datExpr0 <- datExpr1
rm(datExpr1,gsg)

# ##Remove bad probes eliminated by good genes function
datExpr1 <- t(datExpr0)
eset.wgcna <- eset[rownames(eset)%in%rownames(datExpr1),]
nrow(eset.wgcna)==nrow(datExpr1)
nrow(eset.wgcna) #Check the number of probes remaining
setwd(analysisDir)
save(eset.wgcna,file=paste0(file.base,"_eset_WGCNA.RData"))

#######################################################################################
# Plot PCA RECHECK
setwd(analysisDir)
library(affycoretools)
options(rgl.useNULL=TRUE)
library(rgl)

eset.wgcna$cohort <- factor(eset.wgcna$cohort)
pdf("PCA_Post-filter.pdf",width = 16, height = 9)
plotPCA(eset.wgcna, groups=as.numeric(eset.wgcna$cohort), groupnames=levels(eset.wgcna$cohort), main="Post-Filter",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), addtext=sampleNames(eset.wgcna))
dev.off()
# plotPCA(eset.wgcna, groups=as.numeric(eset.wgcna$cohort), groupnames=levels(eset.wgcna$cohort),
#         pcs=c(1,2,3), plot3d=TRUE, col=c("red","blue"), pch=c(16,17,17))

#######################################################################################
# Draw a quick dendrogram to look for outliers in the filtered set
library(flashClust)
sampleTree <- flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
pdf("Dendrogram_Post-filter.pdf",width = 16, height = 9)
par(cex = 0.9);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Post-Filter", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1.5)
dev.off()

datExpr <- datExpr0
ncol(datExpr)
nGenes <- ncol(datExpr) # number of columns (probes)
nSamples <- nrow(datExpr) # number of rows (arrays)
rm(sampleTree, datExpr0, datExpr1, index)

###########################################################################################
## LOAD CLINICAL TRAIT DATA (Applies to both filtered and unfiltered expression sets)

# Create a patients data frame analogous to expression data that will hold binary designations
# of trait membership, or a numerical measurement of a continuous clinical variable
patients.data <- pData(eset.wgcna)
rownames(patients.data) #Are they the same as the datExpr object, minus any outliers?
patients.data[is.na(patients.data)] <- ""
patients.data <- as.data.frame(patients.data)

colnames(patients.data)
datTraits <- patients.data[,c("cohort","SLEDAI","dsDNA","C3_C4","ESR","Prot","patient_age")]
rownames(datTraits) <- colnames(exprs(eset))

# Create numerical factors of binary designations

# cohort
table(patients.data$cohort)
datTraits$cohort <- gsub("CTL",0,datTraits$cohort)
datTraits$cohort <- gsub("SLE",1,datTraits$cohort)
datTraits$cohort <- as.numeric(datTraits$cohort)

# SLEDAI
table(patients.data$SLEDAI)
datTraits$SLEDAI <- as.numeric(datTraits$SLEDAI)

#dsDNA
table(patients.data$dsDNA)
datTraits$dsDNA <- gsub("normal",0,datTraits$dsDNA)
datTraits$dsDNA <- gsub("increased",1,datTraits$dsDNA)
datTraits$dsDNA <- as.numeric(datTraits$dsDNA)

#C3_C4
table(patients.data$C3_C4)
datTraits$C3_C4 <- gsub("low",0,datTraits$C3_C4)
datTraits$C3_C4 <- gsub("normal",1,datTraits$C3_C4)
datTraits$C3_C4 <- as.numeric(datTraits$C3_C4)

# ESR
table(patients.data$ESR)
datTraits$ESR <- gsub("nd",0,datTraits$ESR)
datTraits$ESR <- gsub(">130",130,datTraits$ESR)
datTraits$ESR <- as.numeric(datTraits$ESR)

# Prot
table(patients.data$Prot)
datTraits$Prot <- gsub("no",0,datTraits$Prot)
datTraits$Prot <- gsub("yes",1,datTraits$Prot)
datTraits$Prot <- as.numeric(datTraits$Prot)

# Age
table(patients.data$patient_age)
datTraits$patient_age <- as.numeric(datTraits$patient_age)

#####Clean up the dataset and display final dendrogram with corresponding clinical factors
collectGarbage();

# Re-cluster samples after outlier removal and display dendrogram
sampleTree <- flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath
pdf("Dendrogram_Traits.pdf",width = 16, height = 9)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Post-Filter" )
dev.off()
### SAVE OBJECTS 
save(datExpr, datTraits, eset.wgcna, badGenes, file="part1-final-objects.RData")
