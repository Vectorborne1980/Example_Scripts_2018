# Install packages as needed
source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
# biocLite("genoset")
#biocLite("Biobase")
#biocLite("genefilter")
#biocLite("a4")
#biocLite("siggenes")
#biocLite("affycoretools")
#biocLite("xlsx")
#biocLite("beadarray")

# Load packages
library(Biobase)
library(genoset)
library(limma)
library(affycoretools)
library(rgl)
library(genefilter)
library(a4)
library(siggenes)
library(beadarray)


###########################################################
# READ IN EXPRESSION DATA
# The GEO submitters call this "raw" data, when in fact there's only illumina probe intensities and associated
# detection p value. The Illumina control files haven't been provided, which is needed to run beadarray or lumi packages.
# Instead, we'll resort to creating our own eset and normalizing using neqc.
# Acquire the MINiML formatted family files from GEO. Unzip and untar into a separate folder "GSE72535_family.xml"
# Delete or move the files "GSE72535_family.xml" and "GPL10558-tbl-1.txt". Merge the individual sample files into one compiled data table.
# Modify the non-normalized data frame manually to contain a ".Sample" affixed to the sample names
setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/WGCNA/WGCNA_Standardized_Re-analyses/Skin/GSE72535_New_Skin/ESET/")
file_list <- list.files("./GSE72535_family")
nrow(read.delim("./GSE72535_family/GSM1864221-tbl-1.txt",header = FALSE,stringsAsFactors = FALSE)) #Check the number of rows
dataset <- data.frame(matrix("", ncol = 1, nrow = 47322)) #Create a dataframe to which columns may be added
setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/WGCNA/WGCNA_Standardized_Re-analyses/Skin/GSE72535_New_Skin/ESET/GSE72535_family/")
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header = FALSE, sep="\t",stringsAsFactors = FALSE)
  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header = FALSE, sep="\t",stringsAsFactors = FALSE)
    dataset<-cbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/WGCNA/WGCNA_Standardized_Re-analyses/Skin/GSE72535_New_Skin/ESET/")
dataset <- dataset[,-1] #Get rid of the first empty column
rownames(dataset) <- dataset$V1 #Make the first probe column the rownames
c <- 1:ncol(dataset) #Make an object counting the number of columns
c%%3 #Give columns numeric designations by groups of three
dataset <- dataset[, !(c%%3==1)] #Drop all of the columns with probe names
file_list <- substr(file_list,1,nchar(file_list)-10) #Remove file appendix from file list names
file_list <- paste(file_list,"Sample",sep=".") #Append ".Sample" to the file list names
index <- (2*(1:17))-1 #Make an index of numbers representing every other column beginning with column 1
samples <- dataset[,index] #Collect just the non-normalized sample probe intensity values
colnames(samples) <- file_list #Make the column names of the non-normalized probe values that of file list sample names
file_list <- gsub("Sample","Detection",file_list) #Append "Detection" to the file list names
index <- 2*(1:17) #Make an index of numbers representing every other column beginning with column 2
detection <- dataset[,index] #Collect just the non-normalized sample detection p-values
colnames(detection) <- file_list #Make the column names of the detection p-values that of file list sample names
temp <- cbind(samples,detection) #Join the sample probe values and their detection p-values
n <- ncol(samples) #Get the number of columns
dataset <- temp[,kronecker(1:n, c(0, n), "+") ] #Insert the detection p-value columns behind each of their corresponding probe intensity value columns
dataset <- cbind(ID_Ref = rownames(dataset),dataset) #Add a separate column of the probe IDs onto the dataset
exprs.raw <- dataset #Rename dataset as exprs.raw
rm(c,n,file,file_list,temp,samples,detection,index,dataset)
write.table(exprs.raw,file = "GSE72535_family/GSE72535_family-compiled.txt")
nrow(exprs.raw) # 47,322 probes

#######################################################
# Open fData. According to the supplementary report accompanying the study RNA was hybridized to 
# Illumina HT12 V4 beadchips. The manufacturer's probe description file was downloaded from
# http://support.illumina.com/array/array_kits/humanht-12_v4_expression_beadchip_kit/downloads.html

# Figure out which annotation to use.
# R1 lists 47,231 probes and 887 control probes
# R2 lists 47,323 probes and 887 control probes
# All of the matching probe IDs in R1 and R2 have the same exact annotations listed, ergo we'll use R2
fData <- read.delim("GPL10558_annot.txt",skip=28,sep="\t", stringsAsFactors = FALSE)
fData <- fData[1:48107,]
index <- which(fData$ID%in%rownames(exprs.raw))
fData <- fData[index,]
rownames(fData) <- fData$ID
fData <- fData[order(rownames(fData)), ]
exprs.raw <- exprs.raw[ order(rownames(exprs.raw)), ]
identical(rownames(fData),rownames(exprs.raw))

# Append annotation columns to the exprs.raw frame and write out a new file for read.illmn
exprs.annot <- cbind(exprs.raw, fData)
write.table(exprs.annot, file="GSE72535_family-annotated.txt", sep="\t", quote=FALSE, row.names=FALSE)
colnames(exprs.annot)

#######################################################
# USE THE ILLUMINA FUNCTION IN THE LIMMA PACKAGE
# This will generate an EListRaw object

# New read with annotation (Be sure to leave ID_Ref as a separate column )
annotation.cols <- colnames(fData)
GSE72535.elist <- read.ilmn(files="GSE72535_family-annotated.txt",expr="Sample",other.columns="Detection",
                            probeid="ID_Ref", annotation=annotation.cols, verbose=TRUE)


###########################################################
# NORMALIZE USING NEQC

# Normalize using the neqc package
GSE72535.neqc <- normaliseIllumina(GSE72535.elist,method = "quantile",transform = neqc)
GSE72535.neqc <- neqc(GSE72535.elist)

# GSE72535.exprs <- GSE72535.elist$E
# rownames(GSE72535.exprs) <- rownames(exprs.raw)

GSE72535.neqc.exprs <- GSE72535.neqc$E
rownames(GSE72535.neqc.exprs) <- rownames(exprs.raw)

# dim(GSE72535.exprs)
dim(GSE72535.neqc.exprs)

###########################################################
# RENAME EXPRESSION COLUMNS TO GSM NUMBERS
pData <- read.delim("pdata_complete.txt", stringsAsFactors=FALSE)
index <- 2*(1:17)
colnames(GSE72535.neqc.exprs) <- colnames(exprs.raw[,index]) #Skip if using the family based data file
pData <- pData[ order(pData$sampleid), ]
GSE72535.neqc.exprs <- GSE72535.neqc.exprs[ , order(colnames(GSE72535.neqc.exprs))]
rownames(pData) <- paste(pData$cohort,sep=".",pData$sampleid)
rownames(pData) <- gsub(" ","",rownames(pData))
colnames(GSE72535.neqc.exprs) <- rownames(pData)
identical( rownames(pData),colnames(GSE72535.neqc.exprs) )

fData <- fData[ order(fData$ID), ]
GSE72535.neqc.exprs <- GSE72535.neqc.exprs[ order(rownames(GSE72535.neqc.exprs)), ]
identical(rownames(fData),rownames(GSE72535.neqc.exprs))

###########################################################
# CREATE ESET

# Biobase expects a matrix of expression values to generate an eset
exprs.eset <- as.matrix(GSE72535.neqc.exprs)
rownames(exprs.eset) <- rownames(GSE72535.neqc.exprs)

pData.meta <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData) )
pData.eset <- new("AnnotatedDataFrame", data = data.frame(pData), varMetadata = pData.meta)

fData.meta <- data.frame(labelDescription = colnames(fData), row.names = colnames(fData) )
fData.eset <- new("AnnotatedDataFrame", data = data.frame(fData), varMetadata = fData.meta)

expData <- new("MIAME", name = "Nick",
               lab = "AMPEL", 
               contact = "Nick",title = "Skin",
               abstract = "blah,blah,blah",
               url = "ampel.com",
               other = list(notes = "stuff and things"));

eset <- ExpressionSet(assayData=exprs.eset,
                      phenoData = pData.eset,
                      featureData = fData.eset,
                      experimentData = expData,
                      annotation = "GSE72535 Skin Illumina Eset compiled by AMPEL Biosolutions")

eset
table(eset$cohort)
save(eset, file="GSE72535_Skin_neqc-normalized_NSG_2016-03-03.RData")

#####################################################################################
###Plot PCA

dev.off()
eset$cohort <- factor(eset$cohort)
plotPCA(eset, groups=as.numeric(eset$cohort), groupnames=levels(eset$cohort), main="GSE72535-Skin",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), addtext=sampleNames(eset))

plotPCA(eset, groups=as.numeric(eset$cohort), groupnames=levels(eset$cohort),
        pcs=c(1,2,3), plot3d=TRUE, col=c("red","blue"), pch=c(16,17,17))
