# Load analysis libraries into memory. Run at the beginning of each new R session
library(GEOquery)
library(affy)
library(simpleaffy) 
library(affycoretools) 

# set working directories
dir <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data"
eset <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/eset"
patient_data <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/patient_data"
raw <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/raw"

#--------- LOAD PATIENT DATA ----------------------------------------------------------------

setwd(patient_data)
pdata <- read.table("GSE32591_patient-data.txt", header=TRUE, sep="\t")
rownames(pdata) <- pdata[,1]
pData <- pdata[,-1]

pData$cohort <- factor(pData$cohort)
pData$cohort_full <- factor(pData$cohort_full)
pData$cell_type <- factor(pData$cell_type)
pData$WHO_class <- factor(pData$WHO_class)

#--------- CREATE BATCH FILES ---------------------------------------------------------------
setwd(raw)
rownames(pData) <- paste( rownames(pData),".CEL",sep="" )

# AFFY
library(hgu133a.db)
abatch <- ReadAffy(verbose=TRUE, celfile.path=raw, phenoData=pData)
colnames(abatch)
par(mar=c(7,5,1,1))
boxplot(abatch,las=2,cex.axis=0.5,main="GSE32591 Probe Intensities BEFORE GCRMA_Affymetrix")
# GCRMA (Guanine Cytosine Robust Multi-Array Analysis) adjusts  for  background  intensities  in  A ymetrix  array  data  which  include optical
# noise and non-speci c binding (NSB). The main function gcrma converts background-adjusted probe intensities to expression measures using the
# same normalization and summarization methods as rma.
Aesetgcrma<-gcrma(abatch)
dev.off()
par(mar=c(7,5,1,1))
boxplot(exprs(Aesetgcrma),las=2,cex.axis=0.5,main="GSE32591 Probe Intensities AFTER GCRMA_Affymetrix")


# BRAINARRAY
library(hgu133ahsentrezgcdf)
library(hgu133ahsentrezgprobe)
cdf<-"HGU133A_Hs_ENTREZG" # this from the BrainArray website where the CDFs were obtained
bbatch <- ReadAffy(verbose=TRUE, celfile.path=raw, cdfname=cdf, phenoData=pData)
colnames(bbatch)
dev.off()
par(mar=c(7,5,1,1))
boxplot(bbatch,las=2,cex.axis=0.5,main="GSE32591 Probe Intensities BEFORE GCRMA_Brainarray")
Besetgcrma <- gcrma(bbatch)
dev.off()
par(mar=c(7,5,1,1))
boxplot(exprs(Besetgcrma),las=2,cex.axis=0.5,main="GSE32591 Probe Intensities AFTER GCRMA_Brainarray")


# SAVE THE AFFYBATCH OBJECTS
setwd(eset)
save(abatch,file="GSE32591_abatch_NSG_2016-07-11.RData")
save(bbatch,file="GSE32591_bbatch_NSG_2016-07-11.RData")

#---------- ANNOTATION -------------------------------------------------------------------

# FUNCTION

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=fromKey, col=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

# AFFY

library(hgu133a.db)
fData <- data.frame(probe=rownames(Aesetgcrma)) # Probe IDs
rownames(fData) <- fData$probe
keys<-rownames(fData(Aesetgcrma))
geneSymbol <- data.frame(geneSymbol=convertIDs( keys, "PROBEID", "SYMBOL", hgu133a.db ))
# Ignore the "'select()' returned 1:many mapping between keys and columns" warning message, just ensure the number of objects in each conversion list is the same length
geneName <- data.frame(geneName=convertIDs( keys, "PROBEID", "GENENAME", hgu133a.db ))
geneEntrezID <- data.frame(geneEntrezID=convertIDs( keys, "PROBEID", "ENTREZID", hgu133a.db ))
geneEnsembl <- data.frame(geneEnsembl=convertIDs( keys, "PROBEID", "ENSEMBL", hgu133a.db ))
fDataComplete <- cbind( fData, geneSymbol, geneName, geneEntrezID, geneEnsembl  )

fData(Aesetgcrma) <- fDataComplete

# BRAINARRAY

library(hgu133ahsentrezg.db)
columns(hgu133ahsentrezg.db) # list available gene annotations

fData <- data.frame(probe=rownames(Besetgcrma)) # the probes are obviously the same in all the new esets
rownames(fData) <- fData$probe
keys<-rownames(fData(Besetgcrma))
geneSymbol <- data.frame(geneSymbol=convertIDs( keys, "PROBEID", "SYMBOL", hgu133ahsentrezg.db ))
geneName <- data.frame(geneName=convertIDs( keys, "PROBEID", "GENENAME", hgu133ahsentrezg.db ))
geneEntrezID <- data.frame(geneEntrezID=convertIDs( keys, "PROBEID", "ENTREZID", hgu133ahsentrezg.db ))
geneEnsembl <- data.frame(geneEnsembl=convertIDs( keys, "PROBEID", "ENSEMBL", hgu133ahsentrezg.db ))
fDataComplete <- cbind( fData, geneSymbol, geneName, geneEntrezID, geneEnsembl  )

fData(Besetgcrma) <- fDataComplete

# REMOVE AFFY PROBES

length(rownames(Aesetgcrma)) #22283
index <- which( substr(rownames(Aesetgcrma),1,4)=="AFFX" )
Aesetgcrma <- Aesetgcrma[ -index, ]
length(rownames(Aesetgcrma)) #22215

length(rownames(Besetgcrma))
index <- which( substr(rownames(Besetgcrma),1,4)=="AFFX" )
Besetgcrma <- Besetgcrma[ -index, ]
length(rownames(Besetgcrma))
rm(index)

# SAVE ESETS

setwd(eset)
save(Aesetgcrma,file="GSE32591_Aesetgcrma_NSG_2016-07-11.RData")
save(Besetgcrma,file="GSE32591_Besetgcrma_NSG_2016-07-11.RData")
