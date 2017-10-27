setwd("/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/Standardized_CTL_SLEDAI_NA/")
load("GSE10325_CD19_Standardized_CTL_SLEDAI_NA_analysis_parameters.RData")

load("experimentName.recut.RData")
experimentDir<-paste0(analysisDir,experimentName.recut)
setwd(experimentDir)
load("final_results.DS3.RData")
geneInfo<-geneList.DS3.final.simple
colnames(geneInfo)[2]<-"moduleColor"
moduleCor<-moduleCor.DS3
rm(geneList.DS3.final.simple,moduleCor.DS3)

# you will need to make changes throughout this script

###############################################################################################################
# EXPORT MODULES TO EXCEL
options(java.parameters = "-Xmx16000m")
library(xlsx)
wb <- createWorkbook()

title.stem<-paste0(file.base,"_Active_WGCNA_Final-Report_NSG_25Oct2017")

##############################################
# CREATE SUMMARY PAGE

# Include at least:
# 1. GSE number / ArrayExpress number
# 2. Study title
# 3. PubMed ID (if referenced in a publication)
# 4. Type of sample
# 5. How sample was isolated
# 6. Timepoints in the study if there are more than one
# 7. What is being compared (e.g., Disease to healthy; disease to disease)
# 8. Whether or not there is in vitro stimulation of the cells
# 9. Numbers of samples: divided into healthy control and disease type if applicable
# 10. Ethnicity of patients
# 11. Clinical data (e.g., cohort, SLEDAI, drugs, etc.)
# 12. Platform used (e.g., Affy hgU133plus2)
# 13. ID of outliers (ie GSMxxxx)
# 14. Data normalization steps
# 15. Analysis files location
# 16. Version control (e.g., date of analysis and by whom)
# 17. Deepsplit used (i.e., 1, 2, 3, or 4)
# 18. Soft threshold used (e.g., 20)
# 19. Count of starting probes
# 20. Count of retained, analyzed probes
# 21. Correlation function type (i.e., Pearson, Bicor [biweight midcorrelation], or Spearman)
# 22. Signed or unsigned network
# 23. Were duplicates removed post-modularization, and if so, using what collapseRow method
# 24. Block size used in generating SFT and making TOMs

study <- data.frame(nrow=33, ncol=2)
study[1,] <- c("Study Number","GSE49454")
study[2,] <- c("Analysis_Type","WGCNA")
study[3,] <- c("Cell/Tissue Type","Whole Blood")
study[4,] <- c("Cell Stimulation","None")
study[5,] <- c("Common_Name","Chaussabel Whole Blood")
study[6,] <- c("Disease Comparison","Active SLE vs Healthy Controls")
study[7,] <- c("Cohorts","SLE;CTL")
study[8,] <- c("SLE Females","23")
study[9,] <- c("SLE Males","0")
study[10,] <- c("SLE Unknown Sex","0")
study[11,] <- c("Control Females","10")
study[12,] <- c("Control Males","0")
study[13,] <- c("Control Unknown Sex","0")
study[14,] <- c("SLEDAI Range","6 to 26")
study[15,] <- c("Race(s)","Unknown")
study[16,] <- c("Clinical Variables","cohort, SLEDAI, A.DSD, Age, B.cells, C3.Low, C4.Low, CD4.T.cells, CD8.T.cells, CS.mg.d, CS.mg.kg, Eosinophils, HB, HCQ, Leukocytes, Lymphocytes, Neutrophils, Platelets")
study[17,] <- c("Microarray Type","Illumina HT12v4")
study[18,] <- c("Chip Definition Files used","Illumina HT12v4 (October 2016 annotations from GPL10558 GEO platform page)")
study[19,] <- c("Input ESET File","GSE49454_Chaussabel_Whole-Blood_author-normalized_ALL-samples_NSG_2016-03-17.RData")
study[20,] <- c("Outlier Accessions","GSM1199459;GSM1199431;GSM1199437;GSM1199439;GSM1199473;GSM1199492;GSM1199529;GSM1199531")
study[21,] <- c("Starting Probes","47323")
study[22,] <- c("Probe Intensity Threshold","5.25")
study[23,] <- c("Analyzed Probes","16021")
study[24,] <- c("Network Type","Signed")
study[25,] <- c("Correlation Type","Pearson")
study[26,] <- c("Soft Threshold Power","14")
study[27,] <- c("Block Size","10000")
study[28,] <- c("Block Count","1")
study[29,] <- c("Deep Split","3")
study[30,] <- c("Analysis Performed by","NSG")
study[31,] <- c("Analysis Conducted on","21March2017")
study[32,] <- c("Source Publication PMID(s)","24644022")
study[33,] <- c("Notes","1st visit samples ONLY")

ws <- createSheet(wb=wb, sheetName="WGCNA Summary")
addDataFrame(x=study, sheet=ws, row.names=FALSE,col.names=FALSE)

######################################################################
# ADD A TABLE OF CORRELATIONS FOR MODULE EIGENGENES AND CLINICAL TRAITS

color_index <- levels(factor(geneInfo$moduleColor))
length(color_index)

library(dplyr)

geneInfo.list <- split(geneInfo, f = geneInfo$moduleColor) # Convert to a list object with elements based on module colors
data.list <- lapply(geneInfo.list,function(x){ # Remove all within module duplicates by taking the ones with highest kWithin
  as.data.frame(x) %>%
    group_by(geneEntrezID) %>%
    top_n(1,kWithin)})

library (plyr)

df <- ldply (data.list, data.frame) # Turn list back into a dataframe
colnames(df)[4]<-"ensembl"
df <- df[,-1]
rownames(df) <- df$ensembl #Now it looks just like geneInfo
df <- df[!is.na(df$geneEntrezID),] # Remove all unannotated probes
df<-df[!is.na(df$geneSymbol),]

ws <- createSheet(wb=wb, sheetName="Correlations")
moduleCor <- as.data.frame(moduleCor)
moduleCor <- moduleCor[order(rownames(moduleCor)),]
df$moduleColor <- as.character(df$moduleColor)
mod.counts <- as.data.frame(table(df$moduleColor))
mod.counts <- mod.counts[order(mod.counts$Var1),]
rownames(mod.counts) <- mod.counts$Var1
moduleCor <- moduleCor[which(rownames(moduleCor)!="grey"),]
identical(rownames(moduleCor),rownames(mod.counts))
moduleCor$moduleCount <- mod.counts$Freq

addDataFrame(x=moduleCor, sheet=ws, row.names=TRUE,col.names=TRUE)

######################################################################
# NOW ADD INDIVIDUAL MODULE SHEETS IN WORKBOOK

addDataFrame(x=moduleCor, sheet=ws, row.names=TRUE,col.names=TRUE)
for (i in 1:length(color_index) ) {
  single.mods <- df[ df$moduleColor==color_index[i], ]
  single.mods <- single.mods[ !is.na(single.mods$moduleColor), ]
  single.mods <- single.mods[ order(single.mods$SLEDAI.p), ] ######################MODIFY the primary clinical trait
  nrow(single.mods)
  ws <- createSheet(wb=wb, sheetName=paste(toupper(color_index[i])))
  addDataFrame(x=single.mods, sheet=ws, row.names=FALSE)
}

######################################################################
# SAVE THE WORKBOOK
setwd(experimentDir)
saveWorkbook(wb, paste0(title.stem,".xlsx"))
write.table(geneInfo, paste0(title.stem,".txt"))

