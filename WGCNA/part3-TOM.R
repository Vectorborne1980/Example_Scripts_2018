
# You will need to change the power and update the saveTOMFileBase with the current date

### WGCNA ADJACENCY MATRIX GENERATION AND TOM ENCODING
#
# 07-19-17 RDR: maxBlockSize should always be set to 10,000 despite the resource of the underlying virtual machine
# In general, the only parameters to be changed here are file naming and the power selection.
# Also, this script must be run from the VM terminal, i.e. $sudo Rscript part3-TOM.R

library(WGCNA)
options(stringsAsFactors = FALSE);
allowWGCNAThreads()
setwd("/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/Standardized_CTL_SLEDAI_NA/")
load("GSE10325_CD19_Standardized_CTL_SLEDAI_NA_analysis_parameters.RData")
load("part1-final-objects.RData")

corType = "pearson"
power = 22
networkType = "signed"
TOMType = "signed"
deepSplit = 3
minModuleSize = 100
numericLabels = TRUE
detectCutHeight = 1
mergeCutHeight = 0.2
pamStage = TRUE
pamRespectsDendro = TRUE
saveTOMFileBase = paste0(file.base,"_SP22_10K_25Oct2017")

TOM.settings <- c(power, networkType, TOMType, deepSplit, minModuleSize, numericLabels, detectCutHeight, mergeCutHeight, 
                  pamStage, pamRespectsDendro, saveTOMFileBase)
names(TOM.settings) <- c("power", "networkType", "TOMType", "deepSplit", "minModuleSize", "numericLabels", "detectCutHeight",
                         "mergeCutHeight", "pamStage", "pamRespectsDendro", "saveTOMFileBase")
save(TOM.settings,file=paste0(saveTOMFileBase,"_TOM_settings.RData"))

net = blockwiseModules(datExpr,
                       power = power,
                       corType=corType,
                       deepSplit = deepSplit,
                       networkType = "signed", # ALWAYS SIGNED
                       TOMType = "signed", # ALWAYS SIGNED
                       minModuleSize = minModuleSize,
                       numericLabels = TRUE,
                       detectCutHeight = detectCutHeight,
                       mergeCutHeight = mergeCutHeight,
                       pamStage = pamStage,
                       pamRespectsDendro = pamRespectsDendro, 
                       saveTOMs = TRUE,
                       saveTOMFileBase = saveTOMFileBase,
                       verbose = 5,
                       maxBlockSize = 10000) # ALWAYS 10,000

save(net,file=paste0(saveTOMFileBase,"_NET.RData"))
