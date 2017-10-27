###################################################################################################
# RECUT WGCNA TOM BLOCKS INTO NEW MODULES USING FOUR LEVELS OF DEEP SPLIT
###################################################################################################

# Load and execute the function to set the storage directory for the recut blocks
# The function will leave you in the endpoint directory, so return to the main project directory as needed
setStorageDir <- dget("/mnt/disks/rserver/projects/shared_scripts/wgcna/set.storage.directory.R")
outputStorage <- setStorageDir( paste("recut",experimentName.recut,sep=".") ) # Each element of this vector is a subsequent subdirectory
# setwd(projectDir) # Don't call until end of job!!

####################################################################################################################
###RECUT THE NETWORK BLOCKS WITH A DETECTION CUT HEIGHT OF 1, MERGE HEIGHT OF 0.2, DEEP SPLIT OF 1
TOMFiles <- paste(projectDir,"/TOM/",net$TOMFiles,sep="")

print(paste("Recutting blocks using deepsplit ","1",sep=""))
recutBlocks = recutBlockwiseTrees(datExpr, net$goodSamples, net$goodGenes, net$blocks, TOMFiles, net$dendrograms,
                                  networkType="signed",
                                  corType = "pearson", pamStage=TRUE, pamRespectsDendro = TRUE,
                                  detectCutHeight = 1,
                                  deepSplit = 1, minModuleSize=100, mergeCutHeight=0.2, numericLabels = TRUE)

print(paste("Deepsplit ",i," cuts complete.",sep=""))

MEs.DS1 <- recutBlocks$MEs
rownames(MEs.DS1) <- rownames(datExpr)

MEs.DS1.labels <- as.numeric(substr( colnames(MEs.DS1), 3, nchar(colnames(MEs.DS1)) ))
colnames(MEs.DS1) <- labels2colors(MEs.DS1.labels)

moduleLabels.DS1 <- recutBlocks$colors
moduleColors.DS1 <- labels2colors(moduleLabels.DS1)

colorLevels.DS1 <- data.frame(table(moduleColors.DS1))
colorLevels.DS1 <- colorLevels.DS1[order(colorLevels.DS1[,2], decreasing=TRUE), ]

####################################################################################################################
## RECUT THE NETWORK BLOCKS WITH A DETECTION CUT HEIGHT OF 1, MERGE HEIGHT OF 0.2, DEEP SPLIT OF 2

print(paste("Recutting blocks using deepsplit ","2",sep=""))
recutBlocks = recutBlockwiseTrees(datExpr, net$goodSamples, net$goodGenes, net$blocks, TOMFiles, net$dendrograms,
                                  networkType="signed",
                                  corType="pearson", pamStage=TRUE, pamRespectsDendro = TRUE,
                                  detectCutHeight = 1,
                                  deepSplit=2, minModuleSize=100, mergeCutHeight=0.2, numericLabels = TRUE)

print(paste("Deepsplit ",i," cuts complete.",sep=""))

MEs.DS2 <- recutBlocks$MEs
rownames(MEs.DS2) <- rownames(datExpr)

MEs.DS2.labels <- as.numeric(substr( colnames(MEs.DS2), 3, nchar(colnames(MEs.DS2)) ))
colnames(MEs.DS2) <- labels2colors(MEs.DS2.labels)

moduleLabels.DS2 <- recutBlocks$colors
moduleColors.DS2 = labels2colors(moduleLabels.DS2)

moduleColors.DS2 <- matchLabels(moduleColors.DS2, moduleColors.DS1, pThreshold = 5e-2, na.rm = TRUE,
                                ignoreLabels = if (is.numeric(moduleColors.DS2)) 0 else "grey",
                                extraLabels = if (is.numeric(moduleColors.DS2)) c(1:1000) else standardColors() )

# Regenerate the MEs now that we have new color names
MEs.DS2 <- moduleEigengenes( datExpr, moduleColors.DS2, verbose=2)
MEs.DS2 <- MEs.DS2$eigengenes
rownames(MEs.DS2) <- rownames(datExpr)
colnames(MEs.DS2) <- substr( colnames(MEs.DS2), 3, nchar(colnames(MEs.DS2)) )

colorLevels.DS2 <- data.frame(table(moduleColors.DS2))
colorLevels.DS2 <- colorLevels.DS2[order(colorLevels.DS2[,2], decreasing=TRUE), ]


####################################################################################################################
## RECUT THE NETWORK BLOCKS WITH A DETECTION CUT HEIGHT OF 1, MERGE HEIGHT OF 0.2, DEEP SPLIT OF 3

print(paste("Recutting blocks using deepsplit ","3",sep=""))
recutBlocks = recutBlockwiseTrees(datExpr, net$goodSamples, net$goodGenes, net$blocks, TOMFiles, net$dendrograms,
                                  networkType="signed",
                                  corType="pearson", pamStage=TRUE, pamRespectsDendro = TRUE,
                                  detectCutHeight = 1,
                                  deepSplit=3, minModuleSize=100, mergeCutHeight=0.2, numericLabels = TRUE)

print(paste("Deepsplit ",i," cuts complete.",sep=""))

MEs.DS3 <- recutBlocks$MEs
rownames(MEs.DS3) <- rownames(datExpr)

MEs.DS3.labels <- as.numeric(substr( colnames(MEs.DS3), 3, nchar(colnames(MEs.DS3)) ))
colnames(MEs.DS3) <- labels2colors(MEs.DS3.labels)

moduleLabels.DS3 <- recutBlocks$colors
moduleColors.DS3 = labels2colors(moduleLabels.DS3)

moduleColors.DS3 <- matchLabels(moduleColors.DS3, moduleColors.DS1, pThreshold = 5e-2, na.rm = TRUE,
                                ignoreLabels = if (is.numeric(moduleColors.DS2)) 0 else "grey",
                                extraLabels = if (is.numeric(moduleColors.DS2)) c(1:1000) else standardColors() )

# Regenerate the MEs now that we have new color names
MEs.DS3 <- moduleEigengenes( datExpr, moduleColors.DS3, verbose=2)
MEs.DS3 <- MEs.DS3$eigengenes
rownames(MEs.DS3) <- rownames(datExpr)
colnames(MEs.DS3) <- substr( colnames(MEs.DS3), 3, nchar(colnames(MEs.DS3)) )

colorLevels.DS3 <- data.frame(table(moduleColors.DS3))
colorLevels.DS3 <- colorLevels.DS3[order(colorLevels.DS3[,2], decreasing=TRUE), ]

####################################################################################################################
## RECUT THE NETWORK BLOCKS WITH A DETECTION CUT HEIGHT OF 1, MERGE HEIGHT OF 0.2, DEEP SPLIT OF 4

print(paste("Recutting blocks using deepsplit ","4",sep=""))
recutBlocks = recutBlockwiseTrees(datExpr, net$goodSamples, net$goodGenes, net$blocks, TOMFiles, net$dendrograms,
                                  networkType="signed",
                                  corType="pearson", pamStage=TRUE, pamRespectsDendro = TRUE,
                                  detectCutHeight = 1,
                                  deepSplit=4, minModuleSize=100, mergeCutHeight=0.2, numericLabels = TRUE)

print(paste("Deepsplit ",i," cuts complete.",sep=""))

MEs.DS4 <- recutBlocks$MEs
rownames(MEs.DS4) <- rownames(datExpr)

MEs.DS4.labels <- as.numeric(substr( colnames(MEs.DS4), 3, nchar(colnames(MEs.DS4)) ))
colnames(MEs.DS4) <- labels2colors(MEs.DS4.labels)

moduleLabels.DS4 <- recutBlocks$colors
moduleColors.DS4 = labels2colors(moduleLabels.DS4)

moduleColors.DS4 <- matchLabels(moduleColors.DS4, moduleColors.DS1, pThreshold = 5e-2, na.rm = TRUE,
                                ignoreLabels = if (is.numeric(moduleColors.DS2)) 0 else "grey",
                                extraLabels = if (is.numeric(moduleColors.DS2)) c(1:1000) else standardColors() )

# Regenerate the MEs now that we have new color names
MEs.DS4 <- moduleEigengenes( datExpr, moduleColors.DS4, verbose=2)
MEs.DS4 <- MEs.DS4$eigengenes
rownames(MEs.DS4) <- rownames(datExpr)
colnames(MEs.DS4) <- substr( colnames(MEs.DS4), 3, nchar(colnames(MEs.DS4)) )

colorLevels.DS4 <- data.frame(table(moduleColors.DS4))
colorLevels.DS4 <- colorLevels.DS4[order(colorLevels.DS4[,2], decreasing=TRUE), ]


#######
# Plot a cut as needed diagnostically
# plotDendroAndColors(net$dendrograms[[1]], moduleColors.DS1[net$blockGenes[[1]]],
#                     c("DS1"), main = "TSOKOS CD4 Active SLE Females DS3 PAM SP10 Cut Height 0.995 Block 1/4",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
#######