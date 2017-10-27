###################################################################################################
# RECUT WGCNA TOM BLOCKS INTO NEW MODULES USING FOUR LEVELS OF DEEP SPLIT
###################################################################################################

# Load and execute the function to set the storage directory for the recut blocks
# The function will leave you in the endpoint directory, so return to the main project directory as needed
# setStorageDir <- dget("/mnt/disks/rserver/projects/shared_scripts/wgcna/set.storage.directory.R")
# outputStorage <- setStorageDir( experimentName.recut ) # Each element of this vector is a subsequent subdirectory
# setwd(projectDir) # Don't call until end of job!!

for (i in 1:4) {
  ### CUT WITH DEEPSPLIT 1 THROUGH 4
  recutBlocks = recutBlockwiseTrees(datExpr, net$goodSamples, net$goodGenes, net$blocks,
                                    TOMFiles=paste(projectDir,"/",net$TOMFiles,sep=""), net$dendrograms,
                                    networkType="signed",TOMType = "signed",
                                    corType = recut.corType,
                                    pamStage=recut.pamStage, pamRespectsDendro=recut.pamRespectsDendro,
                                    detectCutHeight=recut.detectCutHeight,
                                    deepSplit=i, minModuleSize=recut.minModuleSize, mergeCutHeight=recut.mergeCutHeight,
                                    numericLabels = TRUE)
  print(paste("Deepsplit ",i," cuts complete.",sep=""))
  
  # Use these in a potentially simplified for-loop instead of the four if(i==..) blocks. 
  # Note usage of assign() and as.name() functions
  # assign( paste("MEs.DS",i,sep=""),recutBlocks$MEs )
  # rownames(as.name( paste("MEs.DS",i,sep="")))  <- rownames(datExpr)
  
  if (i==1) {
    print(paste("Recutting blocks using deepsplit ",i,sep=""))
    MEs.DS1 <- recutBlocks$MEs
    rownames(MEs.DS1) <- rownames(datExpr)
    MEs.DS1.labels <- as.numeric(substr( colnames(MEs.DS1), 3, nchar(colnames(MEs.DS1)) ))
    colnames(MEs.DS1) <- labels2colors(MEs.DS1.labels)
    moduleLabels.DS1 <- recutBlocks$colors
    moduleColors.DS1 <- labels2colors(moduleLabels.DS1)
    colorLevels.DS1 <- data.frame(table(moduleColors.DS1))
    colorLevels.DS1 <- colorLevels.DS1[order(colorLevels.DS1[,2], decreasing=TRUE), ]
    print(paste("Block ",i," recut modules regenerated.",sep=""))
  }
  
  if (i==2) {
    print(paste("Recutting blocks using deepsplit ",i,sep=""))
    MEs.DS2 <- recutBlocks$MEs
    rownames(MEs.DS2) <- rownames(datExpr)
    
    MEs.DS2.labels <- as.numeric(substr( colnames(MEs.DS2), 3, nchar(colnames(MEs.DS2)) ))
    colnames(MEs.DS2) <- labels2colors(MEs.DS2.labels)
    
    moduleLabels.DS2 <- recutBlocks$colors
    moduleColors.DS2 = labels2colors(moduleLabels.DS2)
    
    # Here we being using matchLabels(), where the modules created by DS2 are matched to those created by DS1.
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
    print(paste("Block ",i," recut modules regenerated.",sep=""))
  }
  
  if (i==3) {
    print(paste("Recutting blocks using deepsplit ",i,sep=""))
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
    print(paste("Block ",i," recut modules regenerated.",sep=""))
  }
  
  if (i==4) {
    print(paste("Recutting blocks using deepsplit ",i,sep=""))
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
    print(paste("Block ",i," recut modules regenerated.",sep=""))
  }
}

setwd(projectDir)

#######
# Plot a cut as needed diagnostically
# plotDendroAndColors(net$dendrograms[[1]], moduleColors.DS1[net$blockGenes[[1]]],
#                     c("DS1"), main = "TSOKOS CD4 Active SLE Females DS3 PAM SP10 Cut Height 0.995 Block 1/4",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
#######