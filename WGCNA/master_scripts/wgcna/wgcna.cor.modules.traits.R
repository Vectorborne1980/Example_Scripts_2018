
traitCorr <- function( MEs, datTraits, moduleColors, sigTrait, p_val,
                       title.heatmap,
                       plotMar, cex.lab.x, cex.lab.y, text.size,
                       xLabels, xColorOffset, xLabelsAngle,
                       maxRows,
                       chart.filename, pngWidth, pngHeight, saveFile, deepSplit.cor) {
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  datTraits <- data.matrix(datTraits)
  
  # Use the R cor() function to correlate the module eigengenes to clinical traits along with an associated p-value
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  moduleTraitCor <- data.frame(moduleTraitCor)
  moduleTraitPvalue <- data.frame(moduleTraitPvalue)
  
  corNames <- vector('character')
  moduleCor <- vector('numeric')
  
  for (i in 1:ncol(moduleTraitCor)) {
    moduleCor <- cbind( moduleCor, moduleTraitCor[,i],moduleTraitPvalue[,i] )
    corNames <- c( corNames, colnames(moduleTraitCor)[i], paste(colnames(moduleTraitPvalue)[i],".p",sep="") )
  }
  
  colnames(moduleCor) <- corNames
  rownames(moduleCor) <- rownames(moduleTraitCor)
  
  # Order the correlations frames by decreasing significance of correlation to the primary clinical variable of choice,
  # usually SLEDAI or cohort
  moduleCor <- moduleCor[ order(moduleCor[, paste(sigTrait,".p",sep="")] ), ]
  
  index <- match( rownames(moduleCor),rownames(moduleTraitCor) )
  moduleTraitCor <- moduleTraitCor[index,]
  
  index <- match( rownames(moduleCor),rownames(moduleTraitPvalue) )
  moduleTraitPvalue <- moduleTraitPvalue[ index, ]
  
  moduleTraitCor <- data.matrix(moduleTraitCor)
  moduleTraitPvalue <- data.matrix(moduleTraitPvalue)
  
  # Sum the number of probes in each module and bind this to our correlation matrix
  #table(moduleColors)
  colorFreq <- data.frame(table(factor(moduleColors)))
  rownames(colorFreq) <- colorFreq[,1]
  index <- match(rownames(moduleCor),rownames(colorFreq) )
  colorFreq <- colorFreq[index,]
  moduleCor <- cbind( moduleCor, "moduleCount"=colorFreq$Freq )
  #colnames(moduleCor)
  moduleCor <- moduleCor[ ,c(ncol(moduleCor),1:ncol(moduleCor)-1) ]
  
  index <- which(colnames(moduleCor)==paste(sigTrait,".p",sep=""))
  sig.count <- length(which(moduleCor[,index] < p_val))
  
  # Prepare plot window
  sizeGrWindow(10,6)
  # Color code each association by the correlation value:
  # Display the correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), " (",
                     signif(moduleTraitPvalue, 2), ")", sep = "");
  
  if (saveFile==TRUE) {
    # Generate a PNG of the figure
    png(chart.filename,pngWidth, pngHeight, pointsize=20)
  }
  # dev.off()
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = plotMar); # bottom, left, top, right margins
  
  # Display the correlation values using a heatmap plot
  yLabels <- paste(rownames(moduleTraitCor)," (",moduleCor[,1], ")", sep="")
  ySymbols <- rownames(moduleTraitCor)
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = xLabels,
                 xLabelsAngle = xLabelsAngle,
                 xColorOffset = xColorOffset,
                 yLabels = yLabels,
                 ySymbols = ySymbols,
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = text.size,
                 zlim = c(-1,1),
                 cex.lab.x = cex.lab.x,
                 cex.lab.y = cex.lab.y,
                 main = paste(title.heatmap,
                              ", ", sig.count, " module(s) sig by ",sigTrait," p_val < ",p_val, " (", 
                              round(sig.count/nrow(moduleCor)*100,0),"%)", sep=""))
  if (saveFile==TRUE) {
    # Complete PNG save
    dev.off()
    print(paste( deepSplit.cor,"dendrogram with trait bars saved."))
  }
  
  moduleTraitCor <- data.matrix(moduleTraitCor)
  moduleTraitPvalue <- data.matrix(moduleTraitPvalue)
  
  moduleCor.list <- list(moduleTraitCor,moduleTraitPvalue, moduleCor )
  
  return(moduleCor.list)
  
  plotMar <- c(8, 6, 4, 2) # bottom, left, top, and right figure margins
  
  #rm(plotMar, text.size, p_val, x.size, y.size, maxRows, title.heatmap, chart.filename)
  
}

# Correlate modules created during TOM generation
moduleCor <- traitCorr( MEs=MEs, datTraits=datTraits, moduleColors=moduleColors, sigTrait=sigTrait, p_val=p_val,
                            title.heatmap=paste("DS Original",experimentName.recut,sep=" "),
                            plotMar, cex.lab.x, cex.lab.y, text.size,
                            xLabels, xColorOffset, xLabelsAngle,
                            maxRows=maxRows,chart.filename=paste("DS-orig.corrMap.",experimentName.recut,".png",sep=""),
                            pngWidth, pngHeight,saveFile, deepSplit.cor="DS-orig" )
moduleCor.val <- moduleCor[[1]]
moduleCor.pval <- moduleCor[[2]]
moduleCor <- moduleCor[[3]]

# Correlate deep split 1 modules
moduleCor.DS1 <- traitCorr( MEs=MEs.DS1, datTraits=datTraits, moduleColors=moduleColors.DS1, sigTrait=sigTrait, p_val=p_val,
                            title.heatmap=paste("DS1",experimentName.recut,sep=" "),
                            plotMar, cex.lab.x, cex.lab.y, text.size,
                            xLabels, xColorOffset, xLabelsAngle,
                            maxRows=maxRows,chart.filename=paste("DS1.corrMap.",experimentName.recut,".png",sep=""),
                            pngWidth, pngHeight,saveFile, deepSplit.cor="DS1" )
moduleCor.DS1.val <- moduleCor.DS1[[1]]
moduleCor.DS1.pval <- moduleCor.DS1[[2]]
moduleCor.DS1 <- moduleCor.DS1[[3]]

# Correlate deep split 2 modules
moduleCor.DS2 <- traitCorr( MEs=MEs.DS2, datTraits=datTraits, moduleColors=moduleColors.DS2, sigTrait=sigTrait, p_val=p_val,
                            title.heatmap=paste("DS2",experimentName.recut,sep=" "),
                            plotMar=plotMar, cex.lab.x=cex.lab.x, cex.lab.y = cex.lab.y, text.size=text.size,
                            xLabels=xLabels, xColorOffset=xColorOffset, xLabelsAngle=xLabelsAngle,
                            maxRows=maxRows,
                            chart.filename=paste("DS2.corrMap.",experimentName.recut,".png",sep=""), pngWidth=pngWidth, pngHeight=pngHeight,
                            saveFile=saveFile, deepSplit.cor="DS2"  )
moduleCor.DS2.val <- moduleCor.DS2[[1]]
moduleCor.DS2.pval <- moduleCor.DS2[[2]]
moduleCor.DS2 <- moduleCor.DS2[[3]]

# Correlate deep split 3 modules
moduleCor.DS3 <- traitCorr( MEs=MEs.DS3, datTraits=datTraits, moduleColors=moduleColors.DS3, sigTrait=sigTrait, p_val=p_val,
                            title.heatmap=paste("DS3",experimentName.recut,sep=" "),
                            plotMar=plotMar, cex.lab.x=cex.lab.x, cex.lab.y = cex.lab.y, text.size=text.size,
                            xLabels=xLabels, xColorOffset=xColorOffset, xLabelsAngle=xLabelsAngle,
                            maxRows=maxRows,
                            chart.filename=paste("DS3.corrMap.",experimentName.recut,".png",sep=""), pngWidth=pngWidth, pngHeight=pngHeight,
                            saveFile=saveFile, deepSplit.cor="DS3"  )
moduleCor.DS3.val <- moduleCor.DS3[[1]]
moduleCor.DS3.pval <- moduleCor.DS3[[2]]
moduleCor.DS3 <- moduleCor.DS3[[3]]

# Correlate deep split 4 modules
moduleCor.DS4 <- traitCorr( MEs=MEs.DS4, datTraits=datTraits, moduleColors=moduleColors.DS4, sigTrait=sigTrait, p_val=p_val,
                            title.heatmap=paste("DS4",experimentName.recut,sep=" "),
                            plotMar=plotMar, cex.lab.x=cex.lab.x, cex.lab.y = cex.lab.y, text.size=text.size,
                            xLabels=xLabels, xColorOffset=xColorOffset, xLabelsAngle=xLabelsAngle,
                            maxRows=maxRows,
                            chart.filename=paste("DS4.corrMap.",experimentName.recut,".png",sep=""), pngWidth=pngWidth, pngHeight=pngHeight,
                            saveFile=saveFile, deepSplit.cor="DS4"  )
moduleCor.DS4.val <- moduleCor.DS4[[1]]
moduleCor.DS4.pval <- moduleCor.DS4[[2]]
moduleCor.DS4 <- moduleCor.DS4[[3]]

# Set the proper save directory first. script needs fixing. 
save(moduleCor.DS1, moduleCor.DS1.pval, moduleCor.DS1.val,
     moduleCor.DS2, moduleCor.DS2.pval, moduleCor.DS2.val,
     moduleCor.DS3, moduleCor.DS3.pval, moduleCor.DS3.val,
     moduleCor.DS4, moduleCor.DS4.pval, moduleCor.DS4.val,
     file="recut.DS1-DS4.module.correlations.RData")
