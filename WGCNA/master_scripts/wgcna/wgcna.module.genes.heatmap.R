###################################################################################################
# PLOT A HEATMAP OF TRANSCRIPT EXPRESSION PER MODULE
# 08-02-17 RDR: Function will generate a heatmap with the file name containing deepsplit and pvalue 
# of correlation of most interesting clinical trait to module eigengene. Currently I'm
# not adding a color bar, which we once used for manually-set kmeans number of clusters
###################################################################################################

library(Biobase)

heatmap.module <- function( moduleNames, moduleCor, geneInfo, eset, saveFile, deepSplit,
                            cexRow, cexCol, margins, pngWidth, pngHeight) {
  
  for ( i in 1:length(moduleNames) ) {
    geneInfo.subset <- geneInfo[ geneInfo$moduleColor==moduleNames[i], ]
    index <- which(rownames(eset)%in%geneInfo.subset$probe)
    m.use <- exprs(eset)[index,]
    distance <- dist(m.use,method="euclidean")
    hc <- hclust(distance)
    
    hc.cut = cutree(hc,k=6) # k = number of clusters. 08-02-17 RDR: Do we still need to do this? If so, pass k as parameter
    #table(hc.cut) # show cluster cut results
    hc.cut.levels = levels(as.factor(hc.cut)) # factor cluster level numbers for coloring
    #hc.cut.levels # review cluster level factors
    
    colors = rainbow(length(hc.cut.levels));
    cluster.patient.colors = c()
    
    #for (i in 1:length(hc.cut.levels)) {
      #cluster.patient.colors[names(which(hc.cut==hc.cut.levels[i]))] = colors[i] 
      #}
    #cluster.patient.colors
    
    dist.euclidean = dist.euclidean <- function(x, ...)
      dist(x,method="euclidean")
    label.size = 1
    
    samples = colnames(m.use)
    
    hmcol <- colorRampPalette(c("blue","white","red"))(256) #FOR THE COLORBLIND AMONG US
    
    ## MAKE SURE TO MAKE THE R STUDIO PLOTS FRAME AS LARGE AS POSSIBLE BEFORE RUNNING (SHRINK THE OTHER FRAMES)
    # Plot the result
    if (!.Device=="null device") { dev.off() } # Catches the "cannot shut down device 1 (the null device)" error
    
    # Grab the module's correlation significance to sigTrait
    module.sig <- round(moduleCor[moduleNames[i],paste(sigTrait,".p",sep="")],3)
    
    if (saveFile==TRUE) {
      # Generate a PNG of the figure
      png( paste(deepSplit,".",moduleNames[i],".",sigTrait,".p",module.sig,".png",sep=""), pngWidth, pngHeight, pointsize=20)
    }
    
    # genes.heat = heatmap.2(data.matrix(m.use), col=hmcol, scale="row", distfun = dist2, key=TRUE,
    #                        symkey=TRUE, density.info="none", trace="none",
    #                        ColSideColors=cluster.patient.colors[samples],cexCol=1.2,
    #                        cexRow = 0.1,margins=c(12,9),
    #                        main=paste(moduleNames[i]," module: ",projectName,".", deepSplit, " gene expression", sep=""))
    
    # Plot heatmap without side color bar
    genes.heat = heatmap.2(data.matrix(m.use), col=hmcol, scale="row", distfun = dist.euclidean, key=TRUE,
                           symkey=TRUE, density.info="none", trace="none",
                           cexRow = cexRow, cexCol=cexCol, margins=margins,
                           main=paste(moduleNames[i]," ", sigTrait," p",module.sig,
                                      " module ",deepSplit,".",projectName, sep=""))
    #return(m.use)
    
    if (saveFile==TRUE) {
      # Complete PNG save
      dev.off()
      print(paste(deepSplit,moduleNames[i], "heatmap generated",sep=" "))
    }
  }
  #print(paste(deepSplit," heatmaps of ",length(moduleNames)," modules generated.",sep=""))
}
