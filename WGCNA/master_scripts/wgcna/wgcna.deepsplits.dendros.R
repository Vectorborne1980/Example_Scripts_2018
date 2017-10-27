
plotDendro <- function( datExpr, moduleColors, moduleCor, sigTrait, parMar, p_val, cex, deepSplit,
                        chart.filename, pngWidth, pngHeight, saveFile ) {
  
  ######################################################
  # Prepare p-values of correlation between modules eigengenes and clinical traits
  sigTrait.p <-paste(sigTrait,".p",sep="")
  
  # Prepare module eigengene names
  MEList = moduleEigengenes(datExpr, colors = moduleColors)
  MEs = MEList$eigengenes
  
  ######################################################
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  
  # Cluster module eigengene dissimilarity
  METree = hclust(as.dist(MEDiss), method = "average");
  
  ######################################################
  # Format module names and include members and significance
  METree$labels <- substr( METree$labels, 3, nchar(METree$labels) )
  index <- match(METree$labels, rownames(moduleCor) )
  tree.labels <- data.frame(moduleCor[index,c("moduleCount",sigTrait.p,sigTrait)])
  tree.labels$color <- rownames(moduleCor)
  tree.labels[sigTrait.p] <- round(tree.labels[sigTrait.p],5)
  tree.labels$sign <- "[+]"
  index <- which(tree.labels[sigTrait]<=0)
  tree.labels[index,"sign"] <- "[-]"
  tree.labels$label <- paste("* ",tree.labels$sign," ", rownames(tree.labels)," (",tree.labels$moduleCount," p",tree.labels[,2],")",sep="")
  index <- which(tree.labels[sigTrait.p] > p_val)
  #tree.labels[index,"label"] <- paste(rownames(tree.labels)[index],"(",tree.labels[index,"moduleCount"],")",sep="")
  tree.labels[index,"label"] <- paste( tree.labels[index,"sign"]," ", rownames(tree.labels)[index]," (",tree.labels[index,"moduleCount"],")", sep="" )
  METree$labels <- tree.labels$label
  sigNum <- length(which(tree.labels[,2]<p_val))
  
  ######################################################
  # Plot the result
  if (!.Device=="null device") { dev.off() } # Catches the "cannot shut down device 1 (the null device)" error
  if (saveFile==TRUE) {
    # Generate a PNG of the figure
    png(chart.filename,pngWidth, pngHeight, pointsize=20)
  }
  
  par(oma=c(0,0,0,0),mar=parMar) # bottom, left, top, right margins
  #sizeGrWindow(7, 6)
  dendro.title <- paste(projectName," ",deepSplit," ",experimentName.recut," p_val<",p_val," *",sigNum," sig ",sigTrait,sep="")
  
  # Old vertical style denrogram
  #plot(METree, main = dendro.title, xlab = "", sub = "", cex=cex)
  
  # Plot a nice horizontal plot
  hcd <- as.dendrogram(METree)
  par(oma=c(0,0,0,0),mar=parMar)
  nodePar <- list(lab.cex = cex, pch = c(NA, 21), cex = cex, col = "black")
  plot(hcd,  xlab = "Module Dissimilarity", nodePar = nodePar, horiz = TRUE, main = dendro.title, cex=cex)
  
  # Complete PNG save
  if (!.Device=="null device") { dev.off() }
  
  print(paste(deepSplit," dendrogram of ",ncol(MEs)," modules generated.",sep=""))
  
  ######################################################
  # NEW DENDROGRAM APPROACHES
  # See http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
  
  # Compute distances and hierarchical clustering
  #dd <- dist(MEDiss, method = "euclidean") # note we've scaled this
  # dd <- dist(scale(MEDiss), method = "euclidean") # note we've scaled this
  # hc <- hclust(dd, method = "ward.D2") # consider using other methods such as euclidean, etc.
  # 
  # dev.off()
  # plot(hc)
  # 
  # hc$labels <- tree.labels$label
  # plot(hc)
  # 
  # # Put the labels at the same height: hang = -1
  # plot(hc, hang = -1, cex = 0.6)
  # 
  # # Convert hclust into a dendrogram and plot
  # hcd <- as.dendrogram(hc)
  # 
  # # Default plot
  # dev.off()
  # par(oma=c(0,0,0,0),mar=c(8,0,0,0)) # bottom, left, top, right margins
  # plot(hcd, type = "rectangle", ylab = "Height")
  # 
  # # Triangle plot
  # dev.off()
  # par(oma=c(0,0,0,0),mar=c(8,0,0,0))
  # plot(hcd, type = "triangle", ylab = "Height")
  # 
  # # Zoom in to the first dendrogram
  # plot(hcd, xlim = c(1, 20), ylim = c(1,8))
  # 
  # # Define nodePar
  # nodePar <- list(lab.cex = 0.8, pch = c(NA, 19), cex = 0.9, col = "black")
  # # Customized plot; remove labels
  # plot(hcd, ylab = "Height", nodePar = nodePar, leaflab = "none")
  # 
  # # Horizontal plot
  # library(dendextend)
  # hcd <- as.dendrogram(hc)
  # 
  # dev.off()
  # par(oma=c(0,0,0,0),mar=c(6,4,2,14))
  # nodePar <- list(lab.cex = 0.85, pch = c(NA, 21), cex = 0.5, col = "black")
  # 
  # plot(hcd,  xlab = "Module Dissimilarity", nodePar = nodePar, horiz = TRUE, main = dendro.title)
  # 
  # # labels_colors(hcd) <- tree.labels$color # Not working yet, need to match colors before reorganizing
  # 
  # # The package ape (Analyses of Phylogenetics and Evolution) can be used to produce a more sophisticated dendrogram.
  # library("ape")
  # # Default plot
  # dev.off()
  # par(oma=c(0,0,0,0),mar=c(4,2,2,2))
  # plot(as.phylo(hc), cex = 0.8, label.offset = 0.5)
  # plot(as.phylo(hc), cex = 0.8, label.offset = 0.5, type="cladogram")
  # 
  # # Unrooted phylogenetic tree
  # plot(as.phylo(hc), type = "unrooted", cex = 0.6, no.margin = TRUE)
  # 
  # # Fan
  # plot(as.phylo(hc), type = "fan")
  # 
  # # GGplot2 dendrograms
  # library("ggplot2", "ggdendro")
  # # Visualization using the default theme named theme_dendro()
  # ggdendrogram(hc)
  # # Rotate the plot and remove default theme
  # ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)
  # # Build dendrogram object from hclust results
  # dend <- as.dendrogram(hc)
  # # Extract the data (for rectangular lines)
  # # Type can be "rectangle" or "triangle"
  # dend_data <- dendro_data(dend, type = "rectangle")
  # # What contains dend_data
  # names(dend_data)
  # # Extract data for line segments
  # head(dend_data$segments)
  # # Extract data for labels
  # head(dend_data$labels)
  # # Plot line segments and add labels
  # p <- ggplot(dend_data$segments) + 
  #   geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  #   geom_text(data = dend_data$labels, aes(x, y, label = label),
  #             hjust = 1, angle = 90, size = 3)+
  #   ylim(-3, 15)
  # print(p)
  
  # Phylotools package
  # See http://www.phytools.org/ and https://cran.r-project.org/web/packages/phytools/
  # Unsure of relevancy to distance dendrograms
  
}

