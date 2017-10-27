########################################################################################
## CALCULATE KIMs FOR ALL DEEP SPLITS
# Calculate kIM, or intramodular connectivity, i.e. connectivity of nodes to other nodes
# within the same module. This requires the adjaceny matrix be recreated, so make sure
# to use the same soft power and unsigned or signed designation originally used to create the network

setwd(dir)
kIM = matrix(nrow=ncol(datExpr),ncol=4)
for (i in 1:length(net$TOMFiles)) {
  datExpr.block <- datExpr[ , net$blockGenes[[i]] ]
  adjacency.block = adjacency(datExpr.block, power= as.integer(TOM.settings["power"]), type="signed")
  print( paste("Calculating adjacencies for genes within block ",i," of ",length(net$TOMFiles),sep=""))
  moduleColors.block <- moduleColors[ net$blockGenes[[i]] ]
  print( paste("Calculating scaled kIM intramodular connectivity for adjacencies within block ",i," of ",length(net$TOMFiles),sep=""))
  kIM.block = intramodularConnectivity(adjacency.block, moduleColors.block, scaleByMax = TRUE) # Scale to 1
  kIM[1:nrow(kIM.block),1:4] <- kIM.block # FIX
  
}