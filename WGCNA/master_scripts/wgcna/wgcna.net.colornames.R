###################################################################################################
# ADD COLOR NAMES TO WGCNA MODULES
###################################################################################################

# Load probe/gene numeric module assignments
moduleLabels <- net$colors

# Convert numeric module numbers to color names
moduleColors <- labels2colors(net$colors)

# Extract module eigengenes per sample
MEs <- net$MEs
rownames(MEs) <- rownames(datExpr)

# Rename module numbers to colors 
# CRITICAL!!! You MUST pass a numeric vector to labels2colors
MEs.labels <- as.numeric(substr( colnames(MEs), 3, nchar(colnames(MEs)) ))
colnames(MEs) <- labels2colors(MEs.labels)

# Create frequency table of probes in modules
colorLevels <- data.frame(table(moduleColors))
colorLevels <- colorLevels[order(colorLevels[,2], decreasing=TRUE), ]

# Match row names of datExpr and datTraits
datExpr <- datExpr[ order(rownames(datExpr)), ]
datTraits <- datTraits[ order(rownames(datTraits)), ]
identical( rownames(datExpr), rownames(datTraits) )
