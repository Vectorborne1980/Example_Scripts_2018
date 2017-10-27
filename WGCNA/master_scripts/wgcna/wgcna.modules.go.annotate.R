###################################################################################################
# ANNOTATE WGCNA MODULES WITH TOP GO PATHWAYS
###################################################################################################

# Join some fData to genes with their module color names. Note the gene order in the moduleColors vector is
# still in the same order as rownames of the original datExpr from whence adjacency and subsequent
# TOM encoding was calculated.

geneNames = colnames(datExpr)
geneInfo = data.frame(moduleColor=moduleColors.DS2)

# Annotate the modules with GO functional enrichment
entrezIDs = as.character(geneInfo$geneEntrezID)
GOenr = GOenrichmentAnalysis(geneInfo$moduleColor, allLLIDs, organism = "human", nBestP = 10);
GOenr.terms.all = GOenr$bestPTerms[[4]]$enrichment