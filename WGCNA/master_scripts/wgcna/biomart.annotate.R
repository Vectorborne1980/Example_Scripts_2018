
###################################################################################################
# ANNOTATE ENSEMBL GENES
###################################################################################################

biomart.annotate <- function (IDs, filter, attributes) {
  
  # Query biomaRt for annotation information keyed to IDs (usually Ensembl)
  library('biomaRt')
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  # View available marts, inputs and outputs
  #datasets <- listDatasets(mart)
  #filters <- listFilters(mart)
  #attribs <- listAttributes(mart)
  #ensemblID <- colnames(datExpr)
  
  # Give this a minute or so. Also, remember not all Ensembl areas are actually transcribed and thus assigned an Entrez ID
  print(paste("Querying biomaRt...",sep=""))
  biomart.output <- getBM(filters=filter, attributes=attributes, values=IDs, mart=mart)
  print(paste("Query complete.",sep=""))
  
  return(biomart.output)
}
