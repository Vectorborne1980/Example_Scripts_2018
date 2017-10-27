# A way to get the names of all Gene Ontology terms for Biological Processes at 
# level X as well as the genes (human gene symbols) annotated to each of the level X GO terms.
library(GO.db)
setwd("C:/Users/Nick Geraci/Desktop/UVA/Gene_Ontology/")

getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

# Start with the very top biological process ID (GO:0008150)
# https://gowiki.tamu.edu/wiki/index.php/Category:GO:0008150_!_biological_process
# The immune system process ID is GO:0002376, which is one level below GO:0008150
# The gene ontology immune process terms are curated through http://wiki.geneontology.org/index.php/Immunology
# The gene ontology neurobiology terms are curated through http://wiki.geneontology.org/index.php/Neurobiology_Project

level01_BP_terms <- getAllBPChildren("GO:0008150")      # 27 terms
level02_BP_terms <- getAllBPChildren(level01_BP_terms)  # 399 terms
level03_BP_terms <- getAllBPChildren(level02_BP_terms)  # 3658 terms
level04_BP_terms <- getAllBPChildren(level03_BP_terms)  # 10792 terms
level05_BP_terms <- getAllBPChildren(level04_BP_terms)  # 18101 terms
level06_BP_terms <- getAllBPChildren(level05_BP_terms)  # 20756 terms
level07_BP_terms <- getAllBPChildren(level06_BP_terms)  # 20342 terms
level08_BP_terms <- getAllBPChildren(level07_BP_terms)  # 17868 terms
level09_BP_terms <- getAllBPChildren(level08_BP_terms)  # 14358 terms
level10_BP_terms <- getAllBPChildren(level09_BP_terms)  # 10322 terms
level11_BP_terms <- getAllBPChildren(level10_BP_terms)  # 6828 terms
level12_BP_terms <- getAllBPChildren(level11_BP_terms)  # 4263 terms
level13_BP_terms <- getAllBPChildren(level12_BP_terms)  # 2334 terms
level14_BP_terms <- getAllBPChildren(level13_BP_terms)  # 1136 terms
level15_BP_terms <- getAllBPChildren(level14_BP_terms)  # 500 terms
level16_BP_terms <- getAllBPChildren(level15_BP_terms)  # 192 terms
level17_BP_terms <- getAllBPChildren(level16_BP_terms)  # 72 terms
level18_BP_terms <- getAllBPChildren(level17_BP_terms)  # 21 terms
level19_BP_terms <- getAllBPChildren(level18_BP_terms)  # 2 terms
level20_BP_terms <- getAllBPChildren(level19_BP_terms)  # 0 terms

### HUMAN
# Get the gene Entrez IDs mapped to all of the terms
library(org.Hs.eg.db)
level01_genes <- mget(intersect(level01_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG) # 27 terms

# There is not really a unique "level" of GO terms. Many of your level 5 terms can also be level
# 4 terms, or level 3 terms or level 6 terms, etc. This is due to the acyclic nature of the GO
# terms and the multiple paths possible from one ancestor to one descendent.
# To illustrate this:

length(intersect(level04_BP_terms, level05_BP_terms)) # 8933

# Which means 8933 terms belong to level 4 and 5.
# What many think of as unique is the "minimum level" of a term (i.e., the length of the shortest path 
# between the term and the root of the ontology). In other words, the "minimum level" of a term 
# is its distance to the root.

# NOTE: I, personally, think that levels should be measured from the terminus to a parent.

# If you want the BP terms that are at distance 2 from the root, just do: 

dist02_BP_terms <- setdiff(level02_BP_terms,c(level01_BP_terms))
length(dist02_BP_terms) # 395
level02_genes <- mget(intersect(dist02_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)

# Continue to all levels

#03
dist03_BP_terms <- setdiff(level03_BP_terms,c(level02_BP_terms, level01_BP_terms))
length(dist03_BP_terms) # 3404
level03_genes <- mget(intersect(dist03_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#04
dist04_BP_terms <- setdiff(level04_BP_terms,c(level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist04_BP_terms) # 7987
level04_genes <- mget(intersect(dist04_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#05
dist05_BP_terms <- setdiff(level05_BP_terms,c(level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist05_BP_terms) # 8831
level05_genes <- mget(intersect(dist05_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#06
dist06_BP_terms <- setdiff(level06_BP_terms,c(level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist06_BP_terms) # 5100
level06_genes <- mget(intersect(dist06_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#07
dist07_BP_terms <- setdiff(level07_BP_terms,c(level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist07_BP_terms) # 2608
level07_genes <- mget(intersect(dist07_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#08
dist08_BP_terms <- setdiff(level08_BP_terms,c(level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist08_BP_terms) # 988
level08_genes <- mget(intersect(dist08_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#09
dist09_BP_terms <- setdiff(level09_BP_terms,c(level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist09_BP_terms) # 197
level09_genes <- mget(intersect(dist09_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
#10
dist10_BP_terms <- setdiff(level10_BP_terms,c(level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist10_BP_terms) # 39
level10_genes <- mget(intersect(dist10_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
# 11
dist11_BP_terms <- setdiff(level11_BP_terms,c(level10_BP_terms, level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist11_BP_terms) # 9
level11_genes <- mget(intersect(dist11_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)
# 12
dist12_BP_terms <- setdiff(level11_BP_terms,c(level11_BP_terms, level10_BP_terms, level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist12_BP_terms) # 0
level12_genes <- mget(intersect(dist12_BP_terms,keys(org.Hs.egGO2EG)),org.Hs.egGO2EG)

rm(list = ls(pattern = glob2rx("*BP_terms")))
rm(level12_genes)

# Refromat the list objects into a single data frame with unique GO term and Entrez ID pairings
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
# Get the table of GO IDs and their associated terms
GOid2Terms <- read.delim("Gene_Ontology_Biological_Process_GO-IDs-to-Terms_NSG_28Sept2017.txt", header = TRUE, sep = "\t")
options(warn = -1) # Upon joining the GO IDs to their terms, there are warnings that are better left suppressed
# 01
max.length <- max(sapply(level01_genes, length)) # Get the maximum number of genes mapped to any one GO term
level01_genes <- lapply(level01_genes, function(v) { c(v, rep(NA, max.length-length(v)))}) # Add NA values to list elements
level01_genes <- do.call(rbind, level01_genes) # Convert lists to a matrix through rbind
level01_genes <- as.data.frame(level01_genes) # Force class conversion to data frame
colnames(level01_genes) <- c(1:ncol(level01_genes)) # Provide generic column headers
level01_genes <- t(level01_genes) # Transpose data frame
level01_genes <- stack(level01_genes) # Stack all values with column names into a single data frame
colnames(level01_genes) <- c("row","GO_ID","count","Entrez_ID") # Re-name column headers
level01_genes <- level01_genes[!is.na(level01_genes$Entrez_ID),] # If the row contains 'NA' in the 'Entrez_ID' column, remove it
level01_genes <- level01_genes[,-1] # Remove worthless column
level01_genes <- level01_genes[,-2] # Remove another worthless column
level01_genes <- unique(level01_genes) # Filter only unique GO term and gene Entrez ID mappings --> count: 249
symbols <- getSYMBOL(c(level01_genes$Entrez_ID),data='org.Hs.eg') # Get a vector of gene symbols for each Entrez ID
level01_genes[,"Gene_Symbols"] <- as.character(symbols) # Tack on gene symbols to the data frame
level01_genes <- as.data.frame(level01_genes) # Change to a formal data frame class object
index <- which(GOid2Terms$GO_ID%in%level01_genes$GO_ID) # Select just the GO IDs needed for the terms in the list
selectedGO <- as.data.frame(GOid2Terms[index,]) # Finish selection
index <- which(level01_genes$GO_ID%in%selectedGO$GO_ID) # As GO.db is continually updating there may be more IDs than my most recent gene to term list
level01_genes <- level01_genes[index,] # It is necessary to eliminate new terms until my list is updated on a semi-regular basis (should be quarterly at least)
identical(nrow(selectedGO),length(unique(level01_genes$GO_ID))) # QC check: are all of the GO IDs available?
level01_genes <- left_join(level01_genes,GOid2Terms) # Index and match the GO terms by the GO IDs
level01_genes <- level01_genes[,c(1,4,2,3,5)] # Rearrange the data frame
# 02
max.length <- max(sapply(level02_genes, length))
level02_genes <- lapply(level02_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level02_genes <- do.call(rbind, level02_genes)
level02_genes <- as.data.frame(level02_genes)
colnames(level02_genes) <- c(1:ncol(level02_genes))
level02_genes <- t(level02_genes)
level02_genes <- stack(level02_genes)
colnames(level02_genes) <- c("row","GO_ID","count","Entrez_ID")
level02_genes <- level02_genes[!is.na(level02_genes$Entrez_ID),]
level02_genes <- level02_genes[,-1]
level02_genes <- level02_genes[,-2]
level02_genes <- unique(level02_genes) # count: 4542
symbols <- getSYMBOL(c(level02_genes$Entrez_ID),data='org.Hs.eg')
level02_genes[,"Gene_Symbols"] <- as.character(symbols)
level02_genes <- as.data.frame(level02_genes)
index <- which(GOid2Terms$GO_ID%in%level02_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level02_genes$GO_ID%in%selectedGO$GO_ID)
level02_genes <- level02_genes[index,]
identical(nrow(selectedGO),length(unique(level02_genes$GO_ID)))
level02_genes <- left_join(level02_genes,GOid2Terms)
level02_genes <- level02_genes[,c(1,4,2,3,5)]
# 03
max.length <- max(sapply(level03_genes, length))
level03_genes <- lapply(level03_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level03_genes <- do.call(rbind, level03_genes)
level03_genes <- as.data.frame(level03_genes)
colnames(level03_genes) <- c(1:ncol(level03_genes))
level03_genes <- t(level03_genes)
level03_genes <- stack(level03_genes)
colnames(level03_genes) <- c("row","GO_ID","count","Entrez_ID")
level03_genes <- level03_genes[!is.na(level03_genes$Entrez_ID),]
level03_genes <- level03_genes[,-1]
level03_genes <- level03_genes[,-2]
level03_genes <- unique(level03_genes) # count: 21678
symbols <- getSYMBOL(c(level03_genes$Entrez_ID),data='org.Hs.eg')
level03_genes[,"Gene_Symbols"] <- as.character(symbols)
level03_genes <- as.data.frame(level03_genes)
index <- which(GOid2Terms$GO_ID%in%level03_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level03_genes$GO_ID%in%selectedGO$GO_ID)
level03_genes <- level03_genes[index,]
identical(nrow(selectedGO),length(unique(level03_genes$GO_ID)))
level03_genes <- left_join(level03_genes,GOid2Terms)
level03_genes <- level03_genes[,c(1,4,2,3,5)]
# 04
max.length <- max(sapply(level04_genes, length))
level04_genes <- lapply(level04_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level04_genes <- do.call(rbind, level04_genes)
level04_genes <- as.data.frame(level04_genes)
colnames(level04_genes) <- c(1:ncol(level04_genes))
level04_genes <- t(level04_genes)
level04_genes <- stack(level04_genes)
colnames(level04_genes) <- c("row","GO_ID","count","Entrez_ID")
level04_genes <- level04_genes[!is.na(level04_genes$Entrez_ID),]
level04_genes <- level04_genes[,-1]
level04_genes <- level04_genes[,-2]
level04_genes <- unique(level04_genes) # count: 33134
symbols <- getSYMBOL(c(level04_genes$Entrez_ID),data='org.Hs.eg')
level04_genes[,"Gene_Symbols"] <- as.character(symbols)
level04_genes <- as.data.frame(level04_genes)
index <- which(GOid2Terms$GO_ID%in%level04_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level04_genes$GO_ID%in%selectedGO$GO_ID)
level04_genes <- level04_genes[index,]
identical(nrow(selectedGO),length(unique(level04_genes$GO_ID)))
level04_genes <- left_join(level04_genes,GOid2Terms)
level04_genes <- level04_genes[,c(1,4,2,3,5)]
# 05
max.length <- max(sapply(level05_genes, length))
level05_genes <- lapply(level05_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level05_genes <- do.call(rbind, level05_genes)
level05_genes <- as.data.frame(level05_genes)
colnames(level05_genes) <- c(1:ncol(level05_genes))
level05_genes <- t(level05_genes)
level05_genes <- stack(level05_genes)
colnames(level05_genes) <- c("row","GO_ID","count","Entrez_ID")
level05_genes <- level05_genes[!is.na(level05_genes$Entrez_ID),]
level05_genes <- level05_genes[,-1]
level05_genes <- level05_genes[,-2]
level05_genes <- unique(level05_genes) # count: 33482
symbols <- getSYMBOL(c(level05_genes$Entrez_ID),data='org.Hs.eg')
level05_genes[,"Gene_Symbols"] <- as.character(symbols)
level05_genes <- as.data.frame(level05_genes)
index <- which(GOid2Terms$GO_ID%in%level05_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level05_genes$GO_ID%in%selectedGO$GO_ID)
level05_genes <- level05_genes[index,]
identical(nrow(selectedGO),length(unique(level05_genes$GO_ID)))
level05_genes <- left_join(level05_genes,GOid2Terms)
level05_genes <- level05_genes[,c(1,4,2,3,5)]
# 06
max.length <- max(sapply(level06_genes, length))
level06_genes <- lapply(level06_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level06_genes <- do.call(rbind, level06_genes)
level06_genes <- as.data.frame(level06_genes)
colnames(level06_genes) <- c(1:ncol(level06_genes))
level06_genes <- t(level06_genes)
level06_genes <- stack(level06_genes)
colnames(level06_genes) <- c("row","GO_ID","count","Entrez_ID")
level06_genes <- level06_genes[!is.na(level06_genes$Entrez_ID),]
level06_genes <- level06_genes[,-1]
level06_genes <- level06_genes[,-2]
level06_genes <- unique(level06_genes) # count: 19445
symbols <- getSYMBOL(c(level06_genes$Entrez_ID),data='org.Hs.eg')
level06_genes[,"Gene_Symbols"] <- as.character(symbols)
level06_genes <- as.data.frame(level06_genes)
index <- which(GOid2Terms$GO_ID%in%level06_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level06_genes$GO_ID%in%selectedGO$GO_ID)
level06_genes <- level06_genes[index,]
identical(nrow(selectedGO),length(unique(level06_genes$GO_ID)))
level06_genes <- left_join(level06_genes,GOid2Terms)
level06_genes <- level06_genes[,c(1,4,2,3,5)]
# 07
max.length <- max(sapply(level07_genes, length))
level07_genes <- lapply(level07_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level07_genes <- do.call(rbind, level07_genes)
level07_genes <- as.data.frame(level07_genes)
colnames(level07_genes) <- c(1:ncol(level07_genes))
level07_genes <- t(level07_genes)
level07_genes <- stack(level07_genes)
colnames(level07_genes) <- c("row","GO_ID","count","Entrez_ID")
level07_genes <- level07_genes[!is.na(level07_genes$Entrez_ID),]
level07_genes <- level07_genes[,-1]
level07_genes <- level07_genes[,-2]
level07_genes <- unique(level07_genes) # count: 9040
symbols <- getSYMBOL(c(level07_genes$Entrez_ID),data='org.Hs.eg')
level07_genes[,"Gene_Symbols"] <- as.character(symbols)
level07_genes <- as.data.frame(level07_genes)
index <- which(GOid2Terms$GO_ID%in%level07_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level07_genes$GO_ID%in%selectedGO$GO_ID)
level07_genes <- level07_genes[index,]
identical(nrow(selectedGO),length(unique(level07_genes$GO_ID)))
level07_genes <- left_join(level07_genes,GOid2Terms)
level07_genes <- level07_genes[,c(1,4,2,3,5)]
# 08
max.length <- max(sapply(level08_genes, length))
level08_genes <- lapply(level08_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level08_genes <- do.call(rbind, level08_genes)
level08_genes <- as.data.frame(level08_genes)
colnames(level08_genes) <- c(1:ncol(level08_genes))
level08_genes <- t(level08_genes)
level08_genes <- stack(level08_genes)
colnames(level08_genes) <- c("row","GO_ID","count","Entrez_ID")
level08_genes <- level08_genes[!is.na(level08_genes$Entrez_ID),]
level08_genes <- level08_genes[,-1]
level08_genes <- level08_genes[,-2]
level08_genes <- unique(level08_genes) # count: 2919
symbols <- getSYMBOL(c(level08_genes$Entrez_ID),data='org.Hs.eg')
level08_genes[,"Gene_Symbols"] <- as.character(symbols)
level08_genes <- as.data.frame(level08_genes)
index <- which(GOid2Terms$GO_ID%in%level08_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level08_genes$GO_ID%in%selectedGO$GO_ID)
level08_genes <- level08_genes[index,]
identical(nrow(selectedGO),length(unique(level08_genes$GO_ID)))
level08_genes <- left_join(level08_genes,GOid2Terms)
level08_genes <- level08_genes[,c(1,4,2,3,5)]
# 09
max.length <- max(sapply(level09_genes, length))
level09_genes <- lapply(level09_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level09_genes <- do.call(rbind, level09_genes)
level09_genes <- as.data.frame(level09_genes)
colnames(level09_genes) <- c(1:ncol(level09_genes))
level09_genes <- t(level09_genes)
level09_genes <- stack(level09_genes)
colnames(level09_genes) <- c("row","GO_ID","count","Entrez_ID")
level09_genes <- level09_genes[!is.na(level09_genes$Entrez_ID),]
level09_genes <- level09_genes[,-1]
level09_genes <- level09_genes[,-2]
level09_genes <- unique(level09_genes) # count: 824
symbols <- getSYMBOL(c(level09_genes$Entrez_ID),data='org.Hs.eg')
level09_genes[,"Gene_Symbols"] <- as.character(symbols)
level09_genes <- as.data.frame(level09_genes)
index <- which(GOid2Terms$GO_ID%in%level09_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level09_genes$GO_ID%in%selectedGO$GO_ID)
level09_genes <- level09_genes[index,]
identical(nrow(selectedGO),length(unique(level09_genes$GO_ID)))
level09_genes <- left_join(level09_genes,GOid2Terms)
level09_genes <- level09_genes[,c(1,4,2,3,5)]
# 10
max.length <- max(sapply(level10_genes, length))
level10_genes <- lapply(level10_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level10_genes <- do.call(rbind, level10_genes)
level10_genes <- as.data.frame(level10_genes)
colnames(level10_genes) <- c(1:ncol(level10_genes))
level10_genes <- t(level10_genes)
level10_genes <- stack(level10_genes)
colnames(level10_genes) <- c("row","GO_ID","count","Entrez_ID")
level10_genes <- level10_genes[!is.na(level10_genes$Entrez_ID),]
level10_genes <- level10_genes[,-1]
level10_genes <- level10_genes[,-2]
level10_genes <- unique(level10_genes) # count: 330
symbols <- getSYMBOL(c(level10_genes$Entrez_ID),data='org.Hs.eg')
level10_genes[,"Gene_Symbols"] <- as.character(symbols)
level10_genes <- as.data.frame(level10_genes)
index <- which(GOid2Terms$GO_ID%in%level10_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level10_genes$GO_ID%in%selectedGO$GO_ID)
level10_genes <- level10_genes[index,]
identical(nrow(selectedGO),length(unique(level10_genes$GO_ID)))
level10_genes <- left_join(level10_genes,GOid2Terms)
level10_genes <- level10_genes[,c(1,4,2,3,5)]
# 11
max.length <- max(sapply(level11_genes, length))
level11_genes <- lapply(level11_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level11_genes <- do.call(rbind, level11_genes)
level11_genes <- as.data.frame(level11_genes)
colnames(level11_genes) <- c(1:ncol(level11_genes))
level11_genes <- t(level11_genes)
level11_genes <- stack(level11_genes)
colnames(level11_genes) <- c("row","GO_ID","count","Entrez_ID")
level11_genes <- level11_genes[!is.na(level11_genes$Entrez_ID),]
level11_genes <- level11_genes[,-1]
level11_genes <- level11_genes[,-2]
level11_genes <- unique(level11_genes) # count: 13
symbols <- getSYMBOL(c(level11_genes$Entrez_ID),data='org.Hs.eg')
level11_genes[,"Gene_Symbols"] <- as.character(symbols)
level11_genes <- as.data.frame(level11_genes)
index <- which(GOid2Terms$GO_ID%in%level11_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level11_genes$GO_ID%in%selectedGO$GO_ID)
level11_genes <- level11_genes[index,]
identical(nrow(selectedGO),length(unique(level11_genes$GO_ID)))
level11_genes <- left_join(level11_genes,GOid2Terms)
level11_genes <- level11_genes[,c(1,4,2,3,5)]

rm(index,selectedGO,GOid2Terms,symbols,max.length)

options(warn = 0) # Turn warnings back on

level01_genes$GO_Level <- c(rep(c(1),times = NROW(level01_genes)))
level02_genes$GO_Level <- c(rep(c(2),times = NROW(level02_genes)))
level03_genes$GO_Level <- c(rep(c(3),times = NROW(level03_genes)))
level04_genes$GO_Level <- c(rep(c(4),times = NROW(level04_genes)))
level05_genes$GO_Level <- c(rep(c(5),times = NROW(level05_genes)))
level06_genes$GO_Level <- c(rep(c(6),times = NROW(level06_genes)))
level07_genes$GO_Level <- c(rep(c(7),times = NROW(level07_genes)))
level08_genes$GO_Level <- c(rep(c(8),times = NROW(level08_genes)))
level09_genes$GO_Level <- c(rep(c(9),times = NROW(level09_genes)))
level10_genes$GO_Level <- c(rep(c(10),times = NROW(level10_genes)))
level11_genes$GO_Level <- c(rep(c(11),times = NROW(level11_genes)))

GO.level.genes.human <- rbind(level01_genes,level02_genes,level03_genes,level04_genes,
                              level05_genes,level06_genes,level07_genes,level08_genes,
                              level09_genes,level10_genes,level11_genes)

save(GO.level.genes.human, file="Unique_Human_GO-Entrez_Pairs_NSG_28Sept2017.RData")

rm(list = ls(pattern = glob2rx("*_genes")))

### MOUSE
# Get the gene Entrez IDs mapped to all of the terms
library(org.Mm.eg.db)
level01_genes <- mget(intersect(level01_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG) # 27 terms

# There is not really a unique "level" of GO terms. Many of your level 5 terms can also be level
# 4 terms, or level 3 terms or level 6 terms, etc. This is due to the acyclic nature of the GO
# terms and the multiple paths possible from one ancestor to one descendent.
# To illustrate this:

length(intersect(level04_BP_terms, level05_BP_terms)) # 8933

# Which means 8933 terms belong to level 4 and 5.
# What many think of as unique is the "minimum level" of a term (i.e., the length of the shortest path 
# between the term and the root of the ontology). In other words, the "minimum level" of a term 
# is its distance to the root.

# NOTE: I, personally, think that levels should be measured from the terminus to a parent.

# If you want the BP terms that are at distance 2 from the root, just do: 

dist02_BP_terms <- setdiff(level02_BP_terms,c(level01_BP_terms))
length(dist02_BP_terms) # 395
level02_genes <- mget(intersect(dist02_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)

# Continue to all levels

#03
dist03_BP_terms <- setdiff(level03_BP_terms,c(level02_BP_terms, level01_BP_terms))
length(dist03_BP_terms) # 3404
level03_genes <- mget(intersect(dist03_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#04
dist04_BP_terms <- setdiff(level04_BP_terms,c(level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist04_BP_terms) # 7987
level04_genes <- mget(intersect(dist04_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#05
dist05_BP_terms <- setdiff(level05_BP_terms,c(level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist05_BP_terms) # 8831
level05_genes <- mget(intersect(dist05_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#06
dist06_BP_terms <- setdiff(level06_BP_terms,c(level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist06_BP_terms) # 5100
level06_genes <- mget(intersect(dist06_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#07
dist07_BP_terms <- setdiff(level07_BP_terms,c(level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist07_BP_terms) # 2608
level07_genes <- mget(intersect(dist07_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#08
dist08_BP_terms <- setdiff(level08_BP_terms,c(level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist08_BP_terms) # 988
level08_genes <- mget(intersect(dist08_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#09
dist09_BP_terms <- setdiff(level09_BP_terms,c(level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist09_BP_terms) # 197
level09_genes <- mget(intersect(dist09_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
#10
dist10_BP_terms <- setdiff(level10_BP_terms,c(level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist10_BP_terms) # 39
level10_genes <- mget(intersect(dist10_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
# 11
dist11_BP_terms <- setdiff(level11_BP_terms,c(level10_BP_terms, level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist11_BP_terms) # 9
level11_genes <- mget(intersect(dist11_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)
# 12
dist12_BP_terms <- setdiff(level11_BP_terms,c(level11_BP_terms, level10_BP_terms, level09_BP_terms
                                              , level08_BP_terms, level07_BP_terms, level06_BP_terms
                                              , level05_BP_terms, level04_BP_terms, level03_BP_terms
                                              , level02_BP_terms, level01_BP_terms))
length(dist12_BP_terms) # 0
level12_genes <- mget(intersect(dist12_BP_terms,keys(org.Mm.egGO2EG)),org.Mm.egGO2EG)

rm(list = ls(pattern = glob2rx("*BP_terms")))
rm(level12_genes)

# Refromat the list objects into a single data frame with unique GO term and Entrez ID pairings
library(org.Mm.eg.db)
library(annotate)
library(dplyr)
# Get the table of GO IDs and their associated terms
GOid2Terms <- read.delim("Gene_Ontology_Biological_Process_GO-IDs-to-Terms_NSG_28Sept2017.txt", header = TRUE, sep = "\t")
options(warn = -1) # Upon joining the GO IDs to their terms, there are warnings that are better left suppressed
# 01
max.length <- max(sapply(level01_genes, length)) # Get the maximum number of genes mapped to any one GO term
level01_genes <- lapply(level01_genes, function(v) { c(v, rep(NA, max.length-length(v)))}) # Add NA values to list elements
level01_genes <- do.call(rbind, level01_genes) # Convert lists to a matrix through rbind
level01_genes <- as.data.frame(level01_genes) # Force class conversion to data frame
colnames(level01_genes) <- c(1:ncol(level01_genes)) # Provide generic column headers
level01_genes <- t(level01_genes) # Transpose data frame
level01_genes <- stack(level01_genes) # Stack all values with column names into a single data frame
colnames(level01_genes) <- c("row","GO_ID","count","Entrez_ID") # Re-name column headers
level01_genes <- level01_genes[!is.na(level01_genes$Entrez_ID),] # If the row contains 'NA' in the 'Entrez_ID' column, remove it
level01_genes <- level01_genes[,-1] # Remove worthless column
level01_genes <- level01_genes[,-2] # Remove another worthless column
level01_genes <- unique(level01_genes) # Filter only unique GO term and gene Entrez ID mappings --> count: 1250
symbols <- getSYMBOL(c(level01_genes$Entrez_ID),data='org.Mm.eg') # Get a vector of gene symbols for each Entrez ID
level01_genes[,"Gene_Symbols"] <- as.character(symbols) # Tack on gene symbols to the data frame
level01_genes <- as.data.frame(level01_genes) # Change to a formal data frame class object
index <- which(GOid2Terms$GO_ID%in%level01_genes$GO_ID) # Select just the GO IDs needed for the terms in the list
selectedGO <- as.data.frame(GOid2Terms[index,]) # Finish selection
index <- which(level01_genes$GO_ID%in%selectedGO$GO_ID) # As GO.db is continually updating there may be more IDs than my most recent gene to term list
level01_genes <- level01_genes[index,] # It is necessary to eliminate new terms until my list is updated on a semi-regular basis (should be quarterly at least)
identical(nrow(selectedGO),length(unique(level01_genes$GO_ID))) # QC check: are all of the GO IDs available?
level01_genes <- left_join(level01_genes,GOid2Terms) # Index and match the GO terms by the GO IDs
level01_genes <- level01_genes[,c(1,4,2,3,5)] # Rearrange the data frame
# 02
max.length <- max(sapply(level02_genes, length))
level02_genes <- lapply(level02_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level02_genes <- do.call(rbind, level02_genes)
level02_genes <- as.data.frame(level02_genes)
colnames(level02_genes) <- c(1:ncol(level02_genes))
level02_genes <- t(level02_genes)
level02_genes <- stack(level02_genes)
colnames(level02_genes) <- c("row","GO_ID","count","Entrez_ID")
level02_genes <- level02_genes[!is.na(level02_genes$Entrez_ID),]
level02_genes <- level02_genes[,-1]
level02_genes <- level02_genes[,-2]
level02_genes <- unique(level02_genes) # count: 4451
symbols <- getSYMBOL(c(level02_genes$Entrez_ID),data='org.Mm.eg')
level02_genes[,"Gene_Symbols"] <- as.character(symbols)
level02_genes <- as.data.frame(level02_genes)
index <- which(GOid2Terms$GO_ID%in%level02_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level02_genes$GO_ID%in%selectedGO$GO_ID)
level02_genes <- level02_genes[index,]
identical(nrow(selectedGO),length(unique(level02_genes$GO_ID)))
level02_genes <- left_join(level02_genes,GOid2Terms)
level02_genes <- level02_genes[,c(1,4,2,3,5)]
# 03
max.length <- max(sapply(level03_genes, length))
level03_genes <- lapply(level03_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level03_genes <- do.call(rbind, level03_genes)
level03_genes <- as.data.frame(level03_genes)
colnames(level03_genes) <- c(1:ncol(level03_genes))
level03_genes <- t(level03_genes)
level03_genes <- stack(level03_genes)
colnames(level03_genes) <- c("row","GO_ID","count","Entrez_ID")
level03_genes <- level03_genes[!is.na(level03_genes$Entrez_ID),]
level03_genes <- level03_genes[,-1]
level03_genes <- level03_genes[,-2]
level03_genes <- unique(level03_genes) # count: 24719
symbols <- getSYMBOL(c(level03_genes$Entrez_ID),data='org.Mm.eg')
level03_genes[,"Gene_Symbols"] <- as.character(symbols)
level03_genes <- as.data.frame(level03_genes)
index <- which(GOid2Terms$GO_ID%in%level03_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level03_genes$GO_ID%in%selectedGO$GO_ID)
level03_genes <- level03_genes[index,]
identical(nrow(selectedGO),length(unique(level03_genes$GO_ID)))
level03_genes <- left_join(level03_genes,GOid2Terms)
level03_genes <- level03_genes[,c(1,4,2,3,5)]
# 04
max.length <- max(sapply(level04_genes, length))
level04_genes <- lapply(level04_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level04_genes <- do.call(rbind, level04_genes)
level04_genes <- as.data.frame(level04_genes)
colnames(level04_genes) <- c(1:ncol(level04_genes))
level04_genes <- t(level04_genes)
level04_genes <- stack(level04_genes)
colnames(level04_genes) <- c("row","GO_ID","count","Entrez_ID")
level04_genes <- level04_genes[!is.na(level04_genes$Entrez_ID),]
level04_genes <- level04_genes[,-1]
level04_genes <- level04_genes[,-2]
level04_genes <- unique(level04_genes) # count: 34780
symbols <- getSYMBOL(c(level04_genes$Entrez_ID),data='org.Mm.eg')
level04_genes[,"Gene_Symbols"] <- as.character(symbols)
level04_genes <- as.data.frame(level04_genes)
index <- which(GOid2Terms$GO_ID%in%level04_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level04_genes$GO_ID%in%selectedGO$GO_ID)
level04_genes <- level04_genes[index,]
identical(nrow(selectedGO),length(unique(level04_genes$GO_ID)))
level04_genes <- left_join(level04_genes,GOid2Terms)
level04_genes <- level04_genes[,c(1,4,2,3,5)]
# 05
max.length <- max(sapply(level05_genes, length))
level05_genes <- lapply(level05_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level05_genes <- do.call(rbind, level05_genes)
level05_genes <- as.data.frame(level05_genes)
colnames(level05_genes) <- c(1:ncol(level05_genes))
level05_genes <- t(level05_genes)
level05_genes <- stack(level05_genes)
colnames(level05_genes) <- c("row","GO_ID","count","Entrez_ID")
level05_genes <- level05_genes[!is.na(level05_genes$Entrez_ID),]
level05_genes <- level05_genes[,-1]
level05_genes <- level05_genes[,-2]
level05_genes <- unique(level05_genes) # count: 33615
symbols <- getSYMBOL(c(level05_genes$Entrez_ID),data='org.Mm.eg')
level05_genes[,"Gene_Symbols"] <- as.character(symbols)
level05_genes <- as.data.frame(level05_genes)
index <- which(GOid2Terms$GO_ID%in%level05_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level05_genes$GO_ID%in%selectedGO$GO_ID)
level05_genes <- level05_genes[index,]
identical(nrow(selectedGO),length(unique(level05_genes$GO_ID)))
level05_genes <- left_join(level05_genes,GOid2Terms)
level05_genes <- level05_genes[,c(1,4,2,3,5)]
# 06
max.length <- max(sapply(level06_genes, length))
level06_genes <- lapply(level06_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level06_genes <- do.call(rbind, level06_genes)
level06_genes <- as.data.frame(level06_genes)
colnames(level06_genes) <- c(1:ncol(level06_genes))
level06_genes <- t(level06_genes)
level06_genes <- stack(level06_genes)
colnames(level06_genes) <- c("row","GO_ID","count","Entrez_ID")
level06_genes <- level06_genes[!is.na(level06_genes$Entrez_ID),]
level06_genes <- level06_genes[,-1]
level06_genes <- level06_genes[,-2]
level06_genes <- unique(level06_genes) # count: 17758
symbols <- getSYMBOL(c(level06_genes$Entrez_ID),data='org.Mm.eg')
level06_genes[,"Gene_Symbols"] <- as.character(symbols)
level06_genes <- as.data.frame(level06_genes)
index <- which(GOid2Terms$GO_ID%in%level06_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level06_genes$GO_ID%in%selectedGO$GO_ID)
level06_genes <- level06_genes[index,]
identical(nrow(selectedGO),length(unique(level06_genes$GO_ID)))
level06_genes <- left_join(level06_genes,GOid2Terms)
level06_genes <- level06_genes[,c(1,4,2,3,5)]
# 07
max.length <- max(sapply(level07_genes, length))
level07_genes <- lapply(level07_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level07_genes <- do.call(rbind, level07_genes)
level07_genes <- as.data.frame(level07_genes)
colnames(level07_genes) <- c(1:ncol(level07_genes))
level07_genes <- t(level07_genes)
level07_genes <- stack(level07_genes)
colnames(level07_genes) <- c("row","GO_ID","count","Entrez_ID")
level07_genes <- level07_genes[!is.na(level07_genes$Entrez_ID),]
level07_genes <- level07_genes[,-1]
level07_genes <- level07_genes[,-2]
level07_genes <- unique(level07_genes) # count: 7222
symbols <- getSYMBOL(c(level07_genes$Entrez_ID),data='org.Mm.eg')
level07_genes[,"Gene_Symbols"] <- as.character(symbols)
level07_genes <- as.data.frame(level07_genes)
index <- which(GOid2Terms$GO_ID%in%level07_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level07_genes$GO_ID%in%selectedGO$GO_ID)
level07_genes <- level07_genes[index,]
identical(nrow(selectedGO),length(unique(level07_genes$GO_ID)))
level07_genes <- left_join(level07_genes,GOid2Terms)
level07_genes <- level07_genes[,c(1,4,2,3,5)]
# 08
max.length <- max(sapply(level08_genes, length))
level08_genes <- lapply(level08_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level08_genes <- do.call(rbind, level08_genes)
level08_genes <- as.data.frame(level08_genes)
colnames(level08_genes) <- c(1:ncol(level08_genes))
level08_genes <- t(level08_genes)
level08_genes <- stack(level08_genes)
colnames(level08_genes) <- c("row","GO_ID","count","Entrez_ID")
level08_genes <- level08_genes[!is.na(level08_genes$Entrez_ID),]
level08_genes <- level08_genes[,-1]
level08_genes <- level08_genes[,-2]
level08_genes <- unique(level08_genes) # count: 2480
symbols <- getSYMBOL(c(level08_genes$Entrez_ID),data='org.Mm.eg')
level08_genes[,"Gene_Symbols"] <- as.character(symbols)
level08_genes <- as.data.frame(level08_genes)
index <- which(GOid2Terms$GO_ID%in%level08_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level08_genes$GO_ID%in%selectedGO$GO_ID)
level08_genes <- level08_genes[index,]
identical(nrow(selectedGO),length(unique(level08_genes$GO_ID)))
level08_genes <- left_join(level08_genes,GOid2Terms)
level08_genes <- level08_genes[,c(1,4,2,3,5)]
# 09
max.length <- max(sapply(level09_genes, length))
level09_genes <- lapply(level09_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level09_genes <- do.call(rbind, level09_genes)
level09_genes <- as.data.frame(level09_genes)
colnames(level09_genes) <- c(1:ncol(level09_genes))
level09_genes <- t(level09_genes)
level09_genes <- stack(level09_genes)
colnames(level09_genes) <- c("row","GO_ID","count","Entrez_ID")
level09_genes <- level09_genes[!is.na(level09_genes$Entrez_ID),]
level09_genes <- level09_genes[,-1]
level09_genes <- level09_genes[,-2]
level09_genes <- unique(level09_genes) # count: 707
symbols <- getSYMBOL(c(level09_genes$Entrez_ID),data='org.Mm.eg')
level09_genes[,"Gene_Symbols"] <- as.character(symbols)
level09_genes <- as.data.frame(level09_genes)
index <- which(GOid2Terms$GO_ID%in%level09_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level09_genes$GO_ID%in%selectedGO$GO_ID)
level09_genes <- level09_genes[index,]
identical(nrow(selectedGO),length(unique(level09_genes$GO_ID)))
level09_genes <- left_join(level09_genes,GOid2Terms)
level09_genes <- level09_genes[,c(1,4,2,3,5)]
# 10
max.length <- max(sapply(level10_genes, length))
level10_genes <- lapply(level10_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level10_genes <- do.call(rbind, level10_genes)
level10_genes <- as.data.frame(level10_genes)
colnames(level10_genes) <- c(1:ncol(level10_genes))
level10_genes <- t(level10_genes)
level10_genes <- stack(level10_genes)
colnames(level10_genes) <- c("row","GO_ID","count","Entrez_ID")
level10_genes <- level10_genes[!is.na(level10_genes$Entrez_ID),]
level10_genes <- level10_genes[,-1]
level10_genes <- level10_genes[,-2]
level10_genes <- unique(level10_genes) # count: 351
symbols <- getSYMBOL(c(level10_genes$Entrez_ID),data='org.Mm.eg')
level10_genes[,"Gene_Symbols"] <- as.character(symbols)
level10_genes <- as.data.frame(level10_genes)
index <- which(GOid2Terms$GO_ID%in%level10_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level10_genes$GO_ID%in%selectedGO$GO_ID)
level10_genes <- level10_genes[index,]
identical(nrow(selectedGO),length(unique(level10_genes$GO_ID)))
level10_genes <- left_join(level10_genes,GOid2Terms)
level10_genes <- level10_genes[,c(1,4,2,3,5)]
# 11
max.length <- max(sapply(level11_genes, length))
level11_genes <- lapply(level11_genes, function(v) { c(v, rep(NA, max.length-length(v)))})
level11_genes <- do.call(rbind, level11_genes)
level11_genes <- as.data.frame(level11_genes)
colnames(level11_genes) <- c(1:ncol(level11_genes))
level11_genes <- t(level11_genes)
level11_genes <- stack(level11_genes)
colnames(level11_genes) <- c("row","GO_ID","count","Entrez_ID")
level11_genes <- level11_genes[!is.na(level11_genes$Entrez_ID),]
level11_genes <- level11_genes[,-1]
level11_genes <- level11_genes[,-2]
level11_genes <- unique(level11_genes) # count: 11
symbols <- getSYMBOL(c(level11_genes$Entrez_ID),data='org.Mm.eg')
level11_genes[,"Gene_Symbols"] <- as.character(symbols)
level11_genes <- as.data.frame(level11_genes)
index <- which(GOid2Terms$GO_ID%in%level11_genes$GO_ID)
selectedGO <- as.data.frame(GOid2Terms[index,])
index <- which(level11_genes$GO_ID%in%selectedGO$GO_ID)
level11_genes <- level11_genes[index,]
identical(nrow(selectedGO),length(unique(level11_genes$GO_ID)))
level11_genes <- left_join(level11_genes,GOid2Terms)
level11_genes <- level11_genes[,c(1,4,2,3,5)]

rm(index,selectedGO,GOid2Terms,symbols,max.length)

options(warn = 0) # Turn warnings back on

level01_genes$GO_Level <- c(rep(c(1),times = NROW(level01_genes)))
level02_genes$GO_Level <- c(rep(c(2),times = NROW(level02_genes)))
level03_genes$GO_Level <- c(rep(c(3),times = NROW(level03_genes)))
level04_genes$GO_Level <- c(rep(c(4),times = NROW(level04_genes)))
level05_genes$GO_Level <- c(rep(c(5),times = NROW(level05_genes)))
level06_genes$GO_Level <- c(rep(c(6),times = NROW(level06_genes)))
level07_genes$GO_Level <- c(rep(c(7),times = NROW(level07_genes)))
level08_genes$GO_Level <- c(rep(c(8),times = NROW(level08_genes)))
level09_genes$GO_Level <- c(rep(c(9),times = NROW(level09_genes)))
level10_genes$GO_Level <- c(rep(c(10),times = NROW(level10_genes)))
level11_genes$GO_Level <- c(rep(c(11),times = NROW(level11_genes)))

GO.level.genes.mouse <- rbind(level01_genes,level02_genes,level03_genes,level04_genes,
                              level05_genes,level06_genes,level07_genes,level08_genes,
                              level09_genes,level10_genes,level11_genes)

save(GO.level.genes.mouse, file="Unique_Mouse_GO-Entrez_Pairs_NSG_28Sept2017.RData")

rm(list = ls(pattern = glob2rx("*_genes")))

