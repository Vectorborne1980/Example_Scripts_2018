library(dplyr)
library(GSEABase)
library(GSVA)
library(d3heatmap)
library(htmlwidgets)
library(limma)
setwd("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/GSVA/")

exprs.data <- read.delim("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/Exprs_Data_AdjP-0-05_rLog_NSG_12Oct2017.txt",sep = "\t")
exprs.data <- as.matrix(exprs.data)

### Run "Selection_of_GO_Term_Levels.R" for mouse or human as necessary
load("C:/Users/Nick Geraci/Nick/UVA/Sandbox/Gene_Ontology/Unique_Mouse_GO-Entrez_Pairs_NSG_28Sept2017.RData")
GO.level.genes.mouse[1] <- as.character(GO.level.genes.mouse[,1])

### Convert into usable format for GSVA, GAGE, DAVID, etc.
uniq<-unique(GO.level.genes.mouse[,1]) # get unique GO terms
GO.lists<-vector("list",length=length(uniq)) # set up data structure
names(GO.lists)<-uniq   # name the structure with the unique names
for(i in 1:length(GO.lists)) # loop through getting genes associated with each entry
  GO.lists[[i]]<-GO.level.genes.mouse[which(GO.level.genes.mouse[,1]==uniq[i]),3]
GO.lists[GO.lists==""] <- NA
GO.lists <- lapply(GO.lists,function(x) x[!is.na(x)])
uniq<-unique(uniq<-paste(GO.level.genes.mouse$GO_ID,GO.level.genes.mouse$GO_Term,GO.level.genes.mouse$GO_Definition,GO.level.genes.mouse$GO_Level,sep = "@"))
names(GO.lists)<-uniq
rm(i,uniq)

GO.names <- names(GO.lists)
uniqueList <- lapply(GO.lists, unique)

makeSet <- function(geneIds, GO.names) {
  GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=GO.names)
}

gsList <- mapply(makeSet, uniqueList[], GO.names)
newlist <- gsList %>% lapply(function(x) trimws(x@geneIds))
gsc <- GeneSetCollection(gsList)

Neutrophil.BM <- exprs.data[,c(1:2,5:6)]
Monocyte.BM <- exprs.data[,c(3:4,5:6)]
Mac.PC <- exprs.data[,c(7:8,5:6)]
Mac.Spleen <- exprs.data[,c(9:10,5:6)]
Kup.Liver <- exprs.data[,c(11:12,5:6)]
Mac.Colon <- exprs.data[,c(13:14,5:6)]
Mac.Ileum <- exprs.data[,c(15:16,5:6)]
Mac.Lung <- exprs.data[,c(17:18,5:6)]

### Run GSVA
results.gsva <- gsva(exprs.data, newlist, method = c("gsva"), rnaseq=TRUE,
                     min.sz=1, max.sz=Inf, verbose=TRUE)$es.obs

rm(nGeneSets,progressBar,iGeneSet,gsList,gsc,GO.names,uniqueList,GO.lists,newlist)

### Filter for significant results
design <- model.matrix(~ 0+factor(c(1,1,2,2)))
colnames(design) <- c("Immune_Cell","Microglia")
#Neutrophil.BM
rownames(design) <- c(colnames(results.gsva[,c(1:2,5:6)]))
fit <- lmFit(results.gsva[,c(1:2,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Neutrophil.BM <- allGeneSets[index,]
index <- which(top.sig.Neutrophil.BM$logFC > 0)
top.sig.Neutrophil.BM.UP <- top.sig.Neutrophil.BM[index,]
#Monocyte.BM
rownames(design) <- c(colnames(results.gsva[,c(3:4,5:6)]))
fit <- lmFit(results.gsva[,c(3:4,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Monocyte.BM <- allGeneSets[index,]
index <- which(top.sig.Monocyte.BM$logFC > 0)
top.sig.Monocyte.BM.UP <- top.sig.Monocyte.BM[index,]
#Mac.PC
rownames(design) <- c(colnames(results.gsva[,c(7:8,5:6)]))
fit <- lmFit(results.gsva[,c(7:8,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.PC <- allGeneSets[index,]
index <- which(top.sig.Mac.PC$logFC > 0)
top.sig.Mac.PC.UP <- top.sig.Mac.PC[index,]
#Mac.Spleen
rownames(design) <- c(colnames(results.gsva[,c(9:10,5:6)]))
fit <- lmFit(results.gsva[,c(9:10,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Spleen <- allGeneSets[index,]
index <- which(top.sig.Mac.Spleen$logFC > 0)
top.sig.Mac.Spleen.UP <- top.sig.Mac.Spleen[index,]
#Kup.Liver
rownames(design) <- c(colnames(results.gsva[,c(11:12,5:6)]))
fit <- lmFit(results.gsva[,c(11:12,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Kup.Liver <- allGeneSets[index,]
index <- which(top.sig.Kup.Liver$logFC > 0)
top.sig.Kup.Liver.UP <- top.sig.Kup.Liver[index,]
#Mac.Colon
rownames(design) <- c(colnames(results.gsva[,c(13:14,5:6)]))
fit <- lmFit(results.gsva[,c(13:14,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Colon <- allGeneSets[index,]
index <- which(top.sig.Mac.Colon$logFC > 0)
top.sig.Mac.Colon.UP <- top.sig.Mac.Colon[index,]
#Mac.Ileum
rownames(design) <- c(colnames(results.gsva[,c(15:16,5:6)]))
fit <- lmFit(results.gsva[,c(15:16,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Ileum <- allGeneSets[index,]
index <- which(top.sig.Mac.Ileum$logFC > 0)
top.sig.Mac.Ileum.UP <- top.sig.Mac.Ileum[index,]
#Mac.Lung
rownames(design) <- c(colnames(results.gsva[,c(17:18,5:6)]))
fit <- lmFit(results.gsva[,c(17:18,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Lung <- allGeneSets[index,]
index <- which(top.sig.Mac.Lung$logFC > 0)
top.sig.Mac.Lung.UP <- top.sig.Mac.Lung[index,]

save(results.gsva,
     top.sig.Neutrophil.BM,
     top.sig.Neutrophil.BM.UP,
     top.sig.Monocyte.BM,
     top.sig.Monocyte.BM.UP,
     top.sig.Mac.PC,
     top.sig.Mac.PC.UP,
     top.sig.Kup.Liver,
     top.sig.Kup.Liver.UP,
     top.sig.Mac.Spleen,
     top.sig.Mac.Spleen.UP,
     top.sig.Mac.Colon,
     top.sig.Mac.Colon.UP,
     top.sig.Mac.Ileum,
     top.sig.Mac.Ileum.UP,
     top.sig.Mac.Lung,
     top.sig.Mac.Lung.UP,
     file="GSVA_All-GO_Sig-by-LIMMA_NSG_12Oct2017.RData")
rm(fit,design,index,contrast.matrix,allGeneSets)

######## Make an UpSet plot of the significant GSVA enrichments
list.GO.sig <- list(Neutrophil.BM.UP = c(rownames(top.sig.Neutrophil.BM.UP)),
                    Monocyte.BM.UP = c(rownames(top.sig.Monocyte.BM.UP)),
                    Mac.PC.UP = c(rownames(top.sig.Mac.PC.UP)),
                    Kup.Liver.UP = c(rownames(top.sig.Kup.Liver.UP)),
                    Mac.Spleen.UP = c(rownames(top.sig.Mac.Spleen.UP)),
                    Mac.Colon.UP = c(rownames(top.sig.Mac.Colon.UP)),
                    Mac.Ileum.UP = c(rownames(top.sig.Mac.Ileum.UP)),
                    Mac.Lung.UP = c(rownames(top.sig.Mac.Lung.UP)))
max.length <- max(sapply(list.GO.sig, length))
list.GO.sig <- lapply(list.GO.sig, function(v) { c(v, rep(NA, max.length-length(v)))})
list.GO.sig <- do.call(rbind, list.GO.sig)
list.GO.sig <- t(list.GO.sig)
list.GO.sig <- as.data.frame(list.GO.sig)
list.GO.sig[,1:length(colnames(list.GO.sig))] <- sapply(list.GO.sig[,1:length(colnames(list.GO.sig))],as.character)
list.GO.sig2 <- stack(list.GO.sig)
list.GO.sig2 <- list.GO.sig2[!is.na(list.GO.sig2$values),]
list.GO.sig2 <- list.GO.sig2[order(list.GO.sig2$values),]
list.GO.sig2 <- as.data.frame(unique(list.GO.sig2$values),stringsAsFactors=FALSE)
rownames(list.GO.sig2) <- list.GO.sig2[,1]
list.GO.sig2$Neutrophil.BM.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Neutrophil.BM.UP),nomatch = 0)
list.GO.sig2$Monocyte.BM.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Monocyte.BM.UP),nomatch = 0)
list.GO.sig2$Mac.PC.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.PC.UP),nomatch = 0)
list.GO.sig2$Kup.Liver.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Kup.Liver.UP),nomatch = 0)
list.GO.sig2$Mac.Spleen.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Spleen.UP),nomatch = 0)
list.GO.sig2$Mac.Colon.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Colon.UP),nomatch = 0)
list.GO.sig2$Mac.Ileum.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Ileum.UP),nomatch = 0)
list.GO.sig2$Mac.Lung.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Lung.UP),nomatch = 0)
list.GO.sig2 <- list.GO.sig2[,2:length(colnames(list.GO.sig2))]
list.GO.sig2 <- apply(list.GO.sig2,1,function(x) {ifelse(x>0,1,0)})
list.GO.sig2 <- as.data.frame(t(list.GO.sig2))

library(UpSetR)
pdf("GSE63340_GSVA_GO_sig-UP_UpSet.pdf",height = 20,width = 100,onefile = FALSE)
upset(list.GO.sig2,nsets = length(colnames(list.GO.sig2)),nintersects = NA,sets = c(colnames(list.GO.sig2)),
      empty.intersections = "yes",order.by = "degree", show.numbers = "yes", group.by = "degree",
      mb.ratio = c(0.6, 0.4),number.angles = 0, point.size = 1, line.size = 0.5,
      mainbar.y.label = "Common Enriched GO Terms", sets.x.label = "Enriched GO Terms by Comparison", 
      text.scale = c(2, 1, 2, 1, 1.5, 0.5))
dev.off()

######## Make a heatmap for various significant comparisons of GSVA results
library(reshape)
library(xlsx)
wb <- createWorkbook()

sig_GS <- top.sig.Neutrophil.BM.UP
exprs_hm <- results.gsva[,c(1:2,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Neutrophil.BM-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Monocyte.BM.UP
exprs_hm <- results.gsva[,c(3:4,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Monocyte.BM-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.PC.UP
exprs_hm <- results.gsva[,c(7:8,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.PC-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Kup.Liver.UP
exprs_hm <- results.gsva[,c(9:10,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Kup.Liver-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Spleen.UP
exprs_hm <- results.gsva[,c(11:12,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Spleen-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Colon.UP
exprs_hm <- results.gsva[,c(13:14,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Colon-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Ileum.UP
exprs_hm <- results.gsva[,c(15:16,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Ileum-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Lung.UP
exprs_hm <- results.gsva[,c(17:18,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Lung-over-Microglia_UP")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)


saveWorkbook(wb,file="GSE63340_GSVA_GO_sig-UP_NSG_12Oct2017.xlsx")
rm(wb,ws)

## Take just one GO level if desired
# exprs_hm$GO.Terms <- rownames(exprs_hm) #Turn the row names into a column
# exprs_hm <- cbind(exprs_hm[,1:6],GO.Terms <- 
#       data.frame(do.call('rbind', 
#       strsplit(as.character(exprs_hm$GO.Terms),';',fixed=TRUE)))) #Split the GO terms by the semicolon into separate columns
# index <- which(exprs_hm$X3==2) #Take whichever number level you desire
# exprs_hm <- exprs_hm[index,1:6]

## Get average enrichment value for the two groups
# exprs_hm <- cbind(exprs_hm, avg <- rowMeans(exprs_hm[,1:3]))
# exprs_hm <- cbind(exprs_hm, avg <- rowMeans(exprs_hm[,4:6]))
# colnames(exprs_hm)[7] <- "Plus.Inact.Avg"
# colnames(exprs_hm)[8] <- "Plus.Act.Avg"
# 
# exprs_hm$Plus.Inact <- apply(exprs_hm[,7,drop=F],1,function(x) {ifelse(x>0,1,0)})
# exprs_hm$Plus.Act <- apply(exprs_hm[,8,drop=F],1,function(x) {ifelse(x>0,1,0)})

# index <- which(is.na(rowSums(exprs_hm)))
# if (length(index) > 0) { exprs_hm <- exprs_hm[-index,] }
# 
# distance <- function(x, ...) dist(x,method="euclidean")
# samples = colnames(exprs_hm)
# 
# hmcol <- colorRampPalette(c("blue","white","red"))(256)
# 
# hm <- d3heatmap(data.matrix(exprs_hm),
#                 colors = hmcol,
#                 scale = "row",
#                 dendrogram = "both",
#                 hclustfun = hclust,
#                 distfun = distance,
#                 k_row = 10,
#                 k_col = 2,
#                 show_grid = 1,
#                 height = 10000,
#                 width = 2000,
#                 xaxis_font_size = "15pt",
#                 yaxis_font_size = "2pt",
#                 xaxis_height = 150,
#                 yaxis_width = 1000)
# 
# saveWidget(hm,"Plus-vs-Minus_GSVA_GO_sig-UP.html")



########### Repeat the analyses for KEGG Pathways
setwd("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/GSVA/")

exprs.data <- read.delim("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/Exprs_Data_AdjP-0-05_rLog_NSG_12Oct2017.txt",sep = "\t")
exprs.data <- as.matrix(exprs.data)

mouse.kegg <- read.delim("MouseKEGG.txt",row.names = 1,header = TRUE,sep = "\t")
mouse.kegg$combo.ID <- paste(mouse.kegg$pathway,mouse.kegg$pathway.desc,sep = ";")
mouse.kegg$entrez <- as.character(mouse.kegg$entrez)
k<-mouse.kegg[,c(1,7)] # select the kegg columns we need
uniq<-unique(k[,2]) # get unique kegg descriptions
kegglists<-vector("list",length=length(uniq)) # set up data structure
names(kegglists)<-uniq   # name the structure with the unique names
for(i in 1:length(kegglists)) # loop through getting genes associated with each kegg entry
  kegglists[[i]]<-k[which(k[,2]==uniq[i]),1]
kegglists[kegglists==""] <- NA
kegglists <- lapply(kegglists,function(x) x[!is.na(x)])
KEGG.names <- names(kegglists)
uniqueList <- lapply(kegglists, unique)
makeSet <- function(geneIds, KEGG.names) {
  GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=KEGG.names)
}

gsList <- mapply(makeSet, uniqueList[], KEGG.names)
newlist <- gsList %>% lapply(function(x) trimws(x@geneIds))
gsc <- GeneSetCollection(gsList)

Neutrophil.BM <- exprs.data[,c(1:2,5:6)]
Monocyte.BM <- exprs.data[,c(3:4,5:6)]
Mac.PC <- exprs.data[,c(7:8,5:6)]
Mac.Spleen <- exprs.data[,c(9:10,5:6)]
Kup.Liver <- exprs.data[,c(11:12,5:6)]
Mac.Colon <- exprs.data[,c(13:14,5:6)]
Mac.Ileum <- exprs.data[,c(15:16,5:6)]
Mac.Lung <- exprs.data[,c(17:18,5:6)]

### Run GSVA
results.gsva <- gsva(exprs.data, newlist, method = c("gsva"), rnaseq=TRUE,
                     min.sz=1, max.sz=Inf, verbose=TRUE)$es.obs

rm(nGeneSets,progressBar,iGeneSet,gsList,gsc,KEGG.names,uniqueList,KEGG.lists,newlist)

### Filter for significant results
design <- model.matrix(~ 0+factor(c(1,1,2,2)))
colnames(design) <- c("Immune_Cell","Microglia")
#Neutrophil.BM
rownames(design) <- c(colnames(results.gsva[,c(1:2,5:6)]))
fit <- lmFit(results.gsva[,c(1:2,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Neutrophil.BM <- allGeneSets[index,]
index <- which(top.sig.Neutrophil.BM$logFC > 0)
top.sig.Neutrophil.BM.UP <- top.sig.Neutrophil.BM[index,]
#Monocyte.BM
rownames(design) <- c(colnames(results.gsva[,c(3:4,5:6)]))
fit <- lmFit(results.gsva[,c(3:4,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Monocyte.BM <- allGeneSets[index,]
index <- which(top.sig.Monocyte.BM$logFC > 0)
top.sig.Monocyte.BM.UP <- top.sig.Monocyte.BM[index,]
#Mac.PC
rownames(design) <- c(colnames(results.gsva[,c(7:8,5:6)]))
fit <- lmFit(results.gsva[,c(7:8,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.PC <- allGeneSets[index,]
index <- which(top.sig.Mac.PC$logFC > 0)
top.sig.Mac.PC.UP <- top.sig.Mac.PC[index,]
#Mac.Spleen
rownames(design) <- c(colnames(results.gsva[,c(9:10,5:6)]))
fit <- lmFit(results.gsva[,c(9:10,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Spleen <- allGeneSets[index,]
index <- which(top.sig.Mac.Spleen$logFC > 0)
top.sig.Mac.Spleen.UP <- top.sig.Mac.Spleen[index,]
#Kup.Liver
rownames(design) <- c(colnames(results.gsva[,c(11:12,5:6)]))
fit <- lmFit(results.gsva[,c(11:12,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Kup.Liver <- allGeneSets[index,]
index <- which(top.sig.Kup.Liver$logFC > 0)
top.sig.Kup.Liver.UP <- top.sig.Kup.Liver[index,]
#Mac.Colon
rownames(design) <- c(colnames(results.gsva[,c(13:14,5:6)]))
fit <- lmFit(results.gsva[,c(13:14,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Colon <- allGeneSets[index,]
index <- which(top.sig.Mac.Colon$logFC > 0)
top.sig.Mac.Colon.UP <- top.sig.Mac.Colon[index,]
#Mac.Ileum
rownames(design) <- c(colnames(results.gsva[,c(15:16,5:6)]))
fit <- lmFit(results.gsva[,c(15:16,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Ileum <- allGeneSets[index,]
index <- which(top.sig.Mac.Ileum$logFC > 0)
top.sig.Mac.Ileum.UP <- top.sig.Mac.Ileum[index,]
#Mac.Lung
rownames(design) <- c(colnames(results.gsva[,c(17:18,5:6)]))
fit <- lmFit(results.gsva[,c(17:18,5:6)], design)
contrast.matrix <- makeContrasts(Immune_Cell-Microglia,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Lung <- allGeneSets[index,]
index <- which(top.sig.Mac.Lung$logFC > 0)
top.sig.Mac.Lung.UP <- top.sig.Mac.Lung[index,]

save(results.gsva,
     top.sig.Neutrophil.BM,
     top.sig.Neutrophil.BM.UP,
     top.sig.Monocyte.BM,
     top.sig.Monocyte.BM.UP,
     top.sig.Mac.PC,
     top.sig.Mac.PC.UP,
     top.sig.Kup.Liver,
     top.sig.Kup.Liver.UP,
     top.sig.Mac.Spleen,
     top.sig.Mac.Spleen.UP,
     top.sig.Mac.Colon,
     top.sig.Mac.Colon.UP,
     top.sig.Mac.Ileum,
     top.sig.Mac.Ileum.UP,
     top.sig.Mac.Lung,
     top.sig.Mac.Lung.UP,
     file="GSVA_All-KEGG_Sig-by-LIMMA_NSG_12Oct2017.RData")
rm(fit,design,index,contrast.matrix,allGeneSets)

######## Make an UpSet plot of the significant GSVA enrichments
list.KEGG.sig <- list(Neutrophil.BM.UP = c(rownames(top.sig.Neutrophil.BM.UP)),
                    Monocyte.BM.UP = c(rownames(top.sig.Monocyte.BM.UP)),
                    Mac.PC.UP = c(rownames(top.sig.Mac.PC.UP)),
                    Kup.Liver.UP = c(rownames(top.sig.Kup.Liver.UP)),
                    Mac.Spleen.UP = c(rownames(top.sig.Mac.Spleen.UP)),
                    Mac.Colon.UP = c(rownames(top.sig.Mac.Colon.UP)),
                    Mac.Ileum.UP = c(rownames(top.sig.Mac.Ileum.UP)),
                    Mac.Lung.UP = c(rownames(top.sig.Mac.Lung.UP)))
max.length <- max(sapply(list.KEGG.sig, length))
list.KEGG.sig <- lapply(list.KEGG.sig, function(v) { c(v, rep(NA, max.length-length(v)))})
list.KEGG.sig <- do.call(rbind, list.KEGG.sig)
list.KEGG.sig <- t(list.KEGG.sig)
list.KEGG.sig <- as.data.frame(list.KEGG.sig)
list.KEGG.sig[,1:length(colnames(list.KEGG.sig))] <- sapply(list.KEGG.sig[,1:length(colnames(list.KEGG.sig))],as.character)
list.KEGG.sig2 <- stack(list.KEGG.sig)
list.KEGG.sig2 <- list.KEGG.sig2[!is.na(list.KEGG.sig2$values),]
list.KEGG.sig2 <- list.KEGG.sig2[order(list.KEGG.sig2$values),]
list.KEGG.sig2 <- as.data.frame(unique(list.KEGG.sig2$values),stringsAsFactors=FALSE)
rownames(list.KEGG.sig2) <- list.KEGG.sig2[,1]
list.KEGG.sig2$Neutrophil.BM.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Neutrophil.BM.UP),nomatch = 0)
list.KEGG.sig2$Monocyte.BM.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Monocyte.BM.UP),nomatch = 0)
list.KEGG.sig2$Mac.PC.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.PC.UP),nomatch = 0)
list.KEGG.sig2$Kup.Liver.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Kup.Liver.UP),nomatch = 0)
list.KEGG.sig2$Mac.Spleen.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Spleen.UP),nomatch = 0)
list.KEGG.sig2$Mac.Colon.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Colon.UP),nomatch = 0)
list.KEGG.sig2$Mac.Ileum.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Ileum.UP),nomatch = 0)
list.KEGG.sig2$Mac.Lung.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Lung.UP),nomatch = 0)
list.KEGG.sig2 <- list.KEGG.sig2[,2:length(colnames(list.KEGG.sig2))]
list.KEGG.sig2 <- apply(list.KEGG.sig2,1,function(x) {ifelse(x>0,1,0)})
list.KEGG.sig2 <- as.data.frame(t(list.KEGG.sig2))

library(UpSetR)
pdf("GSE63340_GSVA_KEGG_sig-UP_UpSet.pdf",height = 20,width = 100,onefile = FALSE)
upset(list.KEGG.sig2,nsets = length(colnames(list.KEGG.sig2)),nintersects = NA,sets = c(colnames(list.KEGG.sig2)),
      empty.intersections = "yes",order.by = "degree", show.numbers = "yes", group.by = "degree",
      mb.ratio = c(0.6, 0.4),number.angles = 0, point.size = 1, line.size = 0.5,
      mainbar.y.label = "Common Enriched KEGG Terms", sets.x.label = "Enriched KEGG Terms by Comparison", 
      text.scale = c(2, 1, 2, 1, 1.5, 0.5))
dev.off()

######## Make a heatmap for various significant comparisons of GSVA results
library(reshape)
library(xlsx)
wb <- createWorkbook()

sig_GS <- top.sig.Neutrophil.BM.UP
exprs_hm <- results.gsva[,c(1:2,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Neutrophil.BM-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Monocyte.BM.UP
exprs_hm <- results.gsva[,c(3:4,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Monocyte.BM-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.PC.UP
exprs_hm <- results.gsva[,c(7:8,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.PC-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Kup.Liver.UP
exprs_hm <- results.gsva[,c(9:10,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Kup.Liver-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Spleen.UP
exprs_hm <- results.gsva[,c(11:12,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Spleen-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Colon.UP
exprs_hm <- results.gsva[,c(13:14,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Colon-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Ileum.UP
exprs_hm <- results.gsva[,c(15:16,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Ileum-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Lung.UP
exprs_hm <- results.gsva[,c(17:18,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Mac.Lung-over-Microglia_UP")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)


saveWorkbook(wb,file="GSE63340_GSVA_KEGG_sig-UP_NSG_12Oct2017.xlsx")
rm(wb,ws)

# index <- which(is.na(rowSums(exprs_hm)))
# if (length(index) > 0) { exprs_hm <- exprs_hm[-index,] }
# 
# distance <- dist(t(exprs_hm),method="euclidean")
# hc <- hclust(distance)
# 
# dist2 <- function(x, ...)
#   dist(x,method="euclidean")
# label.size = 1.2
# samples = colnames(exprs_hm)
# 
# hmcol <- colorRampPalette(c("blue","white","red"))(256)
# 
# hm <- d3heatmap(data.matrix(exprs_hm),
#                 colors = hmcol, scale = "row",
#                 dendrogram = "both",
#                 hclustfun = hclust,
#                 distfun = dist2,
#                 k_row = 10,
#                 k_col = 2,
#                 show_grid = 1,
#                 height = 4000,
#                 width = 1000,
#                 xaxis_font_size = "15pt",
#                 yaxis_font_size = "10pt",
#                 xaxis_height = 150,
#                 yaxis_width = 500)
# 
# saveWidget(hm,"Plus-vs-Minus_GSVA_KEGG_sig-UP.html")



######## GAGE
# library(AnnotationDbi)
# library(org.Mm.eg.db)
# library(gage)
# library(gageData)
# library(pathview)
# 
# setwd("C:/Users/Nick Geraci/Desktop/UVA/Bulk_RNA-seq/")
# 
# load("C:/Users/Nick Geraci/Desktop/UVA/Bulk_RNA-seq/gpm_RNAseq_DESeq2_Processed_Data_NSG_26Sept2017.RData")
# 
# columns(org.Mm.eg.db)
# res.Plus.sig$ENTREZID <- mapIds(org.Mm.eg.db,
#                               keys=row.names(res.Plus.sig),
#                               column="ENTREZID",
#                               keytype="ENSEMBL",
#                               multiVals="first")
# 
# res.Plus.sig$SYMBOL <- mapIds(org.Mm.eg.db,
#                             keys=row.names(res.Plus.sig),
#                             column="SYMBOL",
#                             keytype="ENSEMBL",
#                             multiVals="first")
# 
# res.Plus.sig$ENSEMBL <- mapIds(org.Mm.eg.db,
#                              keys=row.names(res.Plus.sig),
#                              column="ENSEMBL",
#                              keytype="ENSEMBL",
#                              multiVals="first")
# 
# data(kegg.sets.mm)
# data(sigmet.idx.mm)
# kegg.sets.mm <- kegg.sets.mm[sigmet.idx.mm]
# head(kegg.sets.mm, 3)
# 
# foldchanges <- res.Plus.sig$log2FoldChange
# names(foldchanges) <- res.Plus.sig$ENTREZID
# head(foldchanges)
# 
# # Get the results
# keggres <- gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
# 
# # Look at both up (greater), down (less), and statatistics.
# lapply(keggres, head)
# 
# keggres.sig <- as.data.frame(keggres$greater)
# keggres.sig <- keggres.sig[keggres.sig$q.val<0.1,]
# keggres.sig <- keggres.sig[1:37,]
# keggrespathways <- as.character(rownames(keggres.sig))
# keggresids <- substr(keggrespathways, start=1, stop=8)
# keggresids
# # Define plotting function for applying later
# plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
# 
# # plot multiple pathways (plots saved to disk and returns a throwaway list object)
# tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
# tmp <- pathview(gene.data=foldchanges, pathway.id="03013", species="mmu")

########################################################################################################################################################################
####Repeat the analyses using Microglia as the numerator################################################################################################################
########################################################################################################################################################################

library(dplyr)
library(GSEABase)
library(GSVA)
library(d3heatmap)
library(htmlwidgets)
library(limma)
setwd("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/GSVA/")

exprs.data <- read.delim("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/Exprs_Data_AdjP-0-05_rLog_NSG_12Oct2017.txt",sep = "\t")
exprs.data <- as.matrix(exprs.data)

### Run "Selection_of_GO_Term_Levels.R" for mouse or human as necessary
load("C:/Users/Nick Geraci/Nick/UVA/Sandbox/Gene_Ontology/Unique_Mouse_GO-Entrez_Pairs_NSG_28Sept2017.RData")
GO.level.genes.mouse[1] <- as.character(GO.level.genes.mouse[,1])

### Convert into usable format for GSVA, GAGE, DAVID, etc.
uniq<-unique(GO.level.genes.mouse[,1]) # get unique GO terms
GO.lists<-vector("list",length=length(uniq)) # set up data structure
names(GO.lists)<-uniq   # name the structure with the unique names
for(i in 1:length(GO.lists)) # loop through getting genes associated with each entry
  GO.lists[[i]]<-GO.level.genes.mouse[which(GO.level.genes.mouse[,1]==uniq[i]),3]
GO.lists[GO.lists==""] <- NA
GO.lists <- lapply(GO.lists,function(x) x[!is.na(x)])
uniq<-unique(uniq<-paste(GO.level.genes.mouse$GO_ID,GO.level.genes.mouse$GO_Term,GO.level.genes.mouse$GO_Definition,GO.level.genes.mouse$GO_Level,sep = "@"))
names(GO.lists)<-uniq
rm(i,uniq)

GO.names <- names(GO.lists)
uniqueList <- lapply(GO.lists, unique)

makeSet <- function(geneIds, GO.names) {
  GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=GO.names)
}

gsList <- mapply(makeSet, uniqueList[], GO.names)
newlist <- gsList %>% lapply(function(x) trimws(x@geneIds))
gsc <- GeneSetCollection(gsList)

Neutrophil.BM <- exprs.data[,c(1:2,5:6)]
Monocyte.BM <- exprs.data[,c(3:4,5:6)]
Mac.PC <- exprs.data[,c(7:8,5:6)]
Mac.Spleen <- exprs.data[,c(9:10,5:6)]
Kup.Liver <- exprs.data[,c(11:12,5:6)]
Mac.Colon <- exprs.data[,c(13:14,5:6)]
Mac.Ileum <- exprs.data[,c(15:16,5:6)]
Mac.Lung <- exprs.data[,c(17:18,5:6)]

### Run GSVA
results.gsva <- gsva(exprs.data, newlist, method = c("gsva"), rnaseq=TRUE,
                     min.sz=1, max.sz=Inf, verbose=TRUE)$es.obs

rm(nGeneSets,progressBar,iGeneSet,gsList,gsc,GO.names,uniqueList,GO.lists,newlist)

### Filter for significant results
design <- model.matrix(~ 0+factor(c(1,1,2,2)))
colnames(design) <- c("Immune_Cell","Microglia")
#Neutrophil.BM
rownames(design) <- c(colnames(results.gsva[,c(1:2,5:6)]))
fit <- lmFit(results.gsva[,c(1:2,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Neutrophil.BM <- allGeneSets[index,]
index <- which(top.sig.Neutrophil.BM$logFC > 0)
top.sig.Neutrophil.BM.UP <- top.sig.Neutrophil.BM[index,]
#Monocyte.BM
rownames(design) <- c(colnames(results.gsva[,c(3:4,5:6)]))
fit <- lmFit(results.gsva[,c(3:4,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Monocyte.BM <- allGeneSets[index,]
index <- which(top.sig.Monocyte.BM$logFC > 0)
top.sig.Monocyte.BM.UP <- top.sig.Monocyte.BM[index,]
#Mac.PC
rownames(design) <- c(colnames(results.gsva[,c(7:8,5:6)]))
fit <- lmFit(results.gsva[,c(7:8,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.PC <- allGeneSets[index,]
index <- which(top.sig.Mac.PC$logFC > 0)
top.sig.Mac.PC.UP <- top.sig.Mac.PC[index,]
#Mac.Spleen
rownames(design) <- c(colnames(results.gsva[,c(9:10,5:6)]))
fit <- lmFit(results.gsva[,c(9:10,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Spleen <- allGeneSets[index,]
index <- which(top.sig.Mac.Spleen$logFC > 0)
top.sig.Mac.Spleen.UP <- top.sig.Mac.Spleen[index,]
#Kup.Liver
rownames(design) <- c(colnames(results.gsva[,c(11:12,5:6)]))
fit <- lmFit(results.gsva[,c(11:12,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Kup.Liver <- allGeneSets[index,]
index <- which(top.sig.Kup.Liver$logFC > 0)
top.sig.Kup.Liver.UP <- top.sig.Kup.Liver[index,]
#Mac.Colon
rownames(design) <- c(colnames(results.gsva[,c(13:14,5:6)]))
fit <- lmFit(results.gsva[,c(13:14,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Colon <- allGeneSets[index,]
index <- which(top.sig.Mac.Colon$logFC > 0)
top.sig.Mac.Colon.UP <- top.sig.Mac.Colon[index,]
#Mac.Ileum
rownames(design) <- c(colnames(results.gsva[,c(15:16,5:6)]))
fit <- lmFit(results.gsva[,c(15:16,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Ileum <- allGeneSets[index,]
index <- which(top.sig.Mac.Ileum$logFC > 0)
top.sig.Mac.Ileum.UP <- top.sig.Mac.Ileum[index,]
#Mac.Lung
rownames(design) <- c(colnames(results.gsva[,c(17:18,5:6)]))
fit <- lmFit(results.gsva[,c(17:18,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Lung <- allGeneSets[index,]
index <- which(top.sig.Mac.Lung$logFC > 0)
top.sig.Mac.Lung.UP <- top.sig.Mac.Lung[index,]

save(results.gsva,
     top.sig.Neutrophil.BM,
     top.sig.Neutrophil.BM.UP,
     top.sig.Monocyte.BM,
     top.sig.Monocyte.BM.UP,
     top.sig.Mac.PC,
     top.sig.Mac.PC.UP,
     top.sig.Kup.Liver,
     top.sig.Kup.Liver.UP,
     top.sig.Mac.Spleen,
     top.sig.Mac.Spleen.UP,
     top.sig.Mac.Colon,
     top.sig.Mac.Colon.UP,
     top.sig.Mac.Ileum,
     top.sig.Mac.Ileum.UP,
     top.sig.Mac.Lung,
     top.sig.Mac.Lung.UP,
     file="GSE63340_GSVA-Rev_All-GO_Sig-by-LIMMA_NSG_12Oct2017.RData")
rm(fit,design,index,contrast.matrix,allGeneSets)

######## Make an UpSet plot of the significant GSVA enrichments
list.GO.sig <- list(Neutrophil.BM.UP = c(rownames(top.sig.Neutrophil.BM.UP)),
                    Monocyte.BM.UP = c(rownames(top.sig.Monocyte.BM.UP)),
                    Mac.PC.UP = c(rownames(top.sig.Mac.PC.UP)),
                    Kup.Liver.UP = c(rownames(top.sig.Kup.Liver.UP)),
                    Mac.Spleen.UP = c(rownames(top.sig.Mac.Spleen.UP)),
                    Mac.Colon.UP = c(rownames(top.sig.Mac.Colon.UP)),
                    Mac.Ileum.UP = c(rownames(top.sig.Mac.Ileum.UP)),
                    Mac.Lung.UP = c(rownames(top.sig.Mac.Lung.UP)))
max.length <- max(sapply(list.GO.sig, length))
list.GO.sig <- lapply(list.GO.sig, function(v) { c(v, rep(NA, max.length-length(v)))})
list.GO.sig <- do.call(rbind, list.GO.sig)
list.GO.sig <- t(list.GO.sig)
list.GO.sig <- as.data.frame(list.GO.sig)
list.GO.sig[,1:length(colnames(list.GO.sig))] <- sapply(list.GO.sig[,1:length(colnames(list.GO.sig))],as.character)
list.GO.sig2 <- stack(list.GO.sig)
list.GO.sig2 <- list.GO.sig2[!is.na(list.GO.sig2$values),]
list.GO.sig2 <- list.GO.sig2[order(list.GO.sig2$values),]
list.GO.sig2 <- as.data.frame(unique(list.GO.sig2$values),stringsAsFactors=FALSE)
rownames(list.GO.sig2) <- list.GO.sig2[,1]
list.GO.sig2$Neutrophil.BM.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Neutrophil.BM.UP),nomatch = 0)
list.GO.sig2$Monocyte.BM.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Monocyte.BM.UP),nomatch = 0)
list.GO.sig2$Mac.PC.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.PC.UP),nomatch = 0)
list.GO.sig2$Kup.Liver.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Kup.Liver.UP),nomatch = 0)
list.GO.sig2$Mac.Spleen.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Spleen.UP),nomatch = 0)
list.GO.sig2$Mac.Colon.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Colon.UP),nomatch = 0)
list.GO.sig2$Mac.Ileum.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Ileum.UP),nomatch = 0)
list.GO.sig2$Mac.Lung.UP <- match(c(rownames(list.GO.sig2)),c(list.GO.sig$Mac.Lung.UP),nomatch = 0)
list.GO.sig2 <- list.GO.sig2[,2:length(colnames(list.GO.sig2))]
list.GO.sig2 <- apply(list.GO.sig2,1,function(x) {ifelse(x>0,1,0)})
list.GO.sig2 <- as.data.frame(t(list.GO.sig2))

library(UpSetR)
pdf("GSE63340_GSVA-Rev_GO_sig-UP_UpSet.pdf",height = 20,width = 100,onefile = FALSE)
upset(list.GO.sig2,nsets = length(colnames(list.GO.sig2)),nintersects = NA,sets = c(colnames(list.GO.sig2)),
      empty.intersections = "yes",order.by = "degree", show.numbers = "yes", group.by = "degree",
      mb.ratio = c(0.6, 0.4),number.angles = 0, point.size = 1, line.size = 0.5,
      mainbar.y.label = "Common Enriched GO Terms", sets.x.label = "Enriched GO Terms by Comparison", 
      text.scale = c(2, 1, 2, 1, 1.5, 0.5))
dev.off()

######## Make a heatmap for various significant comparisons of GSVA results
library(reshape)
library(xlsx)
wb <- createWorkbook()

sig_GS <- top.sig.Neutrophil.BM.UP
exprs_hm <- results.gsva[,c(1:2,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Neutrophil.BM")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Monocyte.BM.UP
exprs_hm <- results.gsva[,c(3:4,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Monocyte.BM")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.PC.UP
exprs_hm <- results.gsva[,c(7:8,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.PC")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Kup.Liver.UP
exprs_hm <- results.gsva[,c(9:10,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Kup.Liver")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Spleen.UP
exprs_hm <- results.gsva[,c(11:12,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Spleen")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Colon.UP
exprs_hm <- results.gsva[,c(13:14,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Colon")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Ileum.UP
exprs_hm <- results.gsva[,c(15:16,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Ileum")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Lung.UP
exprs_hm <- results.gsva[,c(17:18,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Lung")
df <- transform(exprs_hm, GO = colsplit(rownames(exprs_hm), split = "\\@", names = c('ID', 'Term', 'Definition','Level')))
df <- df[,c(5:8,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)


saveWorkbook(wb,file="GSE63340_GSVA-Rev_GO_sig-UP_NSG_12Oct2017.xlsx")
rm(wb,ws)

## Take just one GO level if desired
# exprs_hm$GO.Terms <- rownames(exprs_hm) #Turn the row names into a column
# exprs_hm <- cbind(exprs_hm[,1:6],GO.Terms <- 
#       data.frame(do.call('rbind', 
#       strsplit(as.character(exprs_hm$GO.Terms),';',fixed=TRUE)))) #Split the GO terms by the semicolon into separate columns
# index <- which(exprs_hm$X3==2) #Take whichever number level you desire
# exprs_hm <- exprs_hm[index,1:6]

## Get average enrichment value for the two groups
# exprs_hm <- cbind(exprs_hm, avg <- rowMeans(exprs_hm[,1:3]))
# exprs_hm <- cbind(exprs_hm, avg <- rowMeans(exprs_hm[,4:6]))
# colnames(exprs_hm)[7] <- "Plus.Inact.Avg"
# colnames(exprs_hm)[8] <- "Plus.Act.Avg"
# 
# exprs_hm$Plus.Inact <- apply(exprs_hm[,7,drop=F],1,function(x) {ifelse(x>0,1,0)})
# exprs_hm$Plus.Act <- apply(exprs_hm[,8,drop=F],1,function(x) {ifelse(x>0,1,0)})

# index <- which(is.na(rowSums(exprs_hm)))
# if (length(index) > 0) { exprs_hm <- exprs_hm[-index,] }
# 
# distance <- function(x, ...) dist(x,method="euclidean")
# samples = colnames(exprs_hm)
# 
# hmcol <- colorRampPalette(c("blue","white","red"))(256)
# 
# hm <- d3heatmap(data.matrix(exprs_hm),
#                 colors = hmcol,
#                 scale = "row",
#                 dendrogram = "both",
#                 hclustfun = hclust,
#                 distfun = distance,
#                 k_row = 10,
#                 k_col = 2,
#                 show_grid = 1,
#                 height = 10000,
#                 width = 2000,
#                 xaxis_font_size = "15pt",
#                 yaxis_font_size = "2pt",
#                 xaxis_height = 150,
#                 yaxis_width = 1000)
# 
# saveWidget(hm,"Plus-vs-Minus_GSVA_GO_sig-UP.html")



########### Repeat the analyses for KEGG Pathways
setwd("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/GSVA/")

exprs.data <- read.delim("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/Exprs_Data_AdjP-0-05_rLog_NSG_12Oct2017.txt",sep = "\t")
exprs.data <- as.matrix(exprs.data)

mouse.kegg <- read.delim("MouseKEGG.txt",row.names = 1,header = TRUE,sep = "\t")
mouse.kegg$combo.ID <- paste(mouse.kegg$pathway,mouse.kegg$pathway.desc,sep = ";")
mouse.kegg$entrez <- as.character(mouse.kegg$entrez)
k<-mouse.kegg[,c(1,7)] # select the kegg columns we need
uniq<-unique(k[,2]) # get unique kegg descriptions
kegglists<-vector("list",length=length(uniq)) # set up data structure
names(kegglists)<-uniq   # name the structure with the unique names
for(i in 1:length(kegglists)) # loop through getting genes associated with each kegg entry
  kegglists[[i]]<-k[which(k[,2]==uniq[i]),1]
kegglists[kegglists==""] <- NA
kegglists <- lapply(kegglists,function(x) x[!is.na(x)])
KEGG.names <- names(kegglists)
uniqueList <- lapply(kegglists, unique)
makeSet <- function(geneIds, KEGG.names) {
  GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=KEGG.names)
}

gsList <- mapply(makeSet, uniqueList[], KEGG.names)
newlist <- gsList %>% lapply(function(x) trimws(x@geneIds))
gsc <- GeneSetCollection(gsList)

Neutrophil.BM <- exprs.data[,c(1:2,5:6)]
Monocyte.BM <- exprs.data[,c(3:4,5:6)]
Mac.PC <- exprs.data[,c(7:8,5:6)]
Mac.Spleen <- exprs.data[,c(9:10,5:6)]
Kup.Liver <- exprs.data[,c(11:12,5:6)]
Mac.Colon <- exprs.data[,c(13:14,5:6)]
Mac.Ileum <- exprs.data[,c(15:16,5:6)]
Mac.Lung <- exprs.data[,c(17:18,5:6)]

### Run GSVA
results.gsva <- gsva(exprs.data, newlist, method = c("gsva"), rnaseq=TRUE,
                     min.sz=1, max.sz=Inf, verbose=TRUE)$es.obs

rm(nGeneSets,progressBar,iGeneSet,gsList,gsc,KEGG.names,uniqueList,KEGG.lists,newlist)

### Filter for significant results
design <- model.matrix(~ 0+factor(c(1,1,2,2)))
colnames(design) <- c("Immune_Cell","Microglia")
#Neutrophil.BM
rownames(design) <- c(colnames(results.gsva[,c(1:2,5:6)]))
fit <- lmFit(results.gsva[,c(1:2,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Neutrophil.BM <- allGeneSets[index,]
index <- which(top.sig.Neutrophil.BM$logFC > 0)
top.sig.Neutrophil.BM.UP <- top.sig.Neutrophil.BM[index,]
#Monocyte.BM
rownames(design) <- c(colnames(results.gsva[,c(3:4,5:6)]))
fit <- lmFit(results.gsva[,c(3:4,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Monocyte.BM <- allGeneSets[index,]
index <- which(top.sig.Monocyte.BM$logFC > 0)
top.sig.Monocyte.BM.UP <- top.sig.Monocyte.BM[index,]
#Mac.PC
rownames(design) <- c(colnames(results.gsva[,c(7:8,5:6)]))
fit <- lmFit(results.gsva[,c(7:8,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.PC <- allGeneSets[index,]
index <- which(top.sig.Mac.PC$logFC > 0)
top.sig.Mac.PC.UP <- top.sig.Mac.PC[index,]
#Mac.Spleen
rownames(design) <- c(colnames(results.gsva[,c(9:10,5:6)]))
fit <- lmFit(results.gsva[,c(9:10,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Spleen <- allGeneSets[index,]
index <- which(top.sig.Mac.Spleen$logFC > 0)
top.sig.Mac.Spleen.UP <- top.sig.Mac.Spleen[index,]
#Kup.Liver
rownames(design) <- c(colnames(results.gsva[,c(11:12,5:6)]))
fit <- lmFit(results.gsva[,c(11:12,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Kup.Liver <- allGeneSets[index,]
index <- which(top.sig.Kup.Liver$logFC > 0)
top.sig.Kup.Liver.UP <- top.sig.Kup.Liver[index,]
#Mac.Colon
rownames(design) <- c(colnames(results.gsva[,c(13:14,5:6)]))
fit <- lmFit(results.gsva[,c(13:14,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Colon <- allGeneSets[index,]
index <- which(top.sig.Mac.Colon$logFC > 0)
top.sig.Mac.Colon.UP <- top.sig.Mac.Colon[index,]
#Mac.Ileum
rownames(design) <- c(colnames(results.gsva[,c(15:16,5:6)]))
fit <- lmFit(results.gsva[,c(15:16,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Ileum <- allGeneSets[index,]
index <- which(top.sig.Mac.Ileum$logFC > 0)
top.sig.Mac.Ileum.UP <- top.sig.Mac.Ileum[index,]
#Mac.Lung
rownames(design) <- c(colnames(results.gsva[,c(17:18,5:6)]))
fit <- lmFit(results.gsva[,c(17:18,5:6)], design)
contrast.matrix <- makeContrasts(Microglia-Immune_Cell,levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=1, adjust="BH", n=nrow(results.gsva))
index <- which(allGeneSets$adj.P.Val < 0.05)
top.sig.Mac.Lung <- allGeneSets[index,]
index <- which(top.sig.Mac.Lung$logFC > 0)
top.sig.Mac.Lung.UP <- top.sig.Mac.Lung[index,]

save(results.gsva,
     top.sig.Neutrophil.BM,
     top.sig.Neutrophil.BM.UP,
     top.sig.Monocyte.BM,
     top.sig.Monocyte.BM.UP,
     top.sig.Mac.PC,
     top.sig.Mac.PC.UP,
     top.sig.Kup.Liver,
     top.sig.Kup.Liver.UP,
     top.sig.Mac.Spleen,
     top.sig.Mac.Spleen.UP,
     top.sig.Mac.Colon,
     top.sig.Mac.Colon.UP,
     top.sig.Mac.Ileum,
     top.sig.Mac.Ileum.UP,
     top.sig.Mac.Lung,
     top.sig.Mac.Lung.UP,
     file="GSVA_All-KEGG-Rev_Sig-by-LIMMA_NSG_12Oct2017.RData")
rm(fit,design,index,contrast.matrix,allGeneSets)

######## Make an UpSet plot of the significant GSVA enrichments
list.KEGG.sig <- list(Neutrophil.BM.UP = c(rownames(top.sig.Neutrophil.BM.UP)),
                      Monocyte.BM.UP = c(rownames(top.sig.Monocyte.BM.UP)),
                      Mac.PC.UP = c(rownames(top.sig.Mac.PC.UP)),
                      Kup.Liver.UP = c(rownames(top.sig.Kup.Liver.UP)),
                      Mac.Spleen.UP = c(rownames(top.sig.Mac.Spleen.UP)),
                      Mac.Colon.UP = c(rownames(top.sig.Mac.Colon.UP)),
                      Mac.Ileum.UP = c(rownames(top.sig.Mac.Ileum.UP)),
                      Mac.Lung.UP = c(rownames(top.sig.Mac.Lung.UP)))
max.length <- max(sapply(list.KEGG.sig, length))
list.KEGG.sig <- lapply(list.KEGG.sig, function(v) { c(v, rep(NA, max.length-length(v)))})
list.KEGG.sig <- do.call(rbind, list.KEGG.sig)
list.KEGG.sig <- t(list.KEGG.sig)
list.KEGG.sig <- as.data.frame(list.KEGG.sig)
list.KEGG.sig[,1:length(colnames(list.KEGG.sig))] <- sapply(list.KEGG.sig[,1:length(colnames(list.KEGG.sig))],as.character)
list.KEGG.sig2 <- stack(list.KEGG.sig)
list.KEGG.sig2 <- list.KEGG.sig2[!is.na(list.KEGG.sig2$values),]
list.KEGG.sig2 <- list.KEGG.sig2[order(list.KEGG.sig2$values),]
list.KEGG.sig2 <- as.data.frame(unique(list.KEGG.sig2$values),stringsAsFactors=FALSE)
rownames(list.KEGG.sig2) <- list.KEGG.sig2[,1]
list.KEGG.sig2$Neutrophil.BM.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Neutrophil.BM.UP),nomatch = 0)
list.KEGG.sig2$Monocyte.BM.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Monocyte.BM.UP),nomatch = 0)
list.KEGG.sig2$Mac.PC.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.PC.UP),nomatch = 0)
list.KEGG.sig2$Kup.Liver.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Kup.Liver.UP),nomatch = 0)
list.KEGG.sig2$Mac.Spleen.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Spleen.UP),nomatch = 0)
list.KEGG.sig2$Mac.Colon.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Colon.UP),nomatch = 0)
list.KEGG.sig2$Mac.Ileum.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Ileum.UP),nomatch = 0)
list.KEGG.sig2$Mac.Lung.UP <- match(c(rownames(list.KEGG.sig2)),c(list.KEGG.sig$Mac.Lung.UP),nomatch = 0)
list.KEGG.sig2 <- list.KEGG.sig2[,2:length(colnames(list.KEGG.sig2))]
list.KEGG.sig2 <- apply(list.KEGG.sig2,1,function(x) {ifelse(x>0,1,0)})
list.KEGG.sig2 <- as.data.frame(t(list.KEGG.sig2))

library(UpSetR)
pdf("GSE63340_GSVA_KEGG-Rev_sig-UP_UpSet.pdf",height = 20,width = 100,onefile = FALSE)
upset(list.KEGG.sig2,nsets = length(colnames(list.KEGG.sig2)),nintersects = NA,sets = c(colnames(list.KEGG.sig2)),
      empty.intersections = "yes",order.by = "degree", show.numbers = "yes", group.by = "degree",
      mb.ratio = c(0.6, 0.4),number.angles = 0, point.size = 1, line.size = 0.5,
      mainbar.y.label = "Common Enriched KEGG Terms", sets.x.label = "Enriched KEGG Terms by Comparison", 
      text.scale = c(2, 1, 2, 1, 1.5, 0.5))
dev.off()

######## Make a heatmap for various significant comparisons of GSVA results
library(reshape)
library(xlsx)
wb <- createWorkbook()

sig_GS <- top.sig.Neutrophil.BM.UP
exprs_hm <- results.gsva[,c(1:2,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Neutrophil.BM")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Monocyte.BM.UP
exprs_hm <- results.gsva[,c(3:4,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Monocyte.BM")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.PC.UP
exprs_hm <- results.gsva[,c(7:8,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.PC")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Kup.Liver.UP
exprs_hm <- results.gsva[,c(9:10,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Kup.Liver")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Spleen.UP
exprs_hm <- results.gsva[,c(11:12,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Spleen")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Colon.UP
exprs_hm <- results.gsva[,c(13:14,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Colon")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Ileum.UP
exprs_hm <- results.gsva[,c(15:16,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Ileum")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)

sig_GS <- top.sig.Mac.Lung.UP
exprs_hm <- results.gsva[,c(17:18,5:6)] 
index <- which(row.names(exprs_hm)%in%row.names(sig_GS))
exprs_hm <- as.data.frame(exprs_hm[index,])
ws <- createSheet(wb=wb, sheetName = "Microglia_UP-over-Mac.Lung")
df <- transform(exprs_hm, KEGG = colsplit(rownames(exprs_hm), split = "\\;", names = c('ID', 'Term')))
df <- df[,c(5:6,1:4)]
df <- df[match(c(row.names(sig_GS)),row.names(df)),]
identical(rownames(sig_GS),rownames(df))
df <- cbind(df,sig_GS)
addDataFrame(x=df, sheet=ws, row.names=FALSE)


saveWorkbook(wb,file="GSE63340_GSVA-Rev_KEGG_sig-UP_NSG_12Oct2017.xlsx")
rm(wb,ws)

# index <- which(is.na(rowSums(exprs_hm)))
# if (length(index) > 0) { exprs_hm <- exprs_hm[-index,] }
# 
# distance <- dist(t(exprs_hm),method="euclidean")
# hc <- hclust(distance)
# 
# dist2 <- function(x, ...)
#   dist(x,method="euclidean")
# label.size = 1.2
# samples = colnames(exprs_hm)
# 
# hmcol <- colorRampPalette(c("blue","white","red"))(256)
# 
# hm <- d3heatmap(data.matrix(exprs_hm),
#                 colors = hmcol, scale = "row",
#                 dendrogram = "both",
#                 hclustfun = hclust,
#                 distfun = dist2,
#                 k_row = 10,
#                 k_col = 2,
#                 show_grid = 1,
#                 height = 4000,
#                 width = 1000,
#                 xaxis_font_size = "15pt",
#                 yaxis_font_size = "10pt",
#                 xaxis_height = 150,
#                 yaxis_width = 500)
# 
# saveWidget(hm,"Plus-vs-Minus_GSVA_KEGG_sig-UP.html")