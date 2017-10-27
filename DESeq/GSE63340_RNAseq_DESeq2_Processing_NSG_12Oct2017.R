#Note: importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(vsn)
library(biomaRt)
library(Biobase)
library(limma)
library(edgeR)
setwd("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/")

###Import the count data
gse63340.deseq <- readRDS("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/gse63340.rds")
dds <- gse63340.deseq
colnames(dds) <- c(paste(dds$Condition,c(rep(c(1,2),9)),sep = "."))
# gse6340.SE <- SummarizedExperiment(gse63340.deseq)
# emf.deseq <- readRDS("C:/Users/Nick Geraci/Nick/UVA/Sandbox/gse63340/Bulk_RNA-seq/emf.rds")
# emf.SE <- SummarizedExperiment(emf.deseq)
# # If the SummarizedExperiment objects have the same rowData you can use cbind() to 
# # combine the data into a single SE object. This combines the data in info() and assays()
# # so the counts of each SE object end up in the assays()$counts matrix as separate columns.
# # In this way you can create one SE object that holds all counts for the different experiments.
# # cmb <- cbind(SE1, SE2) 
# index <- which(names(gse6340.SE)%in%names(emf.SE))
# gse6340.SE <- gse6340.SE[index,]
# index <- which(names(emf.SE)%in%names(gse6340.SE))
# emf.SE <- emf.SE[index,]
# gse6340.SE <- gse6340.SE[order(names(gse6340.SE))]
# emf.SE <- emf.SE[order(names(emf.SE)),]
# identical(names(gse6340.SE),names(emf.SE))
# gse63340.counts <- assay(gse6340.SE)
# gse63340.emf.SE.combo <- cbind(gse6340.SE,emf.SE)
# ###Run DESeq2
# #Set up the DESeqDataSet Object and run the DESeq pipeline
# dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metastuff, design = ~Condition)
dds <- DESeq(dds)
pdf(file="Dispersion_Estimates.pdf",height = 9,width = 16)
plotDispEsts( dds, ylim = c(1e-8, 1e1) ) #Check dispersion estimates 
dev.off()

###Differential Analyses
# Compare the Conditions for "Neutrophil_Bone_Marrow" over "Microglia_Brain"
res.Neutrophil_Bone_Marrow.Microglia_Brain <- results(dds, contrast=c("Condition", "Neutrophil_Bone_Marrow", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Neutrophil_Bone_Marrow.Microglia_Brain <- res.Neutrophil_Bone_Marrow.Microglia_Brain[order(res.Neutrophil_Bone_Marrow.Microglia_Brain$padj),]
summary(res.Neutrophil_Bone_Marrow.Microglia_Brain)
index <- which(res.Neutrophil_Bone_Marrow.Microglia_Brain$padj<0.05)
res.Neutrophil_Bone_Marrow.Microglia_Brain.sig <- res.Neutrophil_Bone_Marrow.Microglia_Brain[index,]
index <- which(res.Neutrophil_Bone_Marrow.Microglia_Brain.sig$log2FoldChange>0)
res.Neutrophil_Bone_Marrow.Microglia_Brain.sig.UP <- res.Neutrophil_Bone_Marrow.Microglia_Brain.sig[index,]
# Compare the Conditions for "Monocyte_Bone_Marrow" over "Microglia_Brain"
res.Monocyte_Bone_Marrow.Microglia_Brain <- results(dds, contrast=c("Condition", "Monocyte_Bone_Marrow", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Monocyte_Bone_Marrow.Microglia_Brain <- res.Monocyte_Bone_Marrow.Microglia_Brain[order(res.Monocyte_Bone_Marrow.Microglia_Brain$padj),]
summary(res.Monocyte_Bone_Marrow.Microglia_Brain)
index <- which(res.Monocyte_Bone_Marrow.Microglia_Brain$padj<0.05)
res.Monocyte_Bone_Marrow.Microglia_Brain.sig <- res.Monocyte_Bone_Marrow.Microglia_Brain[index,]
index <- which(res.Monocyte_Bone_Marrow.Microglia_Brain.sig$log2FoldChange>0)
res.Monocyte_Bone_Marrow.Microglia_Brain.sig.UP <- res.Monocyte_Bone_Marrow.Microglia_Brain.sig[index,]
# Compare the Conditions for "Mac_Peritoneal_Cavity" over "Microglia_Brain"
res.Mac_Peritoneal_Cavity.Microglia_Brain <- results(dds, contrast=c("Condition", "Mac_Peritoneal_Cavity", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Mac_Peritoneal_Cavity.Microglia_Brain <- res.Mac_Peritoneal_Cavity.Microglia_Brain[order(res.Mac_Peritoneal_Cavity.Microglia_Brain$padj),]
summary(res.Mac_Peritoneal_Cavity.Microglia_Brain)
index <- which(res.Mac_Peritoneal_Cavity.Microglia_Brain$padj<0.05)
res.Mac_Peritoneal_Cavity.Microglia_Brain.sig <- res.Mac_Peritoneal_Cavity.Microglia_Brain[index,]
index <- which(res.Mac_Peritoneal_Cavity.Microglia_Brain.sig$log2FoldChange>0)
res.Mac_Peritoneal_Cavity.Microglia_Brain.sig.UP <- res.Mac_Peritoneal_Cavity.Microglia_Brain.sig[index,]
# Compare the Conditions for "Mac_Spleen" over "Microglia_Brain"
res.Mac_Spleen.Microglia_Brain <- results(dds, contrast=c("Condition", "Mac_Spleen", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Mac_Spleen.Microglia_Brain <- res.Mac_Spleen.Microglia_Brain[order(res.Mac_Spleen.Microglia_Brain$padj),]
summary(res.Mac_Spleen.Microglia_Brain)
index <- which(res.Mac_Spleen.Microglia_Brain$padj<0.05)
res.Mac_Spleen.Microglia_Brain.sig <- res.Mac_Spleen.Microglia_Brain[index,]
index <- which(res.Mac_Spleen.Microglia_Brain.sig$log2FoldChange>0)
res.Mac_Spleen.Microglia_Brain.sig.UP <- res.Mac_Spleen.Microglia_Brain.sig[index,]
# Compare the Conditions for "Mac_Colon" over "Microglia_Brain"
res.Mac_Colon.Microglia_Brain <- results(dds, contrast=c("Condition", "Mac_Colon", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Mac_Colon.Microglia_Brain <- res.Mac_Colon.Microglia_Brain[order(res.Mac_Colon.Microglia_Brain$padj),]
summary(res.Mac_Colon.Microglia_Brain)
index <- which(res.Mac_Colon.Microglia_Brain$padj<0.05)
res.Mac_Colon.Microglia_Brain.sig <- res.Mac_Colon.Microglia_Brain[index,]
index <- which(res.Mac_Colon.Microglia_Brain.sig$log2FoldChange>0)
res.Mac_Colon.Microglia_Brain.sig.UP <- res.Mac_Colon.Microglia_Brain.sig[index,]
# Compare the Conditions for "Mac_Ileum" over "Microglia_Brain"
res.Mac_Ileum.Microglia_Brain <- results(dds, contrast=c("Condition", "Mac_Ileum", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Mac_Ileum.Microglia_Brain <- res.Mac_Ileum.Microglia_Brain[order(res.Mac_Ileum.Microglia_Brain$padj),]
summary(res.Mac_Ileum.Microglia_Brain)
index <- which(res.Mac_Ileum.Microglia_Brain$padj<0.05)
res.Mac_Ileum.Microglia_Brain.sig <- res.Mac_Ileum.Microglia_Brain[index,]
index <- which(res.Mac_Ileum.Microglia_Brain.sig$log2FoldChange>0)
res.Mac_Ileum.Microglia_Brain.sig.UP <- res.Mac_Ileum.Microglia_Brain.sig[index,]
# Compare the Conditions for "Mac_Lung" over "Microglia_Brain"
res.Mac_Lung.Microglia_Brain <- results(dds, contrast=c("Condition", "Mac_Lung", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Mac_Lung.Microglia_Brain <- res.Mac_Lung.Microglia_Brain[order(res.Mac_Lung.Microglia_Brain$padj),]
summary(res.Mac_Lung.Microglia_Brain)
index <- which(res.Mac_Lung.Microglia_Brain$padj<0.05)
res.Mac_Lung.Microglia_Brain.sig <- res.Mac_Lung.Microglia_Brain[index,]
index <- which(res.Mac_Lung.Microglia_Brain.sig$log2FoldChange>0)
res.Mac_Lung.Microglia_Brain.sig.UP <- res.Mac_Lung.Microglia_Brain.sig[index,]
# Compare the Conditions for "Kupffer_Liver" over "Microglia_Brain"
res.Kupffer_Liver.Microglia_Brain <- results(dds, contrast=c("Condition", "Kupffer_Liver", "Microglia_Brain"), alpha = 0.05, pAdjustMethod = "BH")
res.Kupffer_Liver.Microglia_Brain <- res.Kupffer_Liver.Microglia_Brain[order(res.Kupffer_Liver.Microglia_Brain$padj),]
summary(res.Kupffer_Liver.Microglia_Brain)
index <- which(res.Kupffer_Liver.Microglia_Brain$padj<0.05)
res.Kupffer_Liver.Microglia_Brain.sig <- res.Kupffer_Liver.Microglia_Brain[index,]
index <- which(res.Kupffer_Liver.Microglia_Brain.sig$log2FoldChange>0)
res.Kupffer_Liver.Microglia_Brain.sig.UP <- res.Kupffer_Liver.Microglia_Brain.sig[index,]


# Check MA plots of differential data
pdf(file="MA_plots.pdf",height = 9,width = 16)
DESeq2::plotMA( dds, ylim = c(-10, 10), main = "Default DESeq2")
DESeq2::plotMA( res.Neutrophil_Bone_Marrow.Microglia_Brain, ylim = c(-10, 10), main = "Neutrophil_Bone_Marrow vs Microglia_Brain" )
DESeq2::plotMA( res.Monocyte_Bone_Marrow.Microglia_Brain, ylim = c(-10, 10), main = "Monocyte_Bone_Marrow vs Microglia_Brain" )
DESeq2::plotMA( res.Mac_Peritoneal_Cavity.Microglia_Brain, ylim = c(-10, 10), main = "Mac_Peritoneal_Cavity vs Microglia_Brain" )
DESeq2::plotMA( res.Mac_Spleen.Microglia_Brain, ylim = c(-10, 10), main = "Mac_Spleen vs Microglia_Brain" )
DESeq2::plotMA( res.Mac_Colon.Microglia_Brain, ylim = c(-10, 10), main = "Mac_Colon vs Microglia_Brain" )
DESeq2::plotMA( res.Mac_Ileum.Microglia_Brain, ylim = c(-10, 10), main = "Mac_Ileum vs Microglia_Brain" )
DESeq2::plotMA( res.Mac_Lung.Microglia_Brain, ylim = c(-10, 10), main = "Mac_Lung vs Microglia_Brain" )
DESeq2::plotMA( res.Kupffer_Liver.Microglia_Brain, ylim = c(-10, 10), main = "Kupffer_Liver vs Microglia_Brain" )
dev.off()

pdf(file="MA_plots_Significant.pdf",height = 9,width = 16)
DESeq2::plotMA( res.Neutrophil_Bone_Marrow.Microglia_Brain.sig, ylim = c(-10, 10), main = "Neutrophil_Bone_Marrow vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Monocyte_Bone_Marrow.Microglia_Brain.sig, ylim = c(-10, 10), main = "Monocyte_Bone_Marrow vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Mac_Peritoneal_Cavity.Microglia_Brain.sig, ylim = c(-10, 10), main = "Mac_Peritoneal_Cavity vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Mac_Spleen.Microglia_Brain.sig, ylim = c(-10, 10), main = "Mac_Spleen vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Mac_Colon.Microglia_Brain.sig, ylim = c(-10, 10), main = "Mac_Colon vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Mac_Ileum.Microglia_Brain.sig, ylim = c(-10, 10), main = "Mac_Ileum vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Mac_Lung.Microglia_Brain.sig, ylim = c(-10, 10), main = "Mac_Lung vs Microglia_Brain, padj < 0.05" )
DESeq2::plotMA( res.Kupffer_Liver.Microglia_Brain.sig, ylim = c(-10, 10), main = "Kupffer_Liver vs Microglia_Brain, padj < 0.05" )

dev.off()

###Data Transformation
# Do Log2 transformation
dds.normT <- normTransform(dds)
dds.counts <- assay(dds.normT)
# # Alternatively:
# dds.counts <- log2(1+counts(dds,normalized=TRUE))

# Do regularized log transformation
rld <- rlogTransformation(dds, blind=TRUE)
rld.counts <- assay(rld)
rld.counts <- as.data.frame(rld.counts)

# Do variance stablizing transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsd.counts <- assay(vsd)
vsd.counts <- as.data.frame(vsd.counts)

## Compare Data Transformation Methods
pdf(file="Transformation_Comparison_1.pdf",height=9,width=16)
par(mfrow=c(1,3))
plot(assay(dds.normT)[,1:4],col="#00000020",pch=20,cex=0.3,main="Log2 Transformed")
plot(assay(rld)[,1:4],col="#00000020",pch=20,cex=0.3,main="Regularized Logarithm (rLog)")
plot(assay(vsd)[,1:4],col="#00000020",pch=20,cex=0.3,main="Variance Stabilization Tranformation (VST)")
dev.off()

pdf(file="Transformation_Comparison_2.pdf",height = 9,width = 16)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
colors <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Log2 Transformed")
heatmap.2(assay(rld)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Regularized Logarithm (rLog)")
heatmap.2(assay(vsd)[select,], col = colors,
          Rowv = FALSE, Colv = FALSE, scale='none',
          dendrogram='none', trace='none', margin=c(9,16), main = "Variance Stabilization Tranformation (VST)")
dev.off()

pdf(file="Transformation_Comparison_3.pdf",height=9,width=16)
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(assay(dds.normT[notAllZero,]))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

# Show flattening of data toward consistent variance
pdf(file="Variance_Flattening.pdf",height = 9,width = 16)
par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'red')
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="p", lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

# Consider the total distance between samples
pdf(file="Sample_Distances_rLog.pdf",height = 9,width = 16)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Condition
colnames(sampleDistMatrix) <- rld$Condition
hc <- hclust(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=colors,margins=c(9,16),
          main = "Regularized Logarithm (rLog)")
dev.off()

# Display PCA of rLog Data
pdf("PCA (rLog).pdf",height = 9, width = 16)
print(plotPCA( rld, intgroup = c("Condition"), colors(distinct = FALSE)))
dev.off()

pdf(file="Sample_Distances_VST.pdf",height = 9,width = 16)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Condition
colnames(sampleDistMatrix) <- vsd$Condition
hc <- hclust(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
heatmap.2(sampleDistMatrix,Rowv=as.dendrogram(hc),symm=TRUE,trace="none",col=colors,margins=c(9,16),
          main = "Variance Stabilization Tranformation (VST)")
dev.off()
rm(colors,last,px,ord,index,hc,sampleDistMatrix,sampleDists,vstcol,select,notAllZero)

### rLog is very obviously the better method for count transformations

# Display the 35 most variant genes of rLog data
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
pdf("Most_Variant_rLog_Genes.pdf",height = 9,width = 16)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="both",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           margins = c(16,9))
dev.off()

# Add additional gene IDs to the normalized results
rld.counts$ensembl <- sapply( strsplit( rownames(rld.counts), split="nn+" ), "[", 1 )
ensembl <- useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene"),
                  filters = "ensembl_gene_id",
                  values = rld.counts$ensembl,
                  mart = ensembl )
idx <- match( rld.counts$ensembl, genemap$ensembl_gene_id )
rld.counts$entrez <- genemap$entrezgene[ idx ]
rld.counts <- rld.counts[!duplicated(rld.counts$entrez), ]
rld.counts <- rld.counts[!is.na(rld.counts$entrez), ]
row.names(rld.counts) <- rld.counts$entrez
rld.counts <- rld.counts[,1:18]
write.table(rld.counts,file="Exprs_Data_AdjP-0-05_rLog_NSG_12Oct2017.txt",row.names = TRUE,col.names = TRUE,sep = "\t")
rm(idx,ensembl,genemap)
###Save results
save(res.Neutrophil_Bone_Marrow.Microglia_Brain.sig,
     res.Neutrophil_Bone_Marrow.Microglia_Brain.sig.UP,
     res.Monocyte_Bone_Marrow.Microglia_Brain.sig,
     res.Monocyte_Bone_Marrow.Microglia_Brain.sig.UP,
     res.Mac_Peritoneal_Cavity.Microglia_Brain.sig,
     res.Mac_Peritoneal_Cavity.Microglia_Brain.sig.UP,
     res.Mac_Spleen.Microglia_Brain.sig,
     res.Mac_Spleen.Microglia_Brain.sig.UP,
     res.Mac_Colon.Microglia_Brain.sig,
     res.Mac_Colon.Microglia_Brain.sig.UP,
     res.Mac_Ileum.Microglia_Brain.sig,
     res.Mac_Ileum.Microglia_Brain.sig.UP,
     res.Mac_Lung.Microglia_Brain.sig,
     res.Mac_Lung.Microglia_Brain.sig.UP,
     res.Kupffer_Liver.Microglia_Brain.sig,
     res.Kupffer_Liver.Microglia_Brain.sig.UP,
     dds,
     dds.normT,
     rld,
     rld.counts,
     vsd,
     vsd.counts,
     file="GSE63340_RNAseq_DESeq2_Processed_NSG_12Oct2017.RData")
