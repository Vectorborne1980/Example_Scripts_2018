# load the esets
load("C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/eset/GSE32591_ESET-Affymetrix_Aesetgcrma_NSG_2016-07-11.RData")

# set working directories
dir <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/"
eset <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/eset/"
patient_data <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/ALL_data/patient_data/"

# load libraries
library(genefilter)
library(limma)
library(a4)
library(siggenes)
library(Biobase)

eset_all <- Aesetgcrma

#-- Patient data
setwd(patient_data)
pdata <- read.table("GSE32591_patient-data.txt", header=TRUE, sep="\t")
rownames(pdata) = pdata[,1]
pData = pdata[,-1]

pData$cohort <- factor(pData$cohort)
pData$cohort_full <- factor(pData$cohort_full)
pData$cohort <- gsub("healthy", "CTL", pData$cohort )
pData$cell_type <- factor(pData$cell_type)
pData$WHO_class <- factor(pData$WHO_class)

colnames(eset_all)
colnames(eset_all) <- gsub(".CEL","",colnames(eset_all)) # REMOVE .CEL FROM SAMPLE NAME
colnames(eset_all)
identical(rownames(pData),colnames(eset_all))

pData(eset_all) <- pData

table(eset_all$cohort)  # CTL 29 , LN 64
table(pData$WHO_class) # Class 2 = 16, Class 3 = 18, Class 4 = 26, Class 5 = 4
table(eset_all$cell_type) #Glomeruli = 46, Tubulointerstitium = 47

# Assign a friendly sample names
# NOTE: all samples are in the form GSM260XXX
eset_all$friendlyName <- paste( eset_all$cohort, substr( sampleNames(eset_all), 7,9),sep="" )
table(eset_all$friendlyName)

# ADD COHORT TO SAMPLE NAMES
sampleNames(eset_all)
cohort = as.character(eset_all$cohort)
sample_name = paste( cohort, sampleNames(eset_all), sep=".")
sampleNames(eset_all) = sample_name

rm(cohort, sample_name)

# create SYMBOL column required by some a4 functions
fData(eset_all)[,"SYMBOL"] <- as.character(fData(eset_all)[,"geneSymbol"])


######################## PCA ALL Data
library(affycoretools)
library(rgl)

eset_all$cohort <- factor(eset_all$cohort)
eset_all$cell_type <- factor(eset_all$cell_type)

# Plot PCA of cohort in 2D
groups <- as.numeric(eset_all$cohort)
groupnames = levels(eset_all$cohort)
labels <-  paste( eset_all$friendlyName, eset_all$WHO_class,sep="." )
plotPCA(eset_all, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: ALL_Data_by_Cohort",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17))

groups <- as.numeric(eset_all$cell_type)
groupnames = levels(eset_all$cell_type)
labels <-  paste( eset_all$friendlyName, eset_all$WHO_class,sep="." )
plotPCA(eset_all, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: ALL_Data_by_Tissue_type",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17) )

######################## OUTLIER REMOVAL
# REMOVE OUTLIERS IDENTIFIED USING PCA
# table(eset_all$pca_outlier)
# table(eset_all$pca_outlier_basis)
# index <- which( eset_all$pca_outlier==TRUE)
# eset_all <- eset_all[ , -index ]

# REMOVE OUTLIER FROM NUSE PLOT OF RAW AFFY DATA (POOR RNA HYBRIDIZATION)
# index <- which(colnames(eset_all)=="GSM260892")
# eset_all <- eset_all[ , -index ]

######################## Sub-set data by tissue type
eset.glom <- eset_all[,eset_all$cell_type=="Glomeruli"]
eset.tub <- eset_all[,eset_all$cell_type=="Tubulointerstitium"]

######################## PCA Tissue Sub-sets
groups <- as.factor(eset.glom$cohort)
groupnames = levels(eset.glom$cohort)
labels <-  paste( eset.glom$friendlyName, eset.glom$WHO_class, sep="." )
plotPCA(eset.glom, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: Glomeruli_ONLY",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )

groups <- as.factor(eset.tub$cohort)
groupnames = levels(eset.tub$cohort)
labels <-  paste( eset.tub$friendlyName, eset.glom$WHO_class, sep="." )
plotPCA(eset.tub, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: Tubulointersitium_ONLY",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )



#---------- Plot intensity histogram

min.probes <- round((min(rowMeans(exprs(eset.glom)))), digits = 2)  #2.17
max.probes <- round((max(rowMeans(exprs(eset.glom)))), digits = 2)  #15.36
mean.probes <- round((mean(rowMeans(exprs(eset.glom)))), digits = 2) #5.18

dev.off()
par(mfrow=c(1,2))

nrow <- format(as.numeric(nrow(eset.glom)),big.mark=",",scientific=F)

hist(rowMeans(exprs(eset.glom)), breaks=100, xlab="Average log2 Probe Intensity",
     main="GSE32591_ESET_GCRMA_Affy: Glomeruli, ALL Probes",
     sub=paste(nrow, " total probes, MIN ",min.probes,", MAX ",max.probes, ", MEAN ",mean.probes,sep=""), xlim=c(2,16), col="green")

hist(rowMeans(exprs(eset.glom)), breaks=5000, xlab="Average log2 Probe Intensity",
     main="GSE32591_ESET_GCRMA_Affy: Glomeruli, Probe Itensities 2.15 - 2.25",
     sub=paste(nrow, " total probes, MIN ",min.probes,", MAX ",max.probes, ", MEAN ",mean.probes,sep=""), xlim=c(2.15,2.25), col="green")

abline(v=2.22,col="red",lty=2)

# # Many of the genes on the chip won't be expressed, or might have only small variability across the samples.
# # We try to remove these genes' corresponding probe sets with an intensity filter (the intensity of a gene should be above log2(100)
# # in at least 25% of the samples), and a variance filter (the interquartile range of log2-intensities should be at least 0.5).
# # We create a new exprSet containing only the probe sets which passed our filter. 
# library(genefilter)
# f1 <- pOverA(0.25,log2(100))
# f2 <- function(x) (IQR(x)>0.5)
# ff <- filterfun(f1,f2)
# selected <- genefilter (eset.glom,ff)
# sum(selected) #5596 retained probes
# eset.glom.genefiltered <- eset.glom[selected,]
# setwd(eset)
# save(eset.glom.genefiltered, file="GSE32591_ESET-GCRMA-Affymetrix_Glomeruli-GeneFiltered-All-Classes_NSG_2016-07-11.Rdata")
# rm(selected,f1,f2,ff)

######################## Sub-set data by tissue type
# eset.glom.LN=eset.glom[,eset.glom$cohort=="LN"]
# eset.glom.CTL=eset.glom[,eset.glom$cohort=="CTL"]
# index1 <- which(rowMeans(exprs(eset.glom.LN))<2.22)
# index2 <- which(rowMeans(exprs(eset.glom.CTL))<2.22)

######################## LOW INTENSITY PROBE FILTERING 
index <- which(rowMeans(exprs(eset.tub))<2.22)
nrow <- nrow(data.frame(index)) # TI: 22215 probes - 4790 probes = 17425 probes remaining; Glom: 22215 - 6321 = 15894
eset.tub.filtered = eset.tub[-index,]
setwd(eset)
save(eset.tub.filtered, file="GSE32591_ESET-GCRMA-Affymetrix_Tubulointerstitium-Filtered-All-Classes_NSG_2016-07-11.Rdata")
rm(index, nrow)

# Examine if CTL vs. SLE distributions are different from each other
capture.output(write.table(sapply(levels(eset.tub.filtered$cohort), function(x) summary(as.vector(exprs(eset.tub.filtered)[,+  which(eset.tub.filtered$cohort == x)])))), file = "GSE32591_LIMMA-Affymetrix_Tubulointerstitium-Filtered-All-Classes_Statistics_NSG_2016-07-11.txt")
# sapply(levels(eset.glom.genefiltered$cohort), function(x) summary(as.vector(exprs(eset.glom.genefiltered)[,+  which(eset.glom.genefiltered$cohort == x)])))

######################## Sub-set data by tissue type
# eset.glom.filtered=eset_all.filtered[,eset_all.filtered$cell_type=="Glomeruli"]
# eset.tub.filtered=eset_all.filtered[,eset_all.filtered$cell_type=="Tubulointerstitium"]

######################## PCA Tissue Sub-sets
dev.off()
groups <- as.factor(eset.glom.filtered$cohort)
groupnames = levels(eset.glom.filtered$cohort)
labels <-  paste( eset.glom.filtered$friendlyName, eset.glom$WHO_class, sep="." )
plotPCA(eset.glom.filtered, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: Glomeruli_ONLY",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )

groups <- as.factor(eset.tub.filtered$cohort)
groupnames = levels(eset.tub.filtered$cohort)
labels <-  paste( eset.tub.filtered$friendlyName, eset.glom$WHO_class, sep="." )
plotPCA(eset.tub.filtered, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy: Tubulointersitium_ONLY",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )
# Filtered and all probes PCAs look almost identical, ergo low intensity probes played very little role in sample segregation (?)


######################## RUN LIMMA ON TISSUE SUBSETS

library(limma)

design <- model.matrix(~0 + eset.tub.filtered$cohort)
colnames(design) <- levels(factor(eset.tub.filtered$cohort))
fit <- lmFit(eset.tub.filtered, design)
names(fit)
contrast.matrix <- makeContrasts(LN-CTL,levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top.limma <- topTable(fit2, coef=1, adjust="BH", n=nrow(eset.tub.filtered))
top.limma$probe <- rownames(top.limma)
colnames(top.limma)
index <- which(top.limma$adj.P.Val<0.2)
top.limma.FDR <- top.limma[-index,]
setwd(eset)
save(top.limma.FDR,file="GSE32591_LIMMA-Affymetrix_Tubulointerstitium-Filtered-All-Classes_Top-limma_NSG_2016-07-11.Rdata")
#write.table(top.limma, file = "/Users/Chinnu/Desktop/GSE10325_inactives/Figures/topsiggenes_inactives_affy1.txt", quote=FALSE, sep = "\t",row.names=FALSE, col.names=TRUE)



########################  SUBSET ESET BASED WHO CLASS TYPE 
# 
# table(eset_all$WHO_class)  # Class2- 8LN, Class3&4- 22LN, Class5 - 2LN
# 
# # -- remove class 3 and 4
# index <- which(eset_all$WHO_class=="3")
# eset <- eset_all[,-index] # removed class 3
# index <- which(eset$WHO_class=="4") 
# eseta <- eset[,-index] #removed class 4
# 
# #--- Class 2
# index <- which(eseta$WHO_class=="5")
# eset2 <- eseta[,-index] #removed class 5
# table(eset2$cohort)  # CTL 15, LN 8
# 
# #--- Class 2: without any class type
# index <- which(eset2$WHO_class_type=="b")
# eset2_noclasstype <- eset2[,-index] #removed class type b from eset2
# table(eset2_noclasstype$cohort)
# 
# #--- Class 2b
# index <- which(eset2$WHO_class_type=="na")
# eset2_type2b <- eset2[,-index] #removed class type b from eset2
# table(eset2_type2b$cohort)
# 
# # ---- Save the esets
# setwd(eset)
# save(eset2_type2b,file="eset2b_affy.RData")
# save(eset2_noclasstype,file="eset2_noclasstype_affy.RData")


# Classes 3 & 4: without any class type
index = which(eset.glom.filtered$WHO_class==3)
index = c(index,which(eset.glom.filtered$WHO_class==4))
index = c(index,which(eset.glom.filtered$cohort=="CTL"))
eset.glom.filtered.3.4 = eset.glom.filtered[,index]

labels <-  paste( eset.glom.filtered.3.4$friendlyName, eset.glom.filtered.3.4$WHO_class, eset.glom.filtered.3.4$WHO_class_type, sep="." )

groups <- as.factor(eset.glom.filtered.3.4$cohort)
groupnames = levels(eset.glom.filtered.3.4$cohort)
plotPCA(eset.glom.filtered.3.4, groups=groups, groupnames=groupnames, main="GSE32591_ESET_GCRMA_Affy_Filtered: Glomeruli ONLY, Classes 3&4",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )

sampleNames(eset.glom.filtered.3.4)

setwd(eset)
save(eset.glom.filtered.3.4, file="GSE32591_ESET-GCRMA-Affymetrix_Glomeruli-Classes-3-4-Filtered_NSG_2016-07-11.RData")
######################## ANALYSES OF CLASS TYPES

# ---------- Plot Dendogram
distance <- dist(t(exprs(eset.tub.3.4)),method="euclidean")
clusters <- hclust(distance)
labels <- paste( eset.tub.3.4$friendlyName, eset.tub.3.4$WHO_class, eset.tub.3.4$WHO_class_type, sep="." )
labels <- gsub( ".control", "", labels)
plot(clusters, cex=1, main="GSE32591 WHO_class3.4, Affy GCRMA,Euclidean distance" ,lwd=1, labels=labels)

# ---------- plot PCA
# library(affycoretools)
# library(rgl)
# 
# eset_all$cohort <- factor(eset_all$cohort)
# labels <-  paste( eset2_noclasstype$friendlyName, eset2_noclasstype$WHO_class, eset2_noclasstype$WHO_class_type, sep="." )

# # Plot PCA of cohort in 2D
# groups <- as.numeric(eset2_noclasstype$cohort)
# groupnames = levels(eset2_noclasstype$cohort)
# plotPCA(eset2_noclasstype, groups=groups, groupnames=groupnames, main="GSE32591 WHO_class:No Class Type Specified, Affy GCRMA PCA",
#         pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )









# ------------------------- GENE FILTERING --------------------------------------------------------------

#------------ Filter genes without gene symbol
# We filter out probes with empty gene names, a type of annotation based filtering
index <- which(is.na(fData(eset2_noclasstype)[,"geneSymbol"])) # 2,269(affy gcrma) probes removed without gene mapping
eset2_noclasstype_noNAS <- eset2_noclasstype[-index,]
eset2_noclasstype_filtered <- eset2_noclasstype[-index,]


#------------ Filter low intensity probes
# Filter out very low intensity probes. We use a threshold of 2.82, which is the highest we can go before
# we lose informative genes (as we found in the accompanying script find_lowpass_threshold.R)
index <- which(rowMeans(exprs(eset2_noclasstype_filtered)) < 2.82) # 7658 "lowpass" probes
eset2_noclasstype_filtered <- eset2_noclasstype_filtered[-index,]
rm(index)


#------------ ---LIMMA Analysis ---------------------------------------------------------------------
# Run LIMMA on LN patients of WHO_CLASS: 2
library(a4)

limmaResult <- limmaTwoLevels(eset2_noclasstype_filtered, "cohort")
propDEgenes(limmaResult) # 48.4% genes are differentially expressed

top.active <- topTable(limmaResult, n = nrow(eset2_noclasstype_filtered), coef = 2)

# plot LIMMA results
dev.off()
volcanoPlot(limmaResult, main="GSE32591 WHO_class2: No type, Affy CDF GCRMA Glomeruli LN and normal LIMMA",
            xlab="Log2 Normalized Probe Intensity. x<0 is downregulated, x>0 is upregulated.",
            ylab="Adjusted p-val of differential significance")

# # RUN SAM ON CLASS2_no class type LN PATIENTS AND MERGE WITH LIMMA RESULTS (NO LFC OR FDR FILTERING YET)
# library(siggenes)
# # assign numeric cohort designation
# SAM.cl = as.character(eset2_noclasstype_filtered$cohort)
# SAM.cl[which(SAM.cl=="CTL")]=0
# SAM.cl[which(SAM.cl=="LN")]=1
# SAM.cl = as.numeric(SAM.cl)
# SAM.cl
# # run SAM
# SAM.out <- sam(exprs(eset2_noclasstype_filtered), SAM.cl, rand = 1234567 )
# rm(SAM.cl)
# SAM.out
# summary(SAM.out)
# # We set the FDR to 0.2 for consistency with previous experiments
# FDR <- findDelta(SAM.out, fdr = 0.2)
# FDR <- FDR[1]
# dev.off()
# plot(SAM.out, FDR)
# #identify(SAM.out) # click on the plot to identify genes
# SAM.summary <- summary(SAM.out, FDR)
# #SAM.summary@row.sig.genes
# #SAM.summary@mat.sig
# 
# SAM <- data.frame(SAM.summary@mat.sig)
# SAM$SAM_order <- rep(1:nrow(SAM))
# 
# # merge SAM and LIMMA results
# SAM <- cbind( SAM_probe=rownames(SAM ), SAM )
# 
# 
# top_merged = merge(top.active, SAM, by.x="probe", by.y="SAM_probe", all.x=TRUE )
# 
# top_merged <- top_merged[,-6]
# 
# colnames(top_merged)
# colnames(top_merged)= c("probe","geneSymbol","geneName","geneEntrezID","geneEnsembl","LIMMA.logFC",
#                         "LIMMA.AveExpr", "LIMMA.t","LIMMA.pVal", "LIMMA.pVal.adj", "LIMMA.B",
#                         "SAM.row", "SAM.dVal", "SAM.stdev", "SAM.rawP", "SAM.qVal", "SAM.Rfold", "SAM_order")
# #top_merged = top_merged[, c(1,2,3,4,6,7,8,9,10,11,12,13,17,18) ]
# top_merged = top_merged[ order(top_merged$LIMMA.pVal.adj), ]
# 
# # assign up or down
# top_merged$fc_direction <- "UP"
# index <- which(top_merged$LIMMA.logFC<=0)
# top_merged[ index, "fc_direction"] <- "DOWN"
# 
# # indicate CDF used
# top_merged$CDF <- "Affy_HGU133A"
# #top_merged$CDF <- "BrainArray_hgu133ahsentrezgcdf"
# 
# GSE32591_Affy_GCRMA_Active <- top_merged
# 
# ## SAVE ACTIVE LUPUS VS. NORMALS MERGED LIMMA/SAM RESULTS
# setwd(eset)
# write.table(GSE32591_Affy_GCRMA_Active ,file="GSE32591_Class2_no class type_Affy_GLOM_LIMMA_SAM_ACTIVE_All.txt",sep="\t",quote=FALSE,row.names=FALSE)
# save(GSE32591_Affy_GCRMA_Active, file="GSE32591_Class2_no class type_Affy_GLOM_LIMMA_SAM_ACTIVE_All.RData")

#--------------------- HEATMAP : CLASS 2 - not class type specified-------------------------------------------------------------------
library(gplots)

# First subset out significant probes
index <- which(top.limma$adj.P.Val < 0.2) 

# top.limma.sig <- top.limma[ index, ]
top.limma.sig <- top.limma[ c(1:1000), ]

# Use the rownames of top.limma.sig as a lookup column for the eset
top <- rownames(top.limma.sig)
index <- which(rownames(eset.glom.3.4.filtered)%in%top)
top <- exprs(eset.glom.3.4.filtered)[ index, ]
m.use = as.matrix(top)
rownames(m.use) <- fData(eset.glom.3.4.filtered)[index,"geneSymbol"]
colnames(m.use) = eset.glom.3.4.filtered$friendlyName
colnames(m.use) = paste( eset.glom.3.4.filtered$friendlyName, eset.glom.3.4.filtered$WHO_class, eset.glom.3.4.filtered$WHO_class_type, sep="." )
m.use <- m.use[ , order(colnames(m.use))]

# Need to match how the array dendrogram was created.
hc = hclust(dist(t(m.use)),method = "complete") 

distance <- dist(t(m.use),method="euclidean")
hc <- hclust(distance)

hc.cut = cutree(hc,k=6) # k = number of clusters
table(hc.cut) # show cluster cut results
hc.cut.levels = levels(as.factor(hc.cut)) # factor cluster level numbers for coloring
hc.cut.levels # review cluster level factors

colors = rainbow(length(hc.cut.levels));colors
cluster.patient.colors = c()
for (i in 1:length(hc.cut.levels)){
  cluster.patient.colors[names(which(hc.cut==hc.cut.levels[i]))] = colors[i] }
cluster.patient.colors

dist2 = dist2 <- function(x, ...)
  dist(x,method="euclidean")
label.size = 1
samples = colnames(m.use)
hmcol <- colorRampPalette(c("blue","white","red"))(256)

## MAKE SURE TO MAKE THE R STUDIO PLOTS FRAME AS LARGE AS POSSIBLE BEFORE RUNNING (SHRINK THE OTHER FRAMES)
dev.off()
genes.heat = heatmap.2(data.matrix(m.use), col=hmcol, scale="row", distfun = dist2, key=TRUE,
                       symkey=TRUE, density.info="none", trace="none",
                       ColSideColors=cluster.patient.colors[samples],cexCol=label.size,
                       cexRow = label.size,margins=c(9,9),
                       main="GSE32591 Class 3&4, Affy CDF GCRMA Glomeruli LN-CTL LIMMA, top 1k")

# EXAMINE CLUSTERS OF PROBES
hc_probes = hclust(dist(m.use),method = "complete")
hc_probes.cut = cutree(hc_probes,k=8) # k = number of clusters
table(hc_probes.cut) # show cluster cut results

# Now extract probes within clusters
clusters <- data.frame(cluster=hc_probes.cut,probe=names(hc_probes.cut))
clusters <- clusters[ order(clusters$cluster, clusters$probe), ]

# Save the clustered probes
GSE32591_Class2_notype_GLOM_GCRMA_Active_clusters <- clusters
setwd(eset)
save(GSE32591_Class2_notype_GLOM_GCRMA_Active_clusters,file="GSE32591_Class2_notype_Affy_GLOM_GCRMA_Active_clusters.RData")
