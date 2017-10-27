# load the esets
load("C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/GLOM/eset/Besetgcrma.RData")
#load("C:/Users/Sulbha/Desktop/AMPEL/my work/myeloid/GSE32591/TUB/analysis class2 type b/eset/Besetgcrma.RData")

# set working directories
dir <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/GLOM"
eset <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/GLOM/eset"
patient_data <- "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/GLOM/patient data"

# load libraries
library(genefilter)
library(limma)
library(a4)
library(siggenes)

eset_all <- Besetgcrma

#-- Patient data
setwd(patient_data)
pdata <- read.table("patient data.txt", header=TRUE, sep="\t")
rownames(pdata) = pdata[,1]
pData = pdata[,-1]
pData$cohort <- factor(pData$cohort)
pData$cohort_full <- factor(pData$cohort_full)
pData$cohort <- gsub("healthy", "CTL", pData$cohort )
pData$cell_type <- factor(pData$cell_type)
colnames(eset_all)
colnames(eset_all) = sub(".CEL","",colnames(eset_all)) # REMOVE .CEL FROM SAMPLE NAME
colnames(eset_all)
identical(rownames(pData),colnames(eset_all))
pData(eset_all) <- pData
table(eset_all$cohort)  # CTL 29 , LN 64
pData$WHO_class <- factor(pData$WHO_class)
table(pData$WHO_class)
table(eset_all$cell_type) # Glom 46 Tub 47

# REMOVE Affy PROBES
index <- which( substr(rownames(eset_all),1,4)=="AFFX" )
eset_all <- eset_all[ -index, ]
rm(index)

# REMOVE .CEL FROM SAMPLE NAME
# sampleNames(eset_all) <- substr( sampleNames(eset_all), 1,10)

# Assign a friendly sample name. Note all samples are in the form GSM260XXX
eset_all$friendlyName <- paste( eset_all$cohort, substr( sampleNames(eset_all), 7,9),sep="" )

# ADD COHORT TO SAMPLE NAMES
sampleNames(eset_all)
cohort = as.character(eset_all$cohort)
cohort <- gsub("healthy", "CTL", cohort)
sample_name = paste( cohort, sampleNames(eset_all), sep=".")
sampleNames(eset_all) = sample_name

rm(cohort, sample_name)

# create SYMBOL column required by some a4 functions
fData(eset_all)[,"SYMBOL"] <- as.character(fData(eset_all)[,"geneSymbol"])

#PCA
library(affycoretools)
library(rgl)

eset_all$cohort <- factor(eset_all$cohort)
labels <-  paste( eset_all$friendlyName, eset_all$WHO_class, eset_all$WHO_class_type, sep="." )

# Plot PCA of cohort in 2D
groups <- as.numeric(eset_all$cohort)
groupnames = levels(eset_all$cohort)
plotPCA(eset_all, groups=groups, groupnames=groupnames, main="GSE32591 WHO_class:No Class Type Specified, Affy GCRMA PCA",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )


groups <- as.factor(eset_all$cell_type)
groupnames = levels(eset_all$cell_type)
plotPCA(eset_all, groups=groups, groupnames=groupnames, main="GSE32591 WHO_class:No Class Type Specified, Affy GCRMA PCA",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17) )

eset.glom=eset_all[,eset_all$cell_type=="Glomeruli"]
eset.tub=eset_all[,eset_all$cell_type=="Tubulointerstitium"]

groups <- as.factor(eset.glom$cohort)
groupnames = levels(eset.glom$cohort)
plotPCA(eset.glom, groups=groups, groupnames=groupnames, main="GSE32591 Glom BA WHO_class:No Class Type Specified, Affy GCRMA PCA",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17))


groups <- as.factor(eset.tub$cohort)
groupnames = levels(eset.tub$cohort)
plotPCA(eset.tub, groups=groups, groupnames=groupnames, main="GSE32591 Tub BA WHO_class:No Class Type Specified, Affy GCRMA PCA",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17) )


# ###########################################################################
# ### OUTLIER REMOVAL
# # REMOVE OUTLIERS IDENTIFIED USING PCA
# table(eset_all$pca_outlier)
# table(eset_all$pca_outlier_basis)
# index <- which( eset_all$pca_outlier==TRUE)
# eset_all <- eset_all[ , -index ]
# # 
# # # REMOVE OUTLIER FROM NUSE PLOT OF RAW AFFY DATA (POOR RNA HYBRIDIZATION)
# # # THIS WAS NOT A B CELL ARRAY
# # index <- which(colnames(eset_all)=="GSM260892")
# # eset_all <- eset_all[ , -index ]




# #----------------- SUBSET ESET BASED WHO CLASS TYPE -------------------------------------------
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
# save(eset2_type2b,file="eset2b_BA.RData")
# save(eset2_noclasstype,file="eset2_noclasstype_BA.RData")
# 

index = which(eset.glom$WHO_class==3)
index = c(index,which(eset.glom$WHO_class==4))
index = c(index,which(eset.glom$cohort=="CTL"))
eset.glom.3.4 = eset.glom[,index]


groups <- as.factor(eset.glom.3.4$cohort)
groupnames = levels(eset.glom.3.4$cohort)
plotPCA(eset.glom.3.4, groups=groups, groupnames=groupnames, main="GSE32591_Glom_3-4 BA GCRMA PCA",
        pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels.default(sampleNames(eset.glom.3.4)) )

sampleNames(eset.glom.3.4)


######## ------------ ######## --------- ######## --------- ######## --------- ########
# Analysis of Tubulointerstitial class 2, WHO class type 2b
######## ------------ ######## --------- ######## --------- ######## --------- ########

# ---------- Plot Dendogram
distance <- dist(t(exprs(eset.glom.3.4)),method="euclidean")
clusters <- hclust(distance)
labels <- paste( eset.glom.3.4$friendlyName, eset.glom.3.4$WHO_class, eset.glom.3.4$WHO_class_type, sep="." )
labels <- gsub( ".control", "", labels)
plot(clusters, cex=1, main="GSE32591 WHO_class3.4, BA GCRMA,Euclidean distance" ,lwd=1, labels=labels)

# # ---------- plot PCA
# library(affycoretools)
# library(rgl)
# 
# eset2_type2b$cohort <- factor(eset2_type2b$cohort)
# labels <-  paste( eset2_type2b$friendlyName, eset2_type2b$WHO_class, eset2_type2b$WHO_class_type, sep="." )

# # Plot PCA of cohort in 2D
# groups <- as.numeric(eset2_type2b$cohort)
# groupnames = levels(eset2_type2b$cohort)
# plotPCA(eset2_type2b, groups=groups, groupnames=groupnames, main="GSE32591 WHO_class:2b BA GCRMA PCA",
#         pcs=c(1,2), plot3d=FALSE, col=c("red","blue"), pch=c(16,17), labels )

#---------- Plot intensity histogram

min(exprs(eset.glom.3.4))  #2.132413
max(exprs(eset.glom.3.4))  #15.88462

dev.off()
par(mfrow=c(1,2))

nrow <- format(as.numeric(nrow(eset.glom.3.4)),big.mark=",",scientific=F)
index <- which(rowMeans(exprs(eset.glom.3.4))<2.23)

hist(rowMeans(exprs(eset.glom.3.4)), breaks=100, xlab="Average log2 Probe Intensity",
     main="GSE32591 WHO_class3&4: No type specified, BA CDF All Intensities",
     sub=paste(nrow, " total probes"), xlim=c(2,12), col="green")

hist(rowMeans(exprs(eset.glom.3.4)), breaks=4000, xlab="Average log2 Probe Intensity",
     main="GSE32591 Intensities 2.13-2.24",
     sub=paste(nrow, " total probes. 6276 < mean 2.23", sep=""), xlim=c(2.13,2.24), col="green")
abline(v=2.23,col="red",lty=2)


nrow(data.frame(index)) 

BA.eset.glom.3.4.filtered = eset.glom.3.4[-index,]

# Examine if CTL vs. SLE distributions are different from each other
# (we see very similar distributions)
sapply(levels(BA.eset.glom.3.4.filtered$cohort), function(x) summary(as.vector(exprs(eset.glom.3.4.filtered)[,+  which(eset.glom.3.4.filtered$cohort == x)])))

###############################################################################################################
# MERGE WITH LIMMA/SAM DEG RESULTS
# NOTE! Use ALL LIMMA/SAM resuts, not just the significant DE probes

library(limma)

design <- model.matrix(~0 + BA.eset.glom.3.4.filtered$cohort)
colnames(design) <- levels(factor(BA.eset.glom.3.4.filtered$cohort))
fit <- lmFit(BA.eset.glom.3.4.filtered, design)
names(fit)
contrast.matrix <- makeContrasts(LN-CTL,levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top.limma <- topTable(fit2, coef=1, adjust="BH", n=nrow(BA.eset.glom.3.4.filtered))
top.limma$probe <- rownames(top.limma)
colnames(top.limma)

setwd(eset)
save(BA.eset.glom.3.4.filtered,file="BA.eset.glom.3.4.filtered.RData")
# # ------------------------- GENE FILTERING --------------------------------------------------------------
# 
# #------------ Filter genes without gene symbol
# # We filter out probes with empty gene names, a type of annotation based filtering
# index <- which(is.na(fData(eset2_type2b)[,"geneSymbol"])) # 9(BA gcrma) probes removed without gene mapping
# eset2b_noNAS <- eset2_type2b[-index,]
# eset2b_filtered <- eset2_type2b[-index,]
# 
# 
# #------------ Filter low intensity probes
# # Filter out very low intensity probes. We use a threshold of 2.82, which is the highest we can go before
# # we lose informative genes (as we found in the accompanying script find_lowpass_threshold.R)
# index <- which(rowMeans(exprs(eset2b_filtered)) < 2.82)  #4169
# eset2b_filtered <- eset2b_filtered[-index,]
# rm(index)
# 
# 
# #------------ ---LIMMA Analysis ---------------------------------------------------------------------
# # Run LIMMA on LN patients of WHO_CLASS: 2
# library(a4)
# 
# limmaResult <- limmaTwoLevels(eset2b_filtered, "cohort")
# propDEgenes(limmaResult) # 49% genes are differentially expressed
# 
# top.active <- topTable(limmaResult, n = nrow(eset2b_filtered), coef = 2)
# 
# # plot LIMMA results
# dev.off()
# volcanoPlot(limmaResult, main="GSE32591 WHO_class:2b BA CDF GCRMA Glomeruli LN and normal LIMMA",
#             xlab="Log2 Normalized Probe Intensity. x<0 is downregulated, x>0 is upregulated.",
#             ylab="Adjusted p-val of differential significance")
# 
# # RUN SAM ON CLASS2b LN PATIENTS AND MERGE WITH LIMMA RESULTS (NO LFC OR FDR FILTERING YET)
# library(siggenes)
# # assign numeric cohort designation
# SAM.cl = as.character(eset2b_filtered$cohort)
# SAM.cl[which(SAM.cl=="CTL")]=0
# SAM.cl[which(SAM.cl=="LN")]=1
# SAM.cl = as.numeric(SAM.cl)
# SAM.cl
# # run SAM
# SAM.out <- sam(exprs(eset2b_filtered), SAM.cl, rand = 1234567 )
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
# #top_merged$CDF <- "BA_HGU133A"
# top_merged$CDF <- "BrainArray_hgu133ahsentrezgcdf"
# 
# GSE32591_BA_GCRMA_Active <- top_merged
# 
# ## SAVE ACTIVE LUPUS VS. NORMALS MERGED LIMMA/SAM RESULTS
# setwd(eset)
# write.table(GSE32591_BA_GCRMA_Active ,file="GSE32591_Class2b_BA_Glom_LIMMA_SAM_ACTIVE_All.txt",sep="\t",quote=FALSE,row.names=FALSE)
# save(GSE32591_BA_GCRMA_Active, file="GSE32591_Class2b_BA_Glom_LIMMA_SAM_ACTIVE_All.RData")


#--------------------- HEATMAP : CLASS 2 - not class type specified-------------------------------------------------------------------
library(gplots)

# First subset out significant probes
index <- which(top.limma$adj.P.Val < 0.2) 

# top.limma.sig <- top.limma[ index, ]
top.limma.sig <- top.limma[ c(1:1000), ]

# Use the rownames of top.limma.sig as a lookup column for the eset
top <- rownames(top.limma.sig)
index <- which(rownames(BA.eset.glom.3.4.filtered)%in%top)
top <- exprs(BA.eset.glom.3.4.filtered)[ index, ]
m.use = as.matrix(top)
rownames(m.use) <- fData(BA.eset.glom.3.4.filtered)[index,"geneSymbol"]
colnames(m.use) = BA.eset.glom.3.4.filtered$friendlyName
colnames(m.use) = paste( BA.eset.glom.3.4.filtered$friendlyName, BA.eset.glom.3.4.filtered$WHO_class, BA.eset.glom.3.4.filtered$WHO_class_type, sep="." )
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
                       main="GSE32591 Class 3&4, BA CDF GCRMA Glomeruli LN-CTL LIMMA, top 1k")

# # EXAMINE CLUSTERS OF PROBES
# hc_probes = hclust(dist(m.use),method = "complete")
# hc_probes.cut = cutree(hc_probes,k=8) # k = number of clusters
# table(hc_probes.cut) # show cluster cut results
# 
# # Now extract probes within clusters
# clusters <- data.frame(cluster=hc_probes.cut,probe=names(hc_probes.cut))
# clusters <- clusters[ order(clusters$cluster, clusters$probe), ]
# 
# # Save the clustered probes
# GSE32591_Class2_notype_GLOM_GCRMA_Active_clusters <- clusters
# setwd(eset)
# save(GSE32591_Class2_notype_GLOM_GCRMA_Active_clusters,file="GSE32591_Class2_notype_Affy_GLOM_GCRMA_Active_clusters.RData")
