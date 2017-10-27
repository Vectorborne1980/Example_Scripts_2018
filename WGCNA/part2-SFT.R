################################################################################################
### WGCNA PART 2

# This script uses the datExpr and datTraits to establish a soft threshold

################################################################################################
### LOAD AND CONFIGURE WGCNA LIBRARY
analysisDir <- "/mnt/disks/rserver-hdd/projects/Nick/B-cells/CD19_B-cells/Standardized_CTL_SLEDAI_NA/"
library(WGCNA)
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

################################################################################################

# Load files
setwd(analysisDir)
load("GSE10325_CD19_Standardized_CTL_SLEDAI_NA_analysis_parameters.RData")
load("part1-final-objects.RData")
# Choose a set of soft-thresholding powers
# First generate a sequence of powers to apply
powers <- c(c(1:10), seq(from = 12, to=30, by=2))

# Call the network topology analysis function
# 07-19-17 RDR: I had to set blockSize to 5,000 when using a Google VM with 4 VCPUs and 15 GB memory
# BlockSize of 10,000 quickly ran on 8 VCPUs and 52 Gb memory
sft <- pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers, verbose = 5, networkType="signed", corFnc = "bicor",
                         moreNetworkConcepts = TRUE, blockSize=10000) #bicor is the only option

pdf("SFT.pdf",width = 16,height = 9)
par(mfrow = c(1,2));
cex1 <- 0.8; # font size
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab=paste(active.base,"STP Selection",sep=" "),ylab=paste(active.base,"Scale Free Topology Model Fit, signed",sep=" "),type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
###SELECT THE THRESHOLD NUMBER AT WHICH THE SFT MODEL FIT PLATEAUS
###NOTE: If graphs are unusual, consider re-filtering, OR if sample numbers are low, use the SFT values on the WGCNA FAQ page table

save(sft, file="part2-output-SFT.RData")

sft.fits <- sft$fitIndices
write.table(sft.fits,file="part2-output-STP.txt", quote=FALSE, sep = "\t",row.names=FALSE, col.names=TRUE)
