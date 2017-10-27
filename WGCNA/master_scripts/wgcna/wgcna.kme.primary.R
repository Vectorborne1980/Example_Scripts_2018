
###############################################################################################
# DS1 KMEs
KME.corr.DS1 <- KME.MEs.DS1[,seq(1, ncol(KME.MEs.DS1), 2)]
colnames(KME.corr.DS1) <- substr( colnames(KME.corr.DS1), 5, nchar(colnames(KME.corr.DS1)))
KME.corr.p.DS1 <- KME.MEs.DS1[,seq(2, ncol(KME.MEs.DS1), 2)]
colnames(KME.corr.p.DS1) <- colnames(KME.corr.DS1)

KME.primary.DS1 <- matrix(nrow=nrow(KME.MEs.DS1),ncol=3)

# RDR: This for loop is terribly inefficient, please recode.
print(paste("Retrieving primary KME values for DS1..."))
for (i in 1:nrow(geneList.DS1)) {
  for (j in 1:ncol(KME.corr.DS1)) {
    if (geneList.DS1[i,"moduleColors.DS1"]==colnames(KME.corr.DS1)[j] ) { 
      KME.primary.DS1[i,1] <- geneList.DS1[i,"moduleColors.DS1"]
      KME.primary.DS1[i,2] <- KME.corr.DS1[i,j]
      KME.primary.DS1[i,3] <- KME.corr.p.DS1[i,j]
    }
  }
}
colnames(KME.primary.DS1) <- c("KME.primary","KME.color","KME.p.color")
rm(KME.corr.DS1,KME.corr.p.DS1)

###############################################################################################
# DS2 KMEs
KME.corr.DS2 <- KME.MEs.DS2[,seq(1, ncol(KME.MEs.DS2), 2)]
colnames(KME.corr.DS2) <- substr( colnames(KME.corr.DS2), 5, nchar(colnames(KME.corr.DS2)))
KME.corr.p.DS2 <- KME.MEs.DS2[,seq(2, ncol(KME.MEs.DS2), 2)]
colnames(KME.corr.p.DS2) <- colnames(KME.corr.DS2)

KME.primary.DS2 <- matrix(nrow=nrow(KME.MEs.DS2),ncol=3)

# RDR: This for loop is terribly inefficient, please recode
print(paste("Retrieving primary KME values for DS2..."))
for (i in 1:nrow(geneList.DS2)) {
  for (j in 1:ncol(KME.corr.DS2)) {
    if (geneList.DS2[i,"moduleColors.DS2"]==colnames(KME.corr.DS2)[j] ) { 
      KME.primary.DS2[i,1] <- geneList.DS2[i,"moduleColors.DS2"]
      KME.primary.DS2[i,2] <- KME.corr.DS2[i,j]
      KME.primary.DS2[i,3] <- KME.corr.p.DS2[i,j]
    }
  }
}
colnames(KME.primary.DS2) <- c("KME.primary","KME.color","KME.p.color")
rm(KME.corr.DS2,KME.corr.p.DS2)

###############################################################################################
# DS3 KMEs
KME.corr.DS3 <- KME.MEs.DS3[,seq(1, ncol(KME.MEs.DS3), 2)]
colnames(KME.corr.DS3) <- substr( colnames(KME.corr.DS3), 5, nchar(colnames(KME.corr.DS3)))
KME.corr.p.DS3 <- KME.MEs.DS3[,seq(2, ncol(KME.MEs.DS3), 2)]
colnames(KME.corr.p.DS3) <- colnames(KME.corr.DS3)

KME.primary.DS3 <- matrix(nrow=nrow(KME.MEs.DS3),ncol=3)

# RDR: This for loop is terribly inefficient, please recode
print(paste("Retrieving primary KME values for DS3..."))
for (i in 1:nrow(geneList.DS3)) {
  for (j in 1:ncol(KME.corr.DS3)) {
    if (geneList.DS3[i,"moduleColors.DS3"]==colnames(KME.corr.DS3)[j] ) { 
      KME.primary.DS3[i,1] <- geneList.DS3[i,"moduleColors.DS3"]
      KME.primary.DS3[i,2] <- KME.corr.DS3[i,j]
      KME.primary.DS3[i,3] <- KME.corr.p.DS3[i,j]
    }
  }
}
colnames(KME.primary.DS3) <- c("KME.primary","KME.color","KME.p.color")
rm(KME.corr.DS3,KME.corr.p.DS3)

###############################################################################################
# DS4 KMEs
KME.corr.DS4 <- KME.MEs.DS4[,seq(1, ncol(KME.MEs.DS4), 2)]
colnames(KME.corr.DS4) <- substr( colnames(KME.corr.DS4), 5, nchar(colnames(KME.corr.DS4)))
KME.corr.p.DS4 <- KME.MEs.DS4[,seq(2, ncol(KME.MEs.DS4), 2)]
colnames(KME.corr.p.DS4) <- colnames(KME.corr.DS4)

KME.primary.DS4 <- matrix(nrow=nrow(KME.MEs.DS4),ncol=3)

# RDR: This for loop is terribly inefficient, please recode
print(paste("Retrieving primary KME values for DS4..."))
for (i in 1:nrow(geneList.DS4)) {
  for (j in 1:ncol(KME.corr.DS4)) {
    if (geneList.DS4[i,"moduleColors.DS4"]==colnames(KME.corr.DS4)[j] ) { 
      KME.primary.DS4[i,1] <- geneList.DS4[i,"moduleColors.DS4"]
      KME.primary.DS4[i,2] <- KME.corr.DS4[i,j]
      KME.primary.DS4[i,3] <- KME.corr.p.DS4[i,j]
    }
  }
}
colnames(KME.primary.DS4) <- c("KME.primary","KME.color","KME.p.color")
rm(KME.corr.DS4,KME.corr.p.DS4)

###############################################################################################

save(KME.primary.DS1,KME.primary.DS2,KME.primary.DS3,KME.primary.DS4,file="gene_KMEs.RData")

