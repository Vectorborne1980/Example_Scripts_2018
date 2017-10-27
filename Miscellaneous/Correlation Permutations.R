setwd("E:\\PNNL\\OMICS-LHV (Lethal Human Viruses)\\Infection Experiments\\ICL102\\miRNomics\\Probe-Based Correlations")
dat.raw <- read.delim("ICL102_Paired_miR-mRNA_Log2_Ratios_PROBES.txt")
d.mir <- read.delim(file = "clipboard", as.is =T)
d.mrna <- read.delim(file = "clipboard", as.is = T)
pairs <- read.delim(file = "clipboard", as.is=T)

# Correlations =======================================================================================

# Function for correlation coefficient of 1 pair -----------------------------------------------------
corGetter <- function(probeID, mirID, mytime, myvirus)
{
  mymir <- d.mir[which(d.mir$miRBase_Accession == mirID)[1],]
  mymir<- mymir[,-c(1:2)]
  if(length(mytime)>0){mymir <- mymir[,grep(mytime, names(mymir))]}
  if(length(myvirus)>0){mymir <- mymir[,grep(myvirus, names(mymir))]}
  
  mymrna <- d.mrna[which(d.mrna$ProbeID == probeID)[1],]
  
  x<-as.vector(as.matrix(mymir))
  y <- as.vector(as.matrix(mymrna[, names(mymrna) %in% names(mymir)]))
  
  return(cor(x,y, use = "pairwise.comp"))
}

# Correlation Permutations ----------------------------------------------------------------------------
N <- 100000
permutCor <- rep(NA, N)

for(i in 1:N)
{
  thisProbe <- as.character(sample(d.mrna$ProbeID, 1, replace = F))
  thisMir <- as.character(sample(d.mir$miRBase_Accession, 1, replace = F))
  
  permutCor[i] <- corGetter(probeID = thisProbe, mirID = thisMir, mytime = "12h", myvirus = "AH1")
  
}
hist(permutCor, breaks = 100)


# Calculate the actual p-value for a pair of interest for one pair ------------------------------------
realCor <- corGetter(probeID = "A_33_P3888365", mirID = "MIMAT0000680", mytime = "12h", myvirus = "AH1")
realPval <- sum((na.omit(permutCor)) < (na.omit(realCor)))/N


# Now for all pairs -----------------------------------------------------------------------------------
# pairs <- pairs.BU[1:5000,]
realPval <- rep(NA, length(pairs$ProbeID))
realCor <- rep(NA, length(pairs$ProbeID))

thisHr <- "12h"
thisVirus = "AH1"

for (i in 1:length(realPval))
{
  thisProbe <- pairs$ProbeID[i]
  thisMir <- pairs$miRBase_Accession[i]
  
  realCor[i] <- corGetter(probeID = thisProbe, mirID = thisMir, mytime = thisHr, myvirus = thisVirus )
  realPval[i] <- sum((na.omit(permutCor))<= (na.omit(realCor[i])))/N
}

capture.output(write.table(realPval, file = "pval_AH1_12h.txt"))
capture.output(write.table(realCor, file = "cor_AH1_12h.txt"))

p.corrected <- p.adjust(realPval, method = "fdr")

######################



# Linear Regressions ================================================================================

# Function for linear regression p-value of 1 pair --------------------------------------------------
pGetter <- function(probeID, mirID, mytime, myvirus)
{
  mymir <- d.mir[which(d.mir$miRBase_Accession == mirID)[1],]
  mymir<- mymir[,-c(1:2)]
  if(length(mytime)>0){mymir <- mymir[,grep(mytime, names(mymir))]}
  if(length(myvirus)>0){mymir <- mymir[,grep(myvirus, names(mymir))]}
  
  mymrna <- d.mrna[which(d.mrna$ProbeID == probeID)[1],]
  
  x<-as.vector(as.matrix(mymir))
  y <- as.vector(as.matrix(mymrna[, names(mymrna) %in% names(mymir)]))
  
  #myBeta <- as.numeric(lm(y~x)$coefficients[2])
  myp <- summary(lm(y~x))$coefficients[2,4]
  
  return(myp)
}

# Regression Permutations ----------------------------------------------------------------------------
N <- 100000
pVec <- rep(NA, N)

for(i in 1:N)
{
  thisProbe <- as.character(sample(d.mrna$ProbeID, 1, replace = T))
  thisMir <- as.character(sample(d.mir$miRBase_Accession, 1, replace = T))
  
  pVec[i] <- pGetter(probeID = thisProbe, mirID = thisMir, mytime = "12hr", myvirus = "AH1")
  
}
hist(pVec, breaks = 100)


# Calculate the actual p-value for a pair of interest for one pair -----------------------------------
realp <- pGetter(probeID = "A_33_P3217393", mirID = "MIMAT0003339", mytime = "12hr", myvirus = "AH1")
realPval <- sum((na.omit(pVec))<=realp)/N

# Now for all pairs ----------------------------------------------------------------------------------
pairs <- pairs.BU[1:5000,]
realVec <- rep(NA, length(pairs$ProbeID))

thisHr <- "12hr"
thisVirus = "AH1"

for (i in 1:length(realVec))
{ 
  thisProbe <- pairs$ProbeID[i]
  thisMir <- pairs$miRBase_Accession[i]
  
  realp[i] <- pGetter(probeID = thisProbe, mirID = thisMir, mytime = thisHr, myvirus = thisVirus )
  realVec[i] <- sum((na.omit(pVec))<=(na.omit(realp[i])))/N
}

p.corrected <- p.adjust(realVec, method = "fdr")



