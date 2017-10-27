#####GSE36700
setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/GAGE/Synovium/vs_BIGC_4-3/SLE-vs-RA/")
load("C:/Users/Nick Geraci/Nick/AMPEL/Software/WGCNA/WGCNA_Standardized_Re-analyses/Synovium/SLE_vs_RA/NSG_scripts/GSE36700_Synovium_Standardized/GSE36700-Synovium_Active_WGCNA_Traits-Expression_NSG_2017-08-21.RData")
exprs.data.synovium <- t(datExpr)
colnames(exprs.data.synovium) <- sub("RA","CTL",colnames(exprs.data.synovium))
write.table(exprs.data.synovium, file = "synovium_wgcna_data.txt",sep = "\t")

####Consensus modules
# syn.mods<-read.delim(file="synovium_modules_probes.txt",header=T,sep="\t")
# syn.mods[synovium.mods==""] <- NA # Convert the empty cells to "NA"
# syn.mods <- as.list(syn.mods[,1:16])
# syn.mods <- lapply(syn.mods,function(x) x[!is.na(x)])
# # Convert Probe List IDs
# genemap<-read.delim("./Conversion_map.txt",header=T,sep="\t",stringsAsFactors=F,row.names = 1)
# # Convert the DE lists to the best matching ID, symbol in this case
# for(i in 1:length(synovium.mods)) {
# syn.mods[[i]] <- genemap[match(syn.mods[[i]],rownames(genemap)),1]
# syn.mods[[i]] <- unique(syn.mods[[i]])
# syn.mods[[i]] <- syn.mods[[i]][which(syn.mods[[i]]!="")]
# syn.mods[[i]] <- syn.mods[[i]][which(!is.na(syn.mods[[i]]))]
# }
# rm(i,genemap)

syn.mods <- data.frame(read.delim(file = "C:/Users/Nick Geraci/Nick/AMPEL/Software/WGCNA/WGCNA_Standardized_Re-analyses/Synovium/SLE_vs_RA/NSG_scripts/GSE36700_Synovium_Standardized/GSE36700-Synovium_Active_WGCNA_DS4-Gene-Info-Minus-Dups_NSG_2017-08-21.txt", header=T,sep="\t",stringsAsFactors=F))
syn.mods <- syn.mods[,c(1,2,4,6)]

uniq<-unique(syn.mods[,4]) # get unique syn.mods categories
syn.mods.lists<-vector("list",length=length(uniq)) # set up data structure
names(syn.mods.lists)<-uniq   # name the structure with the unique names
for(i in 1:length(syn.mods.lists)) # loop through getting genes associated with each entry
  syn.mods.lists[[i]]<-syn.mods[which(syn.mods[,4]==uniq[i]),3]
max.length <- max(sapply(syn.mods.lists, length)) # Get the maximum number of genes mapped to any one BIGC term
syn.mods.lists <- lapply(syn.mods.lists, function(v) { c(v, rep(NA, max.length-length(v)))}) #Add NA values to list elements to make them all the same length
rm(i,uniq,max.length)
syn.mods.df <- do.call(rbind,syn.mods.lists)
syn.mods.df <- t(syn.mods.df)
syn.mods.df <- as.data.frame(syn.mods.df)
syn.mods.df <- syn.mods.df[,order(colnames(syn.mods.df))]
write.table(syn.mods.df,"synovium_modules_entrez.txt",sep = "\t")

uniq<-unique(syn.mods[,4]) # get unique syn.mods categories
syn.mods.lists<-vector("list",length=length(uniq)) # set up data structure
names(syn.mods.lists)<-uniq   # name the structure with the unique names
syn.mods$probe <- as.character(syn.mods$probe)
for(i in 1:length(syn.mods.lists)) # loop through getting genes associated with each entry
  syn.mods.lists[[i]]<-syn.mods[which(syn.mods[,4]==uniq[i]),1]
max.length <- max(sapply(syn.mods.lists, length)) # Get the maximum number of genes mapped to any one BIGC term
syn.mods.lists <- lapply(syn.mods.lists, function(v) { c(v, rep(NA, max.length-length(v)))}) #Add NA values to list elements to make them all the same length
rm(i,uniq,max.length)
syn.mods.df <- do.call(rbind,syn.mods.lists)
syn.mods.df <- t(syn.mods.df)
syn.mods.df <- as.data.frame(syn.mods.df)
syn.mods.df <- syn.mods.df[,order(colnames(syn.mods.df))]
write.table(syn.mods.df,"synovium_modules_probes.txt",sep = "\t")

####BigC.v3 Reference
BigC <- data.frame(read.delim(file = "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/BIG-C/BIGCv4-3_NSG_17Aug2017.txt",
                              header=T,sep="\t",stringsAsFactors=F))
BigC <- BigC[order(BigC$BIG.C_Category),]
BigC <- BigC[which(BigC$Affymetrix_HGU133.Plus2==1),] #MODIFY to define the chip type for the dataset
uniq<-unique(BigC[,3]) # get unique BigC categories
BigC.lists<-vector("list",length=length(uniq)) # set up data structure
names(BigC.lists)<-uniq   # name the structure with the unique names
for(i in 1:length(BigC.lists)) # loop through getting genes associated with each entry
  BigC.lists[[i]]<-BigC[which(BigC[,3]==uniq[i]),2]
# BigC.lists[BigC.lists==""] <- NA #Replace blank spots with NAs
# BigC.lists <- lapply(BigC.lists,function(x) x[!is.na(x)]) #If the entry is NA, remove it
rm(i,uniq)

####GAGE
library(gage)
library(rJava)
library(openxlsx)
library(xlsx)
library(gtools)
#Modify to the specific platform for the dataset
genemap<-read.delim("C:/Users/Nick Geraci/Nick/AMPEL/Datasets/CDFs/hgu133_plus2_Affy_CDF.txt",
                    header=T,sep="\t",stringsAsFactors=F,row.names = 1)
genemap<-genemap[,c(3,1)] #Modify to the symbol and entrez columns

setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/GAGE/Synovium/vs_BIGC_4-3/SLE-vs-RA/module_results/")

for(a in 1:ncol(syn.mods.df)){
  x <- exprs.data.synovium[is.element(rownames(exprs.data.synovium),syn.mods.df[,a]),]
  x2 <- as.list(rownames(x))
  for(b in 1:length(x2)) {
    x2[[b]] <- genemap[match(x2[[b]],rownames(genemap)),1]
    x2[[b]] <- unique(x2[[b]])
    x2[[b]] <- x2[[b]][which(x2[[b]]!="")]
    x2[[b]] <- x2[[b]][which(!is.na(x2[[b]]))]
    }
  rownames(x) <- x2
  rm(b,x2)
  nrow(x)
  x <- as.data.frame(x)
  x <- as.matrix(x)
  cn=colnames(x)
  SLE=grep('SLE',cn, ignore.case =T) 
  CTL=grep('CTL',cn, ignore.case =T)
  gage.test <- gage(x,gsets = BigC.lists, ref = CTL, samp = SLE, compare = 'unpaired', same.dir = F)
  stats <- as.data.frame(gage.test$stats)
  stats <- stats[order(rownames(stats)),]
  greater <- as.data.frame(gage.test$greater)
  greater <- greater[order(rownames(greater)),]
  df <- data.frame(matrix(nrow = length(names(BigC.lists)),ncol = 8))
  colnames(df) <- c("Categories","Count","Overlapping_Module_Genes","Percent_Category_Represented_in_Module",
                    "Pval","(neg)log10_Pval","Qval","(neg)log10_Qval")
  df[,1] <- c(names(BigC.lists))
  for(c in 1:length(names(BigC.lists))){
    df[c,2]<-length(BigC.lists[[c]])
    }
  df[,3] <- greater[,5]
  df[,5] <- greater[,3]
  df[,7] <- greater[,4]
  df[,4] <- (df[,3]/df[,2])*100
  df[,6] <- -log(df[,5],base = 10)
  df[,8] <- -log(df[,7],base = 10)
  df[df==""] <- NA
  wb <- createWorkbook()
  ws <- createSheet(wb=wb, sheetName="GAGE_greater")
  addDataFrame(x=greater, sheet=ws, row.names=TRUE)
  ws <- createSheet(wb=wb, sheetName = "GAGE_stats")
  addDataFrame(x=stats, sheet=ws, row.names=TRUE)
  ws <- createSheet(wb=wb, sheetName = "BIGC_stats")
  addDataFrame(x=df, sheet=ws, row.names=FALSE)
  saveWorkbook(wb,file=paste(colnames(syn.mods.df[a]),"_gage_both-directions.xlsx",sep = ""))
  }

inputFileDIR<-"C:/Users/Nick Geraci/Nick/AMPEL/Software/GAGE/Synovium/vs_BIGC_4-3/SLE-vs-RA/module_results/"
setwd(inputFileDIR)
file.list <- list.files(inputFileDIR)
inputDataFiles<-c(file.list)
wksht<-c("BIGC_stats")
num_files<-length(inputDataFiles)  # number of files to read
DF <- data.frame(matrix(nrow = length(names(BigC.lists))))
for (i in 1:num_files){
  xlsxfile <- paste(inputFileDIR,inputDataFiles[i],sep="")
  bigc <- openxlsx::read.xlsx(xlsxfile,sheet=wksht)
  DF <- cbind(DF,bigc[8])
}
DF[,1] <- c(names(BigC.lists))
colnames(DF) <- c("Modules",colnames(syn.mods.df))
wb <- createWorkbook()
ws <- createSheet(wb=wb, sheetName = "Summary")
addDataFrame(x=DF, sheet=ws, row.names=TRUE)
setwd("C:/Users/Nick Geraci/Nick/AMPEL/Software/GAGE/Synovium/vs_BIGC_4-3/SLE-vs-RA/")
saveWorkbook(wb,file="Synovium_GAGE_Results_Summary_NSG_21Aug2017.xlsx")
