library(rJava)
library(openxlsx)
library(xlsx)
library(gtools)

inputFileDIR <- "C:/Users/Nick Geraci/Desktop/MCODE Eigengene Correlations/"
setwd(inputFileDIR)
# list the names of the input data files (without the .xlsx extension)
file.list <- list.files(inputFileDIR)
file.list <- file.list[substr(file.list,nchar(file.list)-4,nchar(file.list))==".xlsx"]
# data files names without .xlsx extension, and be sure to format the input files according to the styles seen in these examples
file.list <- gsub(".xlsx","",file.list)
inputDataFiles<-c(file.list)
wksht_mcode<-c("MCODE")
num_files<-length(inputDataFiles)  # number of files to read
DF <- data.frame()
for (i in 1:num_files){
  xlsxfile<-paste(inputFileDIR,inputDataFiles[i],".xlsx",sep="")
  mcode<-openxlsx::read.xlsx(xlsxfile,sheet=wksht_mcode)
  df <- cbind(mcode$moduleColor,mcode$MCODE_cluster,mcode$geneSymbol)
  DF <- rbind(DF, df)
}
colnames(DF) <- c("Color","Cluster","Symbol")
DF <- DF[which(!is.na(DF$Cluster)),]

# symbols <- mcode$geneSymbol
# library(dplyr)
# 
# ID <- c(1,1,1,2,2,2,2,3,3)
# Value <- c(2,3,5,2,5,8,17,3,5)
# Event <- c(1,1,2,1,2,1,2,2,2)
# group <- data.frame(Subject=ID, pt=Value, Event=Event)
# group %>% group_by(Subject) %>% top_n(1, pt)
# 
# for (i in 1:length(color_index) ) {
#   single.mods <- geneInfo[ geneInfo$moduleColor==color_index[i], ]
#   single.mods <- single.mods[ !is.na(single.mods$moduleColor), ]
#   single.mods <- single.mods[ ,c(1:14) ]
#   ws <- createSheet(wb=wb, sheetName=paste(toupper(color_index[i])))
#   addDataFrame(x=single.mods, sheet=ws, row.names=FALSE)
# }


library(dplyr)
# require(XLConnect)
# wb1 <- loadWorkbook("C:/Users/Nick Geraci/Desktop/Skin/Skin_Standardized/GSE72535-Skin_WGCNA_Final-Standardized-Analysis-Results_NSG_2017-03-06.xlsx")
# mcode <- readWorksheet(wb1, sheet = getSheets(wb1))
# mcode <- mcode[-c(1,2)]
geneInfo.list <- split(geneInfo, f = geneInfo$moduleColor)
data.list <- lapply(geneInfo.list,function(x){
  as.data.frame(x) %>%
    group_by(Gene.ID) %>%
    top_n(1,kWithin)})
library (plyr)
df <- ldply (data.list, data.frame)
df <- df[,-1]
rownames(df) <- df$probe

options(java.parameters = "-Xmx16000m") ####Boost the the maximum RAM available on your machine (i.e., 16GB is -Xmx16000m)
dir <- "C:/Users/Nick Geraci/Desktop/Skin/Skin_Standardized/"
setwd(dir)

# write_list <-function(data, wb_name=wb) {    
#   writeWorksheet(wb, data, names(data),header=FALSE)
# }
# for (i in 1:length(color_index) ) {
#   nrow(df)
#   ws <- createSheet(wb=wb, sheetName=paste(toupper(color_index[i])))
#   addDataFrame(x=df, sheet=ws, row.names=FALSE)
# }

for (i in 1:length(color_index) ) {
  single.mods <- df[ df$moduleColor==color_index[i], ]
  single.mods <- single.mods[ !is.na(single.mods$moduleColor), ]
  single.mods <- single.mods[ order(single.mods$p.GS.CLASI.A), ] ######################MODIFY the primary clinical trait
  single.mods <- single.mods[ ,c(1:14) ] ######################MODIFY for the number of columns from df
  nrow(single.mods)
  ws <- createSheet(wb=wb, sheetName=paste(toupper(color_index[i])))
  addDataFrame(x=single.mods, sheet=ws, row.names=FALSE)
}



