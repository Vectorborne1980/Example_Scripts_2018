# This script is identifying triplets of miRNA, TF, and target gene instances instance by the type of relationship between
# miRNA and TF (type 1 if TF activates miRNA; Type 2 if miRNA inhibits TF)
# Inputs:
#     input data file folder
#     names of .xlsx files (without the .xlsx extension)
#     
# Outputs: 
#     .csv file for each input data file


# Install packages and set path to Java library

#install.packages("xlsx")
#install.packages("openxlsx")
#Be sure that the version of Java (32 or 64-bit) is the same as the type of R used or else rJava will not load

library(rJava)
library(openxlsx)
library(xlsx)

# set path to the directory containing the input .xlsx data files
inputFileDIR<-"C:\\Users\\gera501\\Desktop\\Magia\\"
# list the names of the input data files (without the .xlsx extension)
inputDataFiles<-c("AH1_7hr_Magia_Output","AH1_12hr_Magia_Output","FM_7hr_Magia_Output","FM_12hr_Magia_Output")  # data files names without .xlsx extension, and be sure to format the input files according to the styles seen in these examples

# set worksheet names (assuming they are the same in all data files)
wksht_miRNA_target<-c("miR-Target")
wksht_tf_target<-c("TF-Target")
wksht_edges<-c("Single Network")


num_files<-length(inputDataFiles)  # number of files to read

for (i in 1:num_files){
  xlsxfile<-paste(inputFileDIR,inputDataFiles[i],".xlsx",sep="")
  miR_targets<-openxlsx::read.xlsx(xlsxfile,sheet=wksht_miRNA_target)
  tf_targets<-openxlsx::read.xlsx(xlsxfile,sheet=wksht_tf_target)
  edge_relations<-openxlsx::read.xlsx(xlsxfile,sheet=wksht_edges)
  
  
  # merge by target common to TF and miRNA
  
  miRNA_target_tf<-merge(miR_targets,tf_targets,by=c("Target"))
  
  
  unique_tf<-unique(miRNA_target_tf$TF)
  unique_miRNA<-unique(miRNA_target_tf$microRNA)
  unique_target<-unique(miRNA_target_tf$Target)
  
  # Triplets of type 1
  #TF targets(activates) miRNA; miRNA targets gene; gene targeted_by TF
  
  
  i1<-edge_relations$Origin %in% unique_tf & edge_relations$Edge %in% "Activation" & 
    edge_relations$Target %in% unique_miRNA
  edges<-edge_relations[i1,4:6]
  names(edges)<-c("TF","Edge","microRNA")
  # type 1 set
  miRNA_targetgenes_tf_type1<-merge(miRNA_target_tf,edges,by=c("microRNA","TF"))
  names(miRNA_targetgenes_tf_type1)[4]<-c("Type")
  miRNA_targetgenes_tf_type1$Type<-1
  
  write.csv(miRNA_targetgenes_tf_type1,paste("type1_output_",inputDataFiles[i],".csv",sep=""))
  
  # Triplets of type 2
  # miRNA inhibits (represses) TF
  i1<-edge_relations$Origin %in% unique_miRNA & edge_relations$Edge %in% "Repression" & 
    edge_relations$Target %in% unique_tf
  
  edges<-edge_relations[i1,4:6]
  names(edges)<-c("microRNA","Edge","TF")
  
  miRNA_targetgenes_tf_type2<-merge(miRNA_target_tf,edges,by=c("microRNA","TF"))
  names(miRNA_targetgenes_tf_type2)[4]<-c("Type")
  miRNA_targetgenes_tf_type2$Type<-2
  
  write.csv(miRNA_targetgenes_tf_type2,paste("type2_output_",inputDataFiles[i],".csv",sep=""))

}
