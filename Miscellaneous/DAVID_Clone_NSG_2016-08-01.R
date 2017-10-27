## Loading lists
# Load the "list of lists" files copied from Excel
Kidney.WGCNA <- data.frame(read.delim(file = "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/DAVID Clone Analysis/data.txt",as.is=T))
# Convert the empty cells to "NA"
Kidney.WGCNA[Kidney.WGCNA==""] <- NA
# Change the data frame to a list with multiple elements, with "1:22" meaning columns 1 to 22, or however many columns there are
Kidney.WGCNA.List <- as.list(Kidney.WGCNA[,1:14])
# Now, remove all of the NA without removing all other values in the rows with NA
Kidney.WGCNA.List <- lapply(Kidney.WGCNA.List,function(x) x[!is.na(x)])

# Convert Probe List IDs
genemap<-read.delim("C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/DAVID Clone Analysis/genemap.txt",header=T,sep="\t",stringsAsFactors=F,row.names = 1)

# Get gene universe
universe <- data.frame(read.delim(file="C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/DAVID Clone Analysis/universe.txt",header = F))
universe <- as.vector(universe$V1)
 
# Convert the DE lists to the best matching ID, symbol in this case
for(i in 1:length(Kidney.WGCNA.List)) {
  Kidney.WGCNA.List[[i]] <- genemap[match(Kidney.WGCNA.List[[i]],rownames(genemap)),1]
  Kidney.WGCNA.List[[i]] <- unique(Kidney.WGCNA.List[[i]])
  Kidney.WGCNA.List[[i]] <- Kidney.WGCNA.List[[i]][which(Kidney.WGCNA.List[[i]]!="")]
  Kidney.WGCNA.List[[i]] <- Kidney.WGCNA.List[[i]][which(!is.na(Kidney.WGCNA.List[[i]]))]
}

# Make the ID conversion in the universe list
universe <- genemap[match(universe,rownames(genemap)),1]
# Take only the unique items
universe <- unique(universe)
# Get rid of match failures
universe <- universe[which(!is.na(universe))]
universe <- universe[which(universe!="")]

# BigC Reference
BigC <- data.frame(read.delim(file = "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/Gene_lists/BigC-v3-0_Unique-Genes-by-Category_NSG_2016-07-28.txt", header=T,sep="\t",stringsAsFactors=F,row.names = 1))
uniq<-unique(BigC[,2]) # get unique BigC categories
BigC.lists<-vector("list",length=length(uniq)) # set up data structure
names(BigC.lists)<-uniq   # name the structure with the unique names
for(i in 1:length(BigC.lists)) # loop through getting genes associated with each kegg entry
  BigC.lists[[i]]<-BigC[which(BigC[,2]==uniq[i]),1]
BigC.lists[BigC.lists==""] <- NA
BigC.lists <- lapply(BigC.lists,function(x) x[!is.na(x)])

# Load the 3 column gene map
map <- data.frame(read.delim(file = "C:/Users/Nick Geraci/Nick/AMPEL/Datasets/GSE32591_Berthier_Kidney/Nick_Revised_GSE32591/DAVID Clone Analysis/map.txt", header=T,sep="\t",stringsAsFactors=F,row.names = 1))

## Run the enrichment
source("C:/Users/Nick Geraci/Nick/PNNL/Software/DAVID Functional Enrichment Clone/newkappa_v4.R")
dave1 <- davidSimple(Kidney.WGCNA.List, BigC.lists, universe)
dave2 <- finish(dave1, geneInfoMap = map)
# capture.output(write.table(dave2$finalMat,file="CD4_Kidney_Davis-Grun_WGCNA_Module-Gene-Lists-UNIQUE-MAP_David-Clone-Results_NSG_2016-07-28.txt"))













