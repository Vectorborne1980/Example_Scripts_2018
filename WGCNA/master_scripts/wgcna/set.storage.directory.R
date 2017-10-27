
# Set directory path to save these plots in, formatted as 
# /mnt/disks/rserver/projects/projectName/wgcna/DSpower.pamBool.minMod-Num/plots/blocks/
# The function will traverse into the final directory, creating them along the way as needed

setStorageDir <- function (storageDir) {
  # Create directories if they don't exist. This is a primitive script that works by moving into or creating 
  # directories using dir.create. Note I'm not able to pass a directory path to this function, such as
  # dir.create("plots/blocks/DS10.pamTRUE.minMod-100"), so have to generate these one at time.
  newDir <- vector('character')
  for (i in 1:length(storageDir)) {
    newDir <- storageDir[i]
    if (!dir.exists(newDir)) {
      dir.create(newDir, showWarnings=TRUE, mode="0777")
      print(paste("Sub-directory",newDir,"created."))
    } else {
      print(paste("Sub-directory",paste(getwd(),"/",newDir," exists.",sep="") ))
    }
    setwd( paste(getwd(),newDir,sep="/") )
  }
  print(paste("Current storage directory: ",getwd(),sep=""))
  rm(storageDir,newDir,i)
}
