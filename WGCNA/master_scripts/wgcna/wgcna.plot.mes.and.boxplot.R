
ME.boxplot<-function(MEs,index) {
  quantile.ME<-quantile(MEs[,index])
  median.ME<-round(quantile.ME[3],2)
  max.ME<-round(quantile.ME[5],2)
  diff.ME<-max.ME-median.ME
  par(mar=c(0,3,4,3))
  boxplot(MEs[,index],main=paste0(colnames(MEs)[index], "\n","Median=", median.ME,". Max=", max.ME, ". Diff=", diff.ME))
}
ME.barplot<-function(MEs,index) {
  barplot(MEs[,index],col=colnames(MEs)[index],names.arg=rownames(MEs),las=2,ylim=c(-1,1),cex.names=0.6)
}

ME_plots <- function(MEs,DS.tag) {
  nmods<-ncol(MEs)
  nplots<-floor(nmods/6)
  pdf(paste("ME_plots.",DS.tag,".pdf"),width=16,height=9)
  for (i in 1:nplots) {
    par(mfrow=c(2,6))
    startIndex<-6*i-5
    endIndex<-6*i
    for (j in startIndex:endIndex) {
      ME.boxplot(MEs=MEs,index=j)
    }
    for (j in startIndex:endIndex) {
      ME.barplot(MEs=MEs,index=j)
    }
  }
  nExtraPlots<-nmods-nplots*6
  if (nExtraPlots>0) {
    par(mfrow=c(2,nExtraPlots))
    startIndex<-nplots*6+1
    endIndex<-nmods
    for (k in startIndex:endIndex) {
      ME.boxplot(MEs=MEs,index=k)
    }
    for (k in startIndex:endIndex) {
      ME.barplot(MEs=MEs,index=k)
    }
  }
  dev.off()
  graphics.off()
}