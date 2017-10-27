library(Hmisc)

plotModulesIM<-function(moduleName,moduleColors,kIM) {
  module<-moduleName
  col<-module
  moduleGenesIndex<-which(moduleColors==module)
  moduleGenes.kIM<-kIM[moduleGenesIndex,"kWithin"]
  moduleGenes.sigTrait<-gene.corr.r_p[moduleGenesIndex,sigTrait]
  if (module=="ivory") {col="gray"} ###If you have very light colors, substitute them for a darker color to be seen on graphs
  if (module=="floralwhite") {col="gray"}
  if (module=="lightcyan") {col="cyan"}
  if (module=="lightcyan1") {col="cyan"}
  if (module=="honeydew1") {col="gray"}
  if (module=="white") {col="gray"}
  if (module=="lightyellow") {col="yellow"}
  if (module=="aliceblue") {col="blue"}
  if (module=="antiquewhite1") {col="gray"}
  if (module=="lavenderblush") {col="gray"}
  if (module=="lavenderblush1") {col="gray"}
  if (module=="lavenderblush2") {col="gray"}
  if (module=="lavenderblush3") {col="gray"}
  if (module=="mistyrose") {col="gray"}
  if (module=="moccasin") {col="gray"}
  if (module=="navajowhite") {col="gray"}
  if (module=="navajowhite1") {col="gray"}
  if (module=="navajowhite2") {col="gray"}
  if (module=="paleturquoise") {col="turquoise"}
  if (module=="thistle1") {col="gray"}
  if (module=="thistle2") {col="gray"}
  if (module=="whitesmoke") {col="gray"}
  
  verboseScatterplot(moduleGenes.kIM,moduleGenes.sigTrait,
       xlab = paste0("Probe intramodular membership in ", module, " module"),
       ylab = paste0("Probe correlation to ",sigTrait),
       main = paste(capitalize(module), " kIM vs. gene significance\n", sep=""),
       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col, ylim=c(-1,1), xlim=c(-1,1))
  abline(h=0,lty=2,lwd=1,col="red")
  abline(v=0,lty=2,lwd=1,col="red")
}

kIMvsCorrPlot<-function(moduleCor,moduleColors,kIM,DS.tag) {
  
  colors.mods <- rownames(moduleCor)
  colors.mods <- colors.mods[colors.mods!="grey"]
  colors.mods <- sort(colors.mods)
  
  nmods<-length(colors.mods)
  nplots<-floor(nmods/6)
  pdf(paste0("kIMvsCorrPlots",DS.tag,".pdf"),width=16,height=9)
  for (i in 1:nplots) {
    par(mfrow=c(2,3))
    startIndex<-6*i-5
    endIndex<-6*i
    for(j in startIndex:endIndex) {
      myModName<-colors.mods[j]
      plotModulesIM(moduleName=myModName,moduleColors=moduleColors,kIM=kIM)
    }
  }
  nExtraPlots<-nmods-nplots*6
  if (nExtraPlots>0) {
    if (nExtraPlots==5) {
      par(mfrow=c(2,3))
    } else {
      par(mfrow=c(2,2))
    }
    startIndex<-nplots*6+1
    endIndex<-nmods
    for (k in startIndex:endIndex) {
      myModName<-colors.mods[k]
      plotModulesIM(moduleName=myModName,moduleColors=moduleColors,kIM=kIM)
    }
  }
  dev.off()
  graphics.off()
}

