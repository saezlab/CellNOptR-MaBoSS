mabossSimulation <- function(x, nameSimIndiv, CNOlist, modelCut, timeMaxi) {
  system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
  #system(paste("perl ./tools/MBSS_TrajectoryFig.py ",nameSimIndiv,"_",x, sep = ""))
  
  timeMaxi <- testOnSteadyState(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
  
  return(timeMaxi)
}