mabossSimulation <- function(x, nameSimIndiv, CNOlist, modelCut, timeMaxi) {
  # ====== Run the simulation with MaBoSS and test the steady state statment ====== #

  system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
  
  timeMaxi <- testOnSteadyState(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
  
  return(timeMaxi)
}