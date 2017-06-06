cfgGenerator <- function (CNOlist, modelCut, treatmt, nameSim=NULL, timeMaxi=NULL, initState=TRUE){
  
  if (is.null(timeMaxi) == TRUE){
    timeMaxi <- max(CNOlist@timepoints)
  }
  
  print(paste("create the file : ",nameSim,"_",treatmt,".cfg", sep=""))
  
  dest_file <- paste(nameSim,"_",treatmt,".cfg", sep="")
  file.copy("basicFile.cfg", dest_file, overwrite = TRUE)
  
  #### Max time of the simulation
  runTime <- paste("max_time = ", timeMaxi, ";", sep="")
  write(runTime, file = dest_file, append = TRUE)
  
  all_spec <- modelCut$namesSpecies
  stim <- colnames(CNOlist@stimuli)
  inhib <- colnames(CNOlist@inhibitors)
  
  ###################
  #### the RATES ####
  ###################
  
  #### All the nodes less those that can be inhibited have their rate set to 1
  for (node in setdiff(all_spec,inhib)){
    upRate <- paste("$u_",node," = 1;", sep="")
    write(upRate, file = dest_file, append = TRUE)
    downRate <- paste("$d_", node, " = 1;", sep="")
    write(downRate, file = dest_file, append = TRUE)
  }
  
  #### Particular case of the inhibited nodes
  ## if the treatment indicates an inhibition, the rates are ste to 0
  ## else set to 1
  for(anInh in inhib) {
    if (CNOlist@inhibitors[treatmt, anInh] == 1){
      upRate <- paste("$u_",anInh," = 0;", sep="")
      write(upRate, file = dest_file, append = TRUE)
      downRate <- paste("$d_", anInh, " = 0;", sep="")
      write(downRate, file = dest_file, append = TRUE)
    } else {
      upRate <- paste("$u_",anInh," = 1;", sep="")
      write(upRate, file = dest_file, append = TRUE)
      downRate <- paste("$d_", anInh, " = 1;", sep="")
      write(downRate, file = dest_file, append = TRUE)
    }
  }
  
  ############################
  #### the Initial States ####
  ############################
  
  #### the stimulated nodes are described in CNOlist
  for (aStim in stim){
    istate <- paste(aStim, ".istate = ", CNOlist@stimuli[treatmt,aStim], ";", sep="")
    write(istate, file = dest_file, append = TRUE)
  }
  
  #### all the nodes that have not inhibitors or stimuli are initialized to 0
  if (initState){
    nodes <- setdiff(all_spec, stim)
    for (aNode in nodes){
      istate <- paste(aNode, ".istate = 0;", sep="")
      write(istate, file = dest_file, append = TRUE)
    }
  } #### only the inhibited nodes will be initialized to 0
  else {
    #### the inhibited nodes are set to 0
    for (anInh in inhib){
      if (CNOlist@inhibitors[treatmt, anInh] == 1){
        istate <- paste(anInh, ".istate = 0;", sep="")
        write(istate, file = dest_file, append = TRUE)
      }
    }
  }
  
  return(timeMaxi)
}