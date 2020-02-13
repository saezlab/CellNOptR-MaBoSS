#
#  This file is part of the CNO software
#
#  Copyright (c) 2020 - Heidelberg University
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

cfgGenerator <- function (CNOlist, modelCut, treatmt, 
                          nameSim=NULL, timeMaxi=NULL, initState=TRUE){
  ### Write the cfg file
  ### Initial conditions depends on the treatment
  
  if (is.null(timeMaxi) == TRUE){
    timeMaxi <- max(CNOlist@timepoints)
  }
  
  ##Create the new file copying a basic file present in the directory
  dest_file <- paste(nameSim,"_",treatmt,".cfg", sep="")
  file.copy("basicFile.cfg", dest_file, overwrite = TRUE)
  
  ##Max time of the simulation
  runTime <- paste("max_time = ", timeMaxi+0.10, ";", sep="")
  write(runTime, file = dest_file, append = TRUE)
  
  all_spec <- modelCut$namesSpecies
  stim <- colnames(CNOlist@stimuli)
  inhib <- colnames(CNOlist@inhibitors)
  measured <- colnames(CNOlist@signals$`0`)
  
  ##All the nodes less those that can be inhibited have their rate set to 1
  for (node in setdiff(all_spec,inhib)){
    upRate <- paste("$u_",node," = 1;", sep="")
    write(upRate, file = dest_file, append = TRUE)
    downRate <- paste("$d_", node, " = 1;", sep="")
    write(downRate, file = dest_file, append = TRUE)
  }
  
  ##Particular case of the inhibited nodes
  ##if the treatment indicates an inhibition, the rates are set to 0
  ##else set to 1
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

  ##the stimulated nodes (inputs) are described in CNOlist
  ##the value 1 means an activation, the value 0 means that the input is not stimulated
  for (aStim in stim){
    istate <- paste(aStim, ".istate = ", CNOlist@stimuli[treatmt,aStim], ";", sep="")
    write(istate, file = dest_file, append = TRUE)
  }
  
  ##part that "switch" the value between 0 or 1
  ##regarding the experimental data at t=0
  for (aNode in measured) {
    if (CNOlist@signals$`0`[treatmt,aNode] > 0.5) {
      istate <- paste(aNode, ".istate = 1;", sep="")
      write(istate, file = dest_file, append = TRUE)
    } else {
      istate <- paste(aNode, ".istate = 0;", sep="")
      write(istate, file = dest_file, append = TRUE)
    }
  }


  ##initial condition for all the other nodes
  if (initState){
    ##all the nodes that have not inhibitors or stimuli are initialized to 0
    
    nodes <- setdiff(all_spec, union(stim, measured))
    for (aNode in nodes){
      istate <- paste(aNode, ".istate = 0;", sep="")
      write(istate, file = dest_file, append = TRUE)
    }
  } else {
    ##only the inhibited nodes will be initialized to 0
    
    for (anInh in inhib){
      if (CNOlist@inhibitors[treatmt, anInh] == 1){
        istate <- paste(anInh, ".istate = 0;", sep="")
        write(istate, file = dest_file, append = TRUE)
      }
    }
  }
  
  return(timeMaxi)
}