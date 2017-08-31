#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
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

#Function that computes the score of a specific bitstring
# todo: this is similar to wha is done in gaBinaryT1. need to do the same for T2
computeScoreT1<-function(CNOlist, model, bString, simList=NULL, indexList=NULL, 
                         sizeFac=0.0001, NAFac=1, timeIndex=2, title=NULL,
                         scoreT0=TRUE, initState=TRUE, multiTP=NULL){
  # simList and indexList are computed inside this function. 
  # However, for back-compatibility, we keep the arguments so that if
  # provided, we can still use them.
  
  # uncomment if required for debugging.
  # if (timeIndex<2){ stop("timeIndex must be >=2")}
  # if (timeIndex>length(CNOlist@timeSignals)){ 
  #      stop(paste("timeIndex must be <= ", length(CNOlist@timeSignals),sep=" "))
  # }
  
  if (is.null(simList)==TRUE){
    simList = prep4sim(model)
  }
  if (is.null(indexList)==TRUE){
    indexList = indexFinder(CNOlist, model)
  }
  
  nameSimIndiv <- paste(title,"_",paste(as.character(bString), collapse = ""), sep = "")
  
  # TC oct 2012. cutModel is a function. It performs the cut to select interMat, 
  # reacID, nameSpecies, and notMat. Here, we need only need the 3  and possibly only 2, so let us 
  # just copy and paste the code. This will save 10% of computational time
  modelCut <- cutModel(model, bString)
  #bs = as.logical(bString)
  #modelCut <- list()
  #modelCut$interMat <- model$interMat[, bs]
  #modelCut$reacID <- model$reacID[bs]
  #modelCut$namesSpecies <- model$namesSpecies
  
  lenTr <- dim(CNOlist@cues)[1]
  #timeMaxi <- max(CNOlist@timepoints)
  


  # ====== Writing of the .bnd file ====== #
  bndGenerator(CNOlist, modelCut, nameSimIndiv)
  print("bnd ecrit")

  # ====== Test to decide if wether does a time course study or steady state study ====== #
  # == Write .cfg files one after the other and runs the MaBoSS simulation == #
  if (multiTP == TRUE) { 
    for (x in 1:lenTr) {
      cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
      system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
      timeMaxi <- NULL
    }
  } else {
    # == Find the time of the MaBoSS simulation to get the SS through simDuration.R == #
    for (x in 1:lenTr) {
      if (exists("timeMaxi") == FALSE) {
        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
      } else {
        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv, timeMaxi)
      }
      timeMaxi <- simDuration(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
    }
  }
  

  # ====== Extract the simulated values and store them in the good format to get the score ===== #
  mbssResults <- mbssResults(CNOlist, modelCut, nameSim=nameSimIndiv, multiTP=multiTP, timeMaxi=timeMaxi)


  # ====== Remove the .bnd and .cfg files ====== #
  removal <- paste("rm -r ",nameSimIndiv,"*", sep = "")
  system(removal)
  for (afile in list.files(path = ".")){
    if (str_detect(afile, nameSimIndiv) == TRUE){
      system(paste("rm -r ", afile, sep=""))
    }
  }
  

  # ====== Former code to run the simulation using the CellNOptR deterministic simulator ====== #
  #simListCut <- cutSimList(simList, bString)
  
  
  # Compute the simulated results
  #nStimuli = length(indexList$stimulated)
  #nInhibitors <- length(indexList$inhibited)
  #nCond <- dim(CNOlist@stimuli)[1]
  #nReacs <- length(modelCut$reacID)
  #nSpecies <- length(model$namesSpecies) # this is correct. No need to get modelCut$namesSpecies 
  #nMaxInputs <- dim(simListCut$finalCube)[2]
  
  
  # simList matrices. C code must handle the matrix indices carefully.
  # This is faster than transforming into a vector as in the previous code.
  #finalCube = as.integer(simListCut$finalCube-1)
  #ixNeg = as.integer(simListCut$ixNeg)
  #ignoreCube = as.integer(simListCut$ignoreCube)
  #maxIx = as.integer(simListCut$maxIx-1)
  
  # index. convertion from R to C indices convention.
  #indexSignals <- as.integer(indexList$signals-1)
  #indexStimuli <- as.integer(indexList$stimulated-1)
  #indexInhibitors <- as.integer(indexList$inhibited-1)
  #nSignals <- length(indexSignals)
  
  
  # cnolist
  #valueInhibitors <- as.integer(CNOlist@inhibitors)
  #valueStimuli <- as.integer(CNOlist@stimuli)
  
  #simResults = .Call("simulatorT1", nStimuli, nInhibitors,
  #	nCond, nReacs, nSpecies, nSignals, nMaxInputs,
  #      finalCube, ixNeg, ignoreCube, maxIx,
  
  #      indexSignals, indexStimuli, indexInhibitors, valueInhibitors,
  #      valueStimuli, as.integer(1))
  
  #simResultsT0 = .Call("simulatorT1", nStimuli, nInhibitors,
  #      nCond, nReacs, nSpecies, nSignals, nMaxInputs,
  #      finalCube, ixNeg, ignoreCube, maxIx,
  #      indexSignals, indexStimuli, indexInhibitors, 
  #      valueInhibitors, valueStimuli, as.integer(0))
  
  
  # Unlike in computeScoreTN, the C code does not cut the simResults so you need to do it
  # since it is done in getFit, no need to do it. 
  # If we cut the simRsults in the C code, we must comment these lines AND set indexlist=NULL in the 
  # get Fit call here below. However, the C code is called elwhere in simulateTN where the expected input to
  # the C code simulateTN expect the full simResults matrix !!!
  #simResultsT0 = simResultsT0[, indexList$signals]
  #simResults = simResults[, indexList$signals]
  #nInTot = length(which(model$interMat == -1))
  
  # ====== END of the CellNOptR deterministic simulation ====== #

  
  #Compute the score
  # == Necessary if there is a time course study, can also compute the score for 2 time points == #
  # == without searching fo the steady state. This implementation is faster than the else part == #
  if (multiTP == TRUE) {
    Score <- getFit_multiTimePoints(
      simResults=mbssResults,
      CNOlist=CNOlist,
      model=modelCut,
      #indexList=indexList,
      timePoint=timeIndex,
      sizeFac=sizeFac,
      NAFac=NAFac,
      nInTot=length(which(model$interMat == -1)),
      simResultsT0=NULL)
  } else {
    # == Compute the score when there is only 2 timepoints == #
    Score <- getFit(
      simResults=mbssResults,
      #simResults=simResults, # == line from  CellNOptR == #
      CNOlist=CNOlist,
      model=modelCut,
      #indexList=indexList,
      timePoint=timeIndex,
      sizeFac=sizeFac,
      NAFac=NAFac,
      nInTot=length(which(model$interMat == -1)),
      simResultsT0=NULL)
  }
    
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  nDataP <- sum(!is.na(CNOlist@signals[[timeIndex]]))
  Score <- Score/nDataP
  
  return(Score)
}
