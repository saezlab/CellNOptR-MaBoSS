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

computeScoreT1_mbs<-function(CNOlist, model, mbsPath, bString, simList=NULL, indexList=NULL, 
                            sizeFac=0.0001, NAFac=1, timeIndex=2, title=NULL,
                            scoreT0=TRUE, initState=TRUE, multiTP=NULL){
  
  ## simList and indexList are computed inside this function. 
  ## However, for back-compatibility, we keep the arguments so that if
  ## provided, we can still use them.
  
  tmpPath = getwd()
  setwd(mbsPath)
  system("mkdir allCond")
  
  if (is.null(simList)==TRUE){
    simList = prep4sim(model)
  }
  if (is.null(indexList)==TRUE){
    indexList = indexFinder(CNOlist, model)
  }
  
  nameSimIndiv <- paste(title,"_",paste(as.character(bString), collapse = ""), sep = "")
  
  modelCut <- cutModel(model, bString)
  lenTr <- dim(CNOlist@cues)[1]

  ##Writing of the .bnd file
  bndGenerator(CNOlist, modelCut, nameSimIndiv)
  print("bnd ecrit")

  ##Test to decide if wether does a time course study or steady state study
  ##Write .cfg files one after the other and runs the MaBoSS simulation
  if (multiTP == TRUE) { 
    for (x in 1:lenTr) {
      cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
      system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",
                   nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
      timeMaxi <- NULL
    }
  } else {
    ##Find the time of the MaBoSS simulation to get the SS through simDuration.R
    for (x in 1:lenTr) {
      if (exists("timeMaxi") == FALSE) {
        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
      } else {
        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv, timeMaxi)
      }
      timeMaxi <- simDuration(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
    }
  }

  ##Extract the simulated values and store them in the good format to get the score
  mbssResults <- mbssResults(CNOlist, modelCut, nameSim=nameSimIndiv, 
                             multiTP=multiTP, timeMaxi=timeMaxi)


  ##Remove the .bnd and .cfg files
  removal <- paste("rm -r ",nameSimIndiv,"*", sep = "")
  system(removal)
  for (afile in list.files(path = ".")){
    if (str_detect(afile, nameSimIndiv) == TRUE){
      system(paste("rm -r ", afile, sep=""))
    }
  }
  
  ##Compute the score
  ##Necessary if there is a time course study, can also compute the score for 2 time points
  ##without searching fo the steady state. This implementation is faster than the else part
  if (multiTP == TRUE) {
    Score <- getFit_multiTimePoints(
      simResults=mbssResults,
      CNOlist=CNOlist,
      model=modelCut,
      timePoint=timeIndex,
      sizeFac=sizeFac,
      NAFac=NAFac,
      nInTot=length(which(model$interMat == -1)),
      simResultsT0=NULL)
  } else {
    ##Compute the score when there is only 2 timepoints
    Score <- getFit_mbs(
      simResults=mbssResults,
      CNOlist=CNOlist,
      model=modelCut,
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
  
  setwd(tmpPath)
  
  return(Score)
}