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

cSimulator <- function(CNOlist, model, simList, indexList, mode=1,scoreT0=TRUE, initState=TRUE) {
	lenTr <- dim(CNOlist@cues)[1]
	nameSimIndiv <- "best_solution"
	bndGenerator(CNOlist, model, nameSimIndiv)

	simDuration <- function(x, nameSimIndiv, CNOlist, model, timeMaxi){
    	simDuration <- mabossSimulation(x, nameSimIndiv, CNOlist, model, timeMaxi)
    	if (timeMaxi != simDuration) {
    		timeMaxi <- simDuration
      
    		### Remove the current simulation
    		nameFolder <- paste(nameSimIndiv,x,sep="_")
    		system(paste("rm -r ",nameFolder, "*", sep = ""))
     		for (afile in list.files(path = ".")){
        		if (str_detect(afile, nameFolder) == TRUE){
          			system(paste("rm -r ", afile, sep=""))
        		}
      		}
      
    		### Recursive step
    		cfgGenerator(CNOlist, model, x, nameSimIndiv, timeMaxi, initState)
    		simDuration(x, nameSimIndiv, CNOlist, model, timeMaxi)
    	} else {
    		return(timeMaxi)
    	}
  	}
  
  
  
	for (x in 1:lenTr) {
    	if (exists("timeMaxi") == FALSE) {
    		timeMaxi <- cfgGenerator(CNOlist, model, x, nameSimIndiv)
    	} else {
    		timeMaxi <- cfgGenerator(CNOlist, model, x, nameSimIndiv, timeMaxi=timeMaxi)
    	}
    	timeMaxi <- simDuration(x, nameSimIndiv, CNOlist, model, timeMaxi=timeMaxi)
  	}
  	#print(timeMaxi)
  	#print(mode)
  	#if (mode == 0) {
	  	res <- mbssResults(CNOlist, model, nameSim=nameSimIndiv)
	#} else if (mode == 1) {
	# 	res <- mbssResults(CNOlist, model, nameSimIndiv, mode=1)
	#}

	removal <- paste("rm -r ",nameSimIndiv,"*", sep = "")
  	system(removal)
  	for (afile in list.files(path = ".")){
    	if (str_detect(afile, nameSimIndiv) == TRUE){
      		system(paste("rm -r ", afile, sep=""))
    	}
  	}

  	#print(res)
	# variables
	#nStimuli <- as.integer(length(indexList$stimulated))
	#nInhibitors <- as.integer(length(indexList$inhibited))
	#nCond <- as.integer(dim(CNOlist@stimuli)[1])
	#nReacs <- as.integer(length(model$reacID))
	#nSpecies <- as.integer(length(model$namesSpecies))
	#nMaxInputs <- as.integer(dim(simList$finalCube)[2])
	
	# simList
	# used to be 
    # >>> finalCube = as.integer(as.vector(t(simList$finalCube))-1)
    # but as.vector(t is slow and can be replaced by just as.integer albeit
    # appropriate C modifications
	#finalCube = as.integer(simList$finalCube-1)
	#ixNeg = as.integer(simList$ixNeg)
	#ignoreCube = as.integer(simList$ignoreCube)
	#maxIx = as.integer(simList$maxIx-1)
	
	# index
	#indexSignals <- as.integer(indexList$signals-1)
	#indexStimuli <- as.integer(indexList$stimulated-1)
	#indexInhibitors <- as.integer(indexList$inhibited-1)
    #nSignals <- length(indexSignals)

	# cnolist
	#valueInhibitors <- as.integer(CNOlist@inhibitors)
	#valueStimuli <- as.integer(CNOlist@stimuli)

	#res = .Call("simulatorT1",
		# variables	
	#	nStimuli,
	#	nInhibitors,
	#	nCond,
	#	nReacs,
	#	nSpecies,
    #   nSignals,
	#	nMaxInputs,
	#	# simList
	#	finalCube,
	#	ixNeg,
	#	ignoreCube,
	#	maxIx,		
	#	# index
	#	indexSignals, 
	#	indexStimuli, 
	#	indexInhibitors,
	#	# cnolist
	#	valueInhibitors,
	#	valueStimuli,
    #    as.integer(mode)
	#)
# should not be cut because it is used in simulateTN as an input 
#    res = res[,indexList$signals]
	return(res)
}
