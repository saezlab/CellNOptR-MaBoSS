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

cSimulator <- function(CNOlist, model, simList, indexList, mode=1,scoreT0=TRUE, initState=TRUE, multiTP=NULL) {

	if (is.null(multiTP) == TRUE) {
    	if (length(CNOlist@timepoints) > 2) {
    	  multiTP = TRUE
    	} else {
    	  multiTP = FALSE
    	}
  	}

	lenTr <- dim(CNOlist@cues)[1]
	nameSimIndiv <- "best_solution"
	bndGenerator(CNOlist, model, nameSimIndiv)
  
  
	if (multiTP == TRUE) {
	    for (x in 1:lenTr) {
	      	cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
	      	system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
	    }
  	} else {
	    for (x in 1:lenTr) {
	      if (exists("timeMaxi") == FALSE) {
	        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv)
	        print(timeMaxi)
	      } else {
	        timeMaxi <- cfgGenerator(CNOlist, modelCut, x, nameSimIndiv, timeMaxi)
	      }
	      timeMaxi <- simDuration(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
	    }
	}

  	#print(timeMaxi)
  	#print(mode)
  	#if (mode == 0) {
	  	res <- mbssResults(CNOlist, model, nameSim=nameSimIndiv, multiTP=multiTP)
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
	print("print res of cSimulator")
	print(res)
	return(res)
}
