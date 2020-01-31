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

getFit_multiTimePoints<-function(
    simResults,
    CNOlist,
    model,
    timePoint=NULL,
    sizeFac=0.0001,
    NAFac=1,
    nInTot,
    simResultsT0=NULL
    ){

    if ((class(CNOlist)=="CNOlist")==FALSE){
         CNOlist = CellNOptR::CNOlist(CNOlist)
     }

    timePoints <- CNOlist@timepoints
    
    ##Calculate the SSE
    for (aTime in timePoints) {
    	tPt <- match(aTime,timePoints)
     	Diff <- simResults[[as.character(aTime)]] - CNOlist@signals[[tPt]]
    	sqrd <- Diff^2
    	if (aTime == 0) {
    		r <- sqrd
    	} else {
    		r <- rbind(r, sqrd)
    	}
    }

    ##Penalties
    deviationPen <- sum(r[!is.na(r)])/(length(timePoints))
    NAPen<-NAFac*length(which(is.na(simResults)))
    nDataPts<-dim(CNOlist@signals[[tPt]])[1]*dim(CNOlist@signals[[tPt]])[2]
    nInputs<-length(which(model$interMat == -1))

    ## nInTot: number of inputs of expanded model
    ## nInputs: number of inputs of cut model
    sizePen<-(nDataPts*sizeFac*nInputs)/nInTot
    score<-deviationPen+NAPen+sizePen
    
    print(paste("MSE : ", deviationPen/sum(!is.na(CNOlist@signals[[tPt]])), sep=""))
    print(paste("NA penalty : ", NAPen, sep=""))
    print(paste("Size penalty : ", sizePen, sep=""))
    return(score)

}

