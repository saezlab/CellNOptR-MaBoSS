#
#  This file is part of the CNO software
#
#  Copyright (c) 2020 - BioQuant Zentrum - Heidelberg University
#
#  File author(s): C. Chevalier, A. Dugourd, E. Gjerga
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

mbssResults <- function(CNOlist, model, nameSim=NULL, multiTP=NULL, timeMaxi=NULL){
  # ====== Extract the results of the simulations ====== #
  # == and make it in the same format as the experimental values are stored == #


  mbssSim <- list()

  # ========================================================= #
  # ====== recuperation of the lines in the simulation ====== #
  # ========================================================= #
  nCues <- dim(CNOlist@cues)[1]

  if (multiTP == TRUE) {
    timePoints <- CNOlist@timepoints  
  } else {
    timePoints <- c(1,2)
  }


  for (i in 1:nCues){
    mbssSimulation <- read.table(paste(nameSim, "_", i, "/", nameSim, "_", i, "_probtraj_table.csv", sep = ""), header = TRUE)

    # ====== Keep the columns of interest ====== #
    colProb <- c()
    for (aColName in colnames(mbssSimulation)) {
      #print(aColName)
      if ((str_detect(aColName,"^Prob")==TRUE) || (str_detect(aColName, "^Time")==TRUE)) {
        colProb <- append(colProb, which(colnames(mbssSimulation) == aColName,arr.ind=TRUE))
      }
    }

    # ====== Selection of the lines ====== #
    cueRow <- paste(as.character(as.numeric(CNOlist@cues[i,])), collapse = "")
    mbssSim[[cueRow]] = list()
    if (multiTP == TRUE) {
      for (aTP in timePoints) {
        mbssSim[[cueRow]][[as.character(aTP)]] <- mbssSimulation[which(mbssSimulation$Time==aTP),colProb]

        for (aColName in colnames(mbssSim[[cueRow]][[as.character(aTP)]])) {
          if (str_detect(aColName, "^Prob") == TRUE) {
            colnames(mbssSim[[cueRow]][[as.character(aTP)]]) [which(aColName == colnames(mbssSim[[cueRow]][[as.character(aTP)]]), arr.ind=TRUE)] <- str_replace(aColName,"^Prob", "")
          }
        }
      }
    } else {
      mbssSim[[cueRow]][["1"]] <- mbssSimulation[which(mbssSimulation$Time == 0),colProb]

      for (aColName in colnames(mbssSim[[cueRow]][["1"]])) {
        if (str_detect(aColName, "^Prob") == TRUE) {
          colnames(mbssSim[[cueRow]][["1"]]) [which(aColName == colnames(mbssSim[[cueRow]][["1"]]), arr.ind=TRUE)] <- str_replace(aColName,"^Prob", "")
        }
      }

      print(which(mbssSimulation$Time == max(mbssSimulation$Time)))
      mbssSim[[cueRow]][["2"]] <- mbssSimulation[which(mbssSimulation$Time == max(mbssSimulation$Time)),colProb]
      for (aColName in colnames(mbssSim[[cueRow]][["2"]])) {
        if (str_detect(aColName, "^Prob") == TRUE) {
          colnames(mbssSim[[cueRow]][["2"]]) [which(aColName == colnames(mbssSim[[cueRow]][["2"]]), arr.ind=TRUE)] <- str_replace(aColName,"^Prob", "")
        }
      }
    }
  }


  # ============================================== #
  # ====== Index of the columns of interest ====== #
  # ============================================== #
  indexRO <- list()
  for (aReadOut in colnames(CNOlist@signals$`0`)) {
    indexRO[[aReadOut]] <- list()
    #print(aReadOut)
    for (i in 1:nCues){
      aCue <- paste(as.character(as.numeric(CNOlist@cues[i,])), collapse = "")
      indexRO[[aReadOut]][[aCue]] <- c()
      for (aTransState in names(mbssSim[[aCue]][[1]])){
        if((str_detect(aTransState, aReadOut)==TRUE) && (str_detect(aTransState,"^Time") == FALSE)) {
          #print(aTransState)
          indexRO[[aReadOut]][[aCue]] <- c(indexRO[[aReadOut]][[aCue]], which(names(mbssSim[[aCue]][[1]])==aTransState))
        }
      }
    }
  }
  # ============================================== #
  

  # ====== Create the matrix where the simulated values will be stored ====== #
  mbssMatrix <- list()
  for (aTP in timePoints) {
    tp <- as.character(aTP)
  	mbssMatrix[[tp]] <- matrix(data = seq(1:(length(colnames(CNOlistToy@signals$`0`))*nCues))+1,
                          nrow = nCues,
                          ncol = length(colnames(CNOlistToy@signals$`0`)),
                          dimnames = list(names(mbssSim),colnames(CNOlistToy@signals$`0`))
    )
  }
  
  # ====== Calculates the probability of each node to be activated ====== #
  addProbabilities <- function(x, mat=NULL, tp=NULL){
    mtxIndex <- which(mat==x, arr.ind = TRUE)
    cueName <- rownames(mat)[mtxIndex[1]]
    rdtoutName <- colnames(mat)[mtxIndex[2]]
    rdtoutIndices <- indexRO[[rdtoutName]][[cueName]]
    
    if (length(rdtoutIndices) != 0){
      mat[cueName,rdtoutName] <- sum(mbssSim[[cueName]][[tp]][rdtoutIndices])
    } else {
      mat[cueName,rdtoutName] <- 0
    }
  }


  for (aTP in timePoints) {
    tp <- as.character(aTP)
    mbssMatrix[[tp]] <- apply(mbssMatrix[[tp]], c(1,2), addProbabilities, mat=mbssMatrix[[tp]], tp=tp)
  }
  
  return(mbssMatrix)
}

