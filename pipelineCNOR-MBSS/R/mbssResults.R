mbssResults <- function(CNOlist, model, nameSim=NULL){
  mbssSim <- list()
  # recuperation of 2 lines in the simulation
  # - at time 0
  # - at the end of the simulation by default, will be improved
  # later with the dynamique time warping algorithm

  nCues <- dim(CNOlist@cues)[1]
  timePoints <- CNOlist@timepoints

  for (i in 1:nCues){
    #readMaboss <- function(x) {
    mbssSimulation <- read.table(paste(nameSim, "_", i, "/", nameSim, "_", i, "_probtraj_table.csv", sep = ""), header = TRUE)
    #print(colnames(mbssSimulation))

    #### index of the columns of interest
    colProb <- c()
    for (aColName in colnames(mbssSimulation)) {
      #print(aColName)
      if ((str_detect(aColName,"^Prob")==TRUE) || (str_detect(aColName, "^Time")==TRUE)) {
        colProb <- append(colProb, which(colnames(mbssSimulation) == aColName,arr.ind=TRUE))
      }
    }

    #### Selection of the lines
    cueRow <- paste(as.character(as.numeric(CNOlist@cues[i,])), collapse = "")
    mbssSim[[cueRow]] = list()
    for (aTP in timePoints) {
      mbssSim[[cueRow]][[as.character(aTP)]] <- mbssSimulation[which(mbssSimulation$Time==aTP),colProb]

      for (aColName in colnames(mbssSim[[cueRow]][[as.character(aTP)]])) {
        if (str_detect(aColName, "^Prob") == TRUE) {
          colnames(mbssSim[[cueRow]][[as.character(aTP)]]) [which(aColName == colnames(mbssSim[[cueRow]][[as.character(aTP)]]), arr.ind=TRUE)] <- str_replace(aColName,"^Prob", "")
        }
      }


    }
    #print(mbssSim[[cueRow]])
  }


  ################################
  indexRO <- list()
  for (aReadOut in colnames(CNOlistToy@signals$`0`)) {
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
  #print(indexRO)
  #############################
  
  mbssMatrix <- list()
  for (aTP in timePoints) {
    tp <- as.character(aTP)
  	mbssMatrix[[tp]] <- matrix(data = seq(1:(length(colnames(CNOlistToy@signals$`0`))*nCues))+1,
                          nrow = nCues,
                          ncol = length(colnames(CNOlistToy@signals$`0`)),
                          dimnames = list(names(mbssSim),colnames(CNOlistToy@signals$`0`))
    )
  }
  
  #############################
  addProbabilities <- function(x, mat=NULL, tp=NULL){
    #print(paste(x, tp, sep="    "))
    mtxIndex <- which(mat==x, arr.ind = TRUE)
    #print(mtxIndex)
    cueName <- rownames(mat)[mtxIndex[1]]
    #print(cueName)
    rdtoutName <- colnames(mat)[mtxIndex[2]]
    #print(rdtoutName)
    rdtoutIndices <- indexRO[[rdtoutName]][[cueName]]
    #print(rdtoutIndices)
    

    if (length(rdtoutIndices) != 0){
      mat[cueName,rdtoutName] <- sum(mbssSim[[cueName]][[tp]][rdtoutIndices])
    } else {
      mat[cueName,rdtoutName] <- 0
    }
    #print(mat[cueName,rdtoutName])
  }


  for (aTP in timePoints) {
    tp <- as.character(aTP)
    #print(is.matrix(mbssMatrix[[tp]]))
    #print(dim(mbssMatrix[[tp]]))
    #print(rownames(mbssMatrix[[tp]]))
    #print(colnames(mbssMatrix[[tp]]))
    print(tp)
    #print(mbssMatrix[[tp]])
    mbssMatrix[[tp]] <- apply(mbssMatrix[[tp]], c(1,2), addProbabilities, mat=mbssMatrix[[tp]], tp=tp)
    print(mbssMatrix[[tp]])
    #print("test lalala")
  }
  ############################
  
  
  
  # for loop in namesSpecies
  #for (aCue in names(mbssSim)) {
  #  for (aName in model$namesSpecies) {
  #    simValue <- 0
  #    for  (aTitle in colnames(mbssSim[[aCue]])){
  #      if ((str_detect(aTitle, "^Prob") == TRUE) & (str_detect(aTitle, aName) == TRUE)){
  #        simValue <- simValue + mbssSim[[aCue]][,aTitle]
  #      }
  #    }
  #    mbssMatrix[aCue,aName] <- simValue
  #  }
  #}
  # detection of prob and name of a species
  # addition of the values
  # /!\ maybe a problem with the NAs
  #print(mbssMatrix)
  return(mbssMatrix)
}

