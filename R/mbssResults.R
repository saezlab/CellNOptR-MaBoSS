mbssResults <- function(CNOlist, model, nCues, nameSim=NULL, mode=1){
  mbssSim <- list()
  # recuperation of 2 lines in the simulation
  # - at time 0
  # - at the end of the simulation by default, will be improved
  # later with the dynamique time warping algorithm
  for (i in 1:nCues){
    #readMaboss <- function(x) {
    mbssSimulation <- read.table(paste(nameSim, "_", i, "/", nameSim, "_", i, "_probtraj_table.csv", sep = ""), header = TRUE)
    
    cueRow <- paste(as.character(as.numeric(CNOlistToy@cues[i,])), collapse = "")
    if (mode == 0) {
      
      mbssSim[[cueRow]] <- mbssSimulation[which(mbssSimulation$Time==0),]
      
    } else {
      
      mbssSim[[cueRow]] <- mbssSimulation[dim(mbssSimulation)[1],]
      
    }
  }
  
  
  ################################
  indexRO <- list()
  
  for (aReadOut in model$namesSpecies){
    indexRO[[aReadOut]] <- list()
    
    for (i in 1:nCues){
      aCue <- paste(as.character(as.numeric(CNOlistToy@cues[i,])), collapse = "")
      
      for (aTransState in names(mbssSim[[aCue]])){
        
        if ((str_detect(aTransState, aReadOut)==TRUE) && (str_detect(aTransState,"^Prob")==TRUE)){
          indexRO[[aReadOut]][[aCue]] <- c(indexRO[[aReadOut]][[aCue]], which(names(mbssSim[[aCue]])==aTransState))
        }
      }
    }
  }
  #############################
  
  mbssMatrix <- matrix(data = 1:(length(model$namesSpecies)*nCues),
                          nrow = nCues,
                          ncol = length(model$namesSpecies),
                          dimnames = list(names(mbssSim),model$namesSpecies)
                       )
  
  #############################
  addProbabilities <- function(x){
    mtxIndex <- which(mbssMatrix==x, arr.ind = TRUE)
    cueName <- rownames(mbssMatrix)[mtxIndex[1]]
    rdtoutName <- colnames(mbssMatrix)[mtxIndex[2]]
    rdtoutIndices <- indexRO[[rdtoutName]][[cueName]]
    
    if (length(rdtoutIndices) != 0){
      rdtvalues <- mbssSim[[cueName]][rdtoutIndices]
      x <- sum(rdtvalues)
    } else {
      x <- 0
    }
  }
  
  mbssMatrix <- apply(mbssMatrix, c(1,2), addProbabilities)
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
  
  return(mbssMatrix)
}

