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

bndGenerator <- function(CNOlist, model, nameSim=NULL){
  ### Does the transcription from a given topology in model to logical rules
  ### Creates a new file
  
  # ====== Creation of a new file ====== #
  fileName <- paste(nameSim,".bnd", sep = "")
  write("", file = fileName)

  # ====== Extraction of informations from the topology ====== #
  nameSpecies <- c(model$namesSpecies)
  nameReac <- c(model$reacID)
  stim <- colnames(CNOlist@stimuli)


  # ====== Write the Boolean rules according to the topology tested ====== # 
  if (length(nameReac) == 0) {
    for (spec in setdiff(nameSpecies,stim)) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write("  logic = (0);", file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    for (spec in stim) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (",spec,");", sep=""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
  # ====== Case when there is only one reaction in the topology ====== #
  else if (length(nameReac) == 1) {
    target <- str_split(model$reacID[1], "=", simplify = TRUE)[2]
    reg <- str_split(model$reacID[1], "=", simplify = TRUE)[1]
    
    write(paste("Node", target, "{", sep = " "), file = fileName, append = TRUE)
    write("  logic = ", file = fileName, append=TRUE)
    

    # ====== Write the Boolean rule for the regulated node ====== #
    # == target regulated by only one node == #
    if (str_detect(reg, "[+]") == FALSE){
      write(paste("(",reg,");", sep=""), file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)

    } else {
      # == Rule with an AND gate == #

      reg <- str_split(reg,"[+]", simplify = TRUE)
      firstElt = TRUE
      write("(", file = fileName, append = TRUE)
      for (aReg in reg) {
        ##### Add of AND condition in case of necessity #####
        if (firstElt != TRUE){
          write(" & ", file = fileName, append = TRUE)
        } else {
          firstElt = FALSE
        }
        write(aReg, file = fileName, append = TRUE)
      }
      write(");", file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    

    # ====== Write the rule for all the other nodes : stimuli and not targeted nodes ====== #
    nameSpecies <- setdiff(nameSpecies,target)
    
    for (spec in stim){
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (", spec, ");", sep = ""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }

    for (spec in setdiff(nameSpecies,union(target,stim))) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (0);", sep = ""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
        file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
        file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }


  # ====== Case when there are more than 1 reaction ====== #
  else {
    for (target in nameSpecies){
      write(paste("Node", target, "{", sep = " "), file = fileName, append = TRUE)
      write("  logic = ", file = fileName, append=TRUE)
      
      # ====== Construction of the boolean rules ====== #
      bRule = FALSE # == indicates if the node is regulated by another or not
      firstReg = TRUE
      for (reac in nameReac){
        if (model$interMat[target,reac][1] == 1){
          bRule = TRUE

          # ====== The node is regulated and requires a Boolean rule ====== #
          # == Add of OR condition in case of necessity == #
          if (firstReg == FALSE){
            write(" | ", file = fileName, append = TRUE)
          } else {
            firstReg = FALSE
          }
          
          reg <- str_split(reac,"=", simplify = TRUE)[1]
          if (str_detect(reg, "[+]") == FALSE){ ### simple edge A --> B
            write(paste("(",reg,")", sep=""), file = fileName, append = TRUE)
          
            } else { ### hyperedge A + C --> B
            reg <- str_split(reg,"[+]", simplify = TRUE)
            firstElt = TRUE
            write("(", file = fileName, append = TRUE)
            for (aReg in reg) {
              ##### Add of AND condition in case of necessity #####
              if (firstElt != TRUE){
                write(" & ", file = fileName, append = TRUE)
              } else {
                firstElt = FALSE
              }
              write(aReg, file = fileName, append = TRUE)
            }
            write(")", file = fileName, append = TRUE)
            }  
        }

        # ====== the node is never targeted, so not regulated ====== #
        # == might be a stimuli or a simple node == #
        else if ((bRule == FALSE) && (reac == tail(nameReac,1))){
          if (target %in% colnames(CNOlist@stimuli)) {
            write(paste("(",target,")", sep = ""), file = fileName, append = TRUE)
          } else {
            write("(0)", file = fileName, append = TRUE)
          }
        }
      }
      write(";", file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
  # ====== End of the creation of the Boolean rules ====== #

  # ====== Combine the lines in order to write the file as it is required by MaBoSS ====== #
  myFile <- readLines(fileName)
  write("", file=fileName)
  for (wd in myFile){
    if (str_detect(wd, "^$") != TRUE){
      if (str_detect(wd, "^Node") == TRUE){ ### Beginning of the rule
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "^  logic") == TRUE){ ### Write the boolean rule
        lgcl <- wd
      } else if (str_detect(wd, "^  rate_up") == TRUE){ ### Write the rules regarding rates
        write(lgcl, file = fileName, append = TRUE)
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "^  rate_down") == TRUE){ ### kind of useless because written in the previous step
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "[}]") == TRUE){ ### End of the rule
        write(paste(wd,"\n", sep = ""), file = fileName, append = TRUE)
      } else {
        lgcl <- paste(lgcl, wd, sep="")
      }
    }
  }
}