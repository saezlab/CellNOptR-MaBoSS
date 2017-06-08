bndGenerator <- function(CNOlist, model, nameSim=NULL){
  ### Does the transcription from a given topology in model to logical rules
  ### Creates a new file
  
  fileName <- paste(nameSim,".bnd", sep = "")
  write("", file = fileName)
  nameSpecies <- c(model$namesSpecies)
  nameReac <- c(model$reacID)
  
  # May be useless, CNOR doesn't allow models with 1 reaction
  if (length(nameReac) == 1){
    target <- str_split(model$reacID, "=", simplify = TRUE)[2]
    reg <- str_split(model$reacID, "=", simplify = TRUE)[1]
    
    write(paste("Node", target, "{", sep = " "), file = fileName, append = TRUE)
    write("  logic = ", file = fileName, append=TRUE)
    
    if (str_detect(reg, "[+]") == FALSE){
      write(paste("(",reg,");", sep=""), file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    } else {
      reg <- str_split(reg,"+", simplify = TRUE)
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
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    
    nameSpecies <- setdiff(nameSpecies,target)
    
    for (spec in nameSpecies){
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (", spec, ");", sep = ""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  } else {
    for (spec in nameSpecies){
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write("  logic = ", file = fileName, append=TRUE)
      
      ##### Construction of the boolean rules #####
      bRule = FALSE
      firstReg = TRUE
      for (reac in nameReac){
        if (model$interMat[spec,reac][1] == 1){
          bRule = TRUE
          ##### Add of OR condition in case of necessity #####
          if (firstReg != TRUE){
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
          
          ### I doubt on this test
        } else if ((bRule == FALSE) && (reac == tail(nameReac,1))){
          if (spec %in% colnames(CNOlist@stimuli)) {
            write(paste("(",spec,")", sep = ""), file = fileName, append = TRUE)
          } else {
            write("(0)", file = fileName, append = TRUE)
          }
        }
      }
      write(";", file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
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