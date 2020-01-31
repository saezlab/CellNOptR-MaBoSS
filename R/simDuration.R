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

simDuration <- function(x, nameSimIndiv, CNOlist, model, timeMaxi){
  
  ##Recursive function necessaryt to run MaBoSS simulations until the steady state
  simDuration <- mabossSimulation(x, nameSimIndiv, CNOlist, model, timeMaxi)
  if (timeMaxi != simDuration) {
    timeMaxi <- simDuration
      
    ##Remove the current simulation
    nameFolder <- paste(nameSimIndiv,x,sep="_")
    system(paste("rm -r ",nameFolder, "*", sep = ""))
    for (afile in list.files(path = ".")){
      if (str_detect(afile, nameFolder) == TRUE){
        system(paste("rm -r ", afile, sep=""))
      }
    }
      
    ##Rewrite the .cfg file with additional time
    cfgGenerator(CNOlist, model, x, nameSimIndiv, timeMaxi)
    simDuration(x, nameSimIndiv, CNOlist, model, timeMaxi)
  } else {
    return(timeMaxi)
  }
  
}