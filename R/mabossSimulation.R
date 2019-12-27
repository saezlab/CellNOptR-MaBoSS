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

mabossSimulation <- function(x, nameSimIndiv, CNOlist, modelCut, timeMaxi) {
  # ====== Run the simulation with MaBoSS and test the steady state statment ====== #

  system(paste("source ./MaBoSS.env ; perl ./tools/MBSS_FormatTable.pl ",nameSimIndiv,".bnd ",nameSimIndiv,"_",x,".cfg", sep=""))
  
  timeMaxi <- testOnSteadyState(x, nameSimIndiv, CNOlist, modelCut, timeMaxi)
  
  return(timeMaxi)
}