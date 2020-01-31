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

densityNodesMBSS <- function(dataMB, nameFolder, timeMaxi, species, qtty, origines) {
  myPath <- getwd()
  timeSeq <- dataMB[,which(colnames(dataMB) == "Time")]
  pdf(paste(myPath,"/allCond/",nameFolder, "_", timeMaxi,".pdf", sep=""))
  par(mfrow=c(1,1))
  vecCol <- rainbow(length(species)-1)
  plot(timeSeq, unlist(qtty[[1]]), type="l", col=1, ylim = c(0,1),
       main = paste(nameFolder,"_", timeMaxi, sep=""))
  for (i in 2:length(species)){
    spec <- species[i]
    lines(timeSeq, unlist(qtty[[spec]]), col=vecCol[i-1])
    abline(origines[i],slopeValues[i], col=vecCol[i-1], lty="dotdash")
  }
  abline(v=c(timeSeq[length(timeSeq)-20]), col="black")
  legend("right", legend = names(qtty), col = c(1,vecCol), lwd = 1)
  dev.off()
}