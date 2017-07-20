#!/usr/bin/env zsh

rm(list=c(ls()))

#setwd(dir = "/Users/celine/MaBoSS-env-2.0")

#source("https://bioconductor.org/biocLite.R")
#biocLite("doParallel")

#library(CNORode2017)
library(CellNOptR)
library(stringr)
#library(parallel)
#library(doParallel)

#registerDoParallel(cores=(detectCores(all.tests = FALSE, logical = TRUE)-1))
#getDoParWorkers()

#system("mkdir probPlots")
#system("mkdir allCond")

###### Working with arguments
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) == 0){
#  stop("At least one argument must be supplied (name of simulation)", call. = FALSE)
#} else if (length(args) > 1) {
#  stop("Only one argument reauired. Please choose a name for the simulation", call. = FALSE)
#}
#nameSim <- args[1]
#system(paste("rm -r ",nameSim,"*",sep = ""))
####### the call to gaBinaryT1 function changes too

####################################################
##  Loading the data and prior knowledge network  ##
####################################################

# experimental data : provides studied species and timepoints
# there is also a list of treatments that does not appear by the simple
# command "> CNOlistToy"

#CNOlistToy = CNOlist("ToyDataMMB.csv")
#CNOlistToy=CNOlist(system.file("doc", "ToyModelMMB_FeedbackAnd.csv",
#                               package="CNORode"))
CNOlistToy=CNOlist("/Users/celine/modelMacNamara2012/ToyModelPBode.csv");
#CNOlistToy=CNOlist("/Users/celine/model2nodes/model2nodes.csv");
#CNOlistToy=CNOlist(system.file("doc", "ToyModelMMB_FeedbackAnd.csv",
#                               package="CNORode"))
#CNOlistToy=CNOlist("/Users/celine/modelMacNamara2012/ToyModelPB.csv");
#CNOlistToy
# Another way to visualize the data and export the plot on a .pdf format
#plot(CNOlistToy)
#plotCNOlistPDF(CNOlist=CNOlistToy,filename="ToyModelGraph.pdf")


# PKN model obtained from CytoScape
#pknmodel<-readSIF("ToyPKNMMB.sif")
#model=readSIF(system.file("doc", "ToyModelMMB_FeedbackAnd.sif",
#         package="CNORode"));
model=readSIF("/Users/celine/modelMacNamara2012/ToyModelPB_TRUE_plus.sif");
#model=readSIF("/Users/celine/model2nodes/model2nodes.sif");
#model=readSIF(system.file("doc", "ToyModelMMB_FeedbackAnd.sif",
#         package="CNORode"));
#model=readSIF("/Users/celine/modelMacNamara2012/ToyModelPB_TRUE_plus.sif");

# Having loaded the data set and corresponding model, we run a check
# to make sure that our data and model were correctly loaded and that our
# data match our model (i.e. that species that were inhibited/stimulated
# measured in our data set are present in our model).
# Then 
#checkSignals(CNOlistToy,pknmodel) # also checked in gaBinaryT1() function

# View of the PKN model previously loaded from a .sif file
# here, the model is not trained on the data
#plotModel(pknmodel, CNOlistToy)


###############################
##  Preprocessing the model  ##
###############################

# 3 steps: removal of non-observable/non-controllable species, compression,
# and expansion
# preprocessing function is available and makes the following 3 steps in 1
# command line
model <- preprocessing(CNOlistToy, model, expansion=FALSE, compression=TRUE, cutNONC=TRUE, verbose=FALSE)
#checkSignals(CNOlistToy,model) # also checked in gaBinaryT1() function
#plotModel(model, CNOlistToy)


#############################
##  Training of the model  ##
#############################

# creation of a vector as long as the number of influences in the model
# obtained after preprocessing
initBstring<-rep(1,length(model$reacID))
#print(initBstring)
# This function is the genetic algorithm to be used to optimise a model by
# fitting to data containing one time
#timeHist = c()
#scoreHist = c()
#for (x in seq(1:15)) {
  startRun <- proc.time()
  ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy, model=model, initBstring=initBstring, popSize=10, maxGens=100,elitism=5,#sizeFac=0,
                       verbose=TRUE, scoreT0=TRUE, initState=TRUE)#, nameSim=nameSim)
  timeExec <- proc.time()-startRun
#  timeHist <- append(timeHist, timeExec[1])
#  scoreHist <- append(scoreHist,ToyT1opt$bScore)
#}
#pdf("evolTime.pdf")
#hist(timeHist)
#dev.off()
#pdf("evolScore.pdf")
#hist(scoreHist)
#hist(scoreHist, breaks = 10)
#hist(scoreHist, breaks = 12)
#dev.off()
print(timeExec)
#ToyT1opt
ToyT1opt
# get an eye of the function Help (section Value) to better underestand the
# returned values

bs <- ToyT1opt$bString
#score <- computeScoreT1(CNOlistToy, model, bs, simList=NULL, indexList=NULL, sizeFac=0.0001, NAFac=1, timeIndex=2)


# This function takes a model and an optimised bitstring, it cuts the
# model according to the bitstring and plots the results of the simulation
# along with the experimental data
#cutAndPlot(model=model, bStrings=list(ToyT1opt$bString), CNOlist=CNOlistToy,plotPDF=TRUE)
modelCut <- cutModel(model, bs)
sim_data = cSimulator(CNOlistToy, modelCut, scoreT0=TRUE, initState=TRUE)
plotParams=list(margin=0.1, width=15, height=12, cmap_scale=1, cex=1.6, ymin=NULL)
plotOptimResultsPan(sim_data, yInterpol=NULL, xCoords=NULL,
             CNOlist=CNOlist(CNOlistToy), formalism="ode", pdf=FALSE,
             plotParams=plotParams,pdfFileName="", tPt=NULL)

# one row for a treatment, each species has its column

# plots the evolution of best fit and average fit against generations on a .pdf file
pdf("evolFitToyT1.pdf")
plotFit(optRes=ToyT1opt)
dev.off()


####################################
##  Plotting the optimised model  ##
####################################

#par(mfrow=c(1,2))
#pdf("bestTopology_PKN.pdf")
#plotModel(model, CNOlistToy, bString=ToyT1opt$bString)
#bs = mapBack(model, pknmodel, ToyT1opt$bString)
#plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
#dev.off()
#par(mfrow=c(1,1))

#par(mfrow=c(1,2))
pdf("bestTopology_PKN.pdf")
plotModel(model, CNOlistToy, bString=ToyT1opt$bString)
#bs = mapBack(model, pknmodel, ToyT1opt$bString)
#plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
dev.off()
#par(mfrow=c(1,1))

# Processed model (left) and original PKN (right)
# The edges are on (black or red) or off (grey or pink) according to the
# best set of parameters found during the optimisation (the best bit string)


############################
##  Writing your results  ##
############################

#writeScaffold(modelComprExpanded=model, optimResT1=ToyT1opt, optimResT2=NA, modelOriginal=pknmodel, CNOlist=CNOlistToy)
#writeNetwork(modelOriginal=pknmodel, modelComprExpanded=model, optimResT1=ToyT1opt, optimResT2=NA, CNOlist=CNOlistToy)
#namesFilesToy<-list(dataPlot="ToyModelGraph.pdf", evolFitT1="evolFitToyT1.pdf", evolFitT2=NA, simResultsT1="SimResultsT1_1.pdf", simResultsT2=NA, scaffold="Scaffold.sif", scaffoldDot="Scaffold.dot", tscaffold="TimesScaffold.EA", wscaffold="weightsScaffold.EA", PKN="PKN.sif", PKNdot="PKN.dot", wPKN="TimesPKN.EA", nPKN="nodesPKN.NA")
#writeReport(modelOriginal=pknmodel, modelOpt=model, optimResT1=ToyT1opt, optimResT2=NA, CNOlist=CNOlistToy, directory="testToy", namesFiles=namesFilesToy, namesData=list(CNOlist="Toy",model="ToyModel"))
