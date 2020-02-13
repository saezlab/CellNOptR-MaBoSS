# CellNOptR-MaBoSS Example
Here we show an application of CellNOptR-MaBoSS over the ToyMMB case study. More details about this case-study can be found [here](https://saezlab.github.io/CellNOptR/5_Models%20and%20Documentation/).

## Loading the package and the data

After installation of [CellNoptR-MaBoSS](https://github.com/saezlab/CellNOptR-MaBoSS/blob/package_development/Readme.md) package, users can load the `CellNOptR` package and the toy case study after reading the network from the *SIF* and the data from the *MIDAS* file:

```R
library(CellNOptR)

CNOlistToy=CNOlist("ToyDataMMB.csv");
model=readSIF("ToyPKNMMB.sif");
```

## Preprocessing and setting useful parameters

Our network can be compressed and we can also set other parameters useful for the CellNOptR analysis as described below:

```R
model <- preprocessing(CNOlistToy, model, expansion=FALSE, compression=TRUE, 
                       cutNONC=TRUE, verbose=FALSE) ##compressing the network
                       
multiTP <- FALSE  ##case-study without multiple time-points, otherwise it should be set multiTP <- TRUE
initBstring<-rep(1,length(model$reacID)) ##setting the initial bit-string to evaluate
mbsPath = "~/Downloads/MaBoSS-env-2.0/" ##setting the download path of MaBoSS
```

## Running CellNOptR-MaBoSS

Now we perfor the CellNoptR-MaBoSS analysis through the *computeMaBoSS* function:

```R
ToyT1opt<-computeMaBoSS(CNOlist=CNOlistToy, model=model, mbsPath, initBstring=initBstring, popSize=10, maxGens=10,elitism=5,
                           verbose=FALSE, scoreT0=TRUE, initState=TRUE, multiTP=multiTP, ttime=startRun)
```

## Obtaining results

Finally we show the results and plot the optimal model

```R
print(ToyT1opt)

pdf("bestTopology_PKN.pdf")
plotModel(model, CNOlistToy, bString=ToyT1opt$bString, graphvizParams = list(fontsize=46, nodeWidth=3))
dev.off()
```
