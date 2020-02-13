CellNOptR-MaBoSS
=========
[![Build Status](https://travis-ci.org/saezlab/CellNOptR.svg?branch=master)](https://travis-ci.org/saezlab/CellNOptR)

Training of boolean logic models of signalling networks using prior knowledge networks and perturbation data with a stochastic simulator.

- Please visit [CellNOptR](https://saezlab.github.io/CellNOptR/) for details about the project (references, news, ...)


## Installation:

Here CellNOptR-MaBoSS has been implemented within the CellNOptR package. Before starting, make sure you have installed the latest version of R. For more information and download of R, please refer to [R project page](http://www.r-project.org/) . For more information about how to  install R packages, please refer to [Installing packages](http://cran.r-project.org/doc/manuals/R-admin.html#Installing-packages).
These packages rely on several Bioconductor package (e.g., RBGL, graph, methods, etc.). As an example, you can
install RGBL package by typing:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RBGL")
```

### Installation from GitHub
Using the `devtools` package you can install the latest version from the GitHub repository:
```
if(!require("devtools")) install.packages("devtools")   # installs devtools package if not already installed
devtools::install_github("saezlab/CellNOptR-MaBoSS")
```

Finally we can call the package:
```
library(CellNOptR)
```

### Example
An example about how to run CellNOptR-MaBoSS is provided [here](https://github.com/saezlab/CellNOptR-MaBoSS/tree/toy_example).

## Feedbacks, bug reports, features
Feedbacks and bugreports are always very welcomed!  
Please use the Github issue page to report bugs or for questions: https://github.com/saezlab/CellNOptR/issues.
Many thanks!
