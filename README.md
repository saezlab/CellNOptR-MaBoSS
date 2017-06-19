# CellNOptR-MaBoSS
CellNOptR with MaBoSS simulation

## Tools installation

* To download MaBoSS, go to this [webpage](https://maboss.curie.fr) and follow the instructions given for the installation.
* Download the CellNOptR-MaBoSS folder on your local machine. The package is contained in the `pipelineCNOR-MBSS` directory.

## Make the pipeline running on the model

* Compress and install the package as a library :  
`tar zcvf pipelineCNOR-MBSS.tar.gz pipelineCNOR-MBSS/`  
`R CMD INSTALL --no-html pipelineCNOR-MBSS.tar.gz`
* Go the MaBoSS directory, and copy the `ToyDataMMB.csv`, `ToyPKNMMB.sif` et `basicFile.cfg` from the CellNOptR-MaBoSS directory.
* Run the pipeline with the command line :  
`nohup Rscript ~/CellNOptR-MaBoSS/runPipeline.R &`

## Change the model

* A `.csv` and a `.sif` files are required. The first file contains the experiemental data and the second one gives the regulatory edges
of the model. For more informations : *Terfve, Camille, et al. "CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms." BMC systems biology 6.1 (2012): 133.*
* Replace files' names in the `runPipeline.R` lines 41 and 49 respectively. If the files are not in the MaBoSS directory, enter
the absolute path of the files.

## The outputs

* A file called `nohup.out` had store  all the printings that can occure during the execution.
* For each generation, a file `resultsGen_1` (e.g. for the first generation) is created. It is the consequence of the line 89 in the `runPipeline.R`, where the parameter `verbose` is TRUE.
* `Rplots.pdf` and `SimResultsT1_1.pdf` represent the fitting scores between the best solution found and the experimental data
* `bestTopology_PKN.pdf` gives first the topology obtained for the best solution
* `evolFitToyT1.pdf` represents the average score through generation and the best score along the time

