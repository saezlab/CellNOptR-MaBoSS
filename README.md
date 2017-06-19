# CellNOptR-MaBoSS
CellNOptR with MaBoSS simulation

## Tools installation

* To download MaBoSS, go to this [webpage](https://maboss.curie.fr) and follow the instructions given for the installation.
* Download the pipeline and the CellNOptR-MaBoSS folder on your local machine. The package is contained in the `pipelineCNOR-MBSS`
directory.

## Make the pipeline running no the model

* Compress and install the CellNOptR-MaBoSS package :  
`tar zcvf CellNOptR-MaBoSS.tar.gz CellNOptR-MaBoSS/`  
`R CMD INSTALL --no-html CellNOptR-MaBoSS.tar.gz`
* Go the MaBoSS directory, and copy the `ToyDataMMB.csv`, `ToyPKNMMB.sif` et `basicFile.cfg` from the CellNOptR-MaBoSS directory.
* Run the pipeline typing :  
`nohup Rscript ~/CellNOptR-MaBoSS/runPipeline.R &`

## Change the model

* A `.csv` and a `.sif` files are required. The first file contains the experiemental data and the second one gives the regulatory edges
of the model. For more informations :
* Replace files' names in the `runPipeline.R` lines 41 and 49 respectively. If the files are not in the CellNOptR-MaBoSS directory, enter
the absolute path of the files.
