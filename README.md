# BayesDeBulk

## Introduction
![alt text](./algorithm_schematic.png)

Characterizing the tumor microenvironment is crucial in order to improve responsiveness to immunotherapy and develop new therapeutic strategies. The fraction of different cell-types in the tumor microenvironment can be estimated based on transcriptomic profiling of bulk RNA data via deconvolution algorithms. Some of these algorithms, defined as reference-based, rely on a reference signature containing gene expression for different cell-types. The limitation of these methods is that such signature is derived from gene expression of pure cells which might not be consistent to the transcriptomic profiling of different cells in solid tumors.sOn the other hand, reference-free methods, such as non-negative matrix factorization, require only a set of cell-specific markers to perform the deconvolution;showever, once estimated the labeling of different components might be problematic. To overcome these limitations, we propose BayesDeBulk - a new referencefree Bayesian method for bulk deconvolution. Given a list of markers expressedsin each cell-type (cell-specific markers), a repulsive prior is placed on the mean ofsgene expression in different cell-types to guarantee that cell-specific markers aresupregulated in a particular component. Contrary to existing reference-free methods, the labeling of different components is decided a priori through a repulsivesprior. On the other hand, the advantage over reference-based algorithms is thatsthe cell fractions as well as gene expression of different cells are estimated from thesdata, simultaneously. Given its flexibility, BayesDeBulk can be utilized to performs the deconvolution based on other data types such as proteomic profiling or thesintegration of both transcriptomic and proteomic profiling.

For more information, please visit or cite the related preprint: [Petralia, F., Calinawan, A. P., Feng, S., Gosline, S., Pugliese, P., Ceccarelli, M., & Wang, P. (2021). BayesDeBulk: A Flexible Bayesian Algorithm for the Deconvolution of Bulk Tumor Data. Cold Spring Harbor Laboratory. doi.org/10.1101/2021.06.25.449763](https://www.biorxiv.org/content/10.1101/2021.06.25.449763v1)

## Running from the command line
* Requires R >= 3.6

The following command will perform tumor deconvolution with an input signature matrix for combined multi-omic data, including a protein abundance file and RNA expression file.

```sh
cd R
Rscript main.R --multiomic=TRUE --abundanceFile='../test_data/proteo_dummy.tsv' --expressionFile='../test_data/RNA_dummy.tsv' --signatureMatrix='../test_data/LM22_combined_cell_types.tsv' --rowMeansImputation=TRUE
```

## Docker image

The code and test data are available as a Docker image tagged [cptacdream/bayesdebulk](https://hub.docker.com/repository/docker/cptacdream/bayesdebulk).


## Inputs

| flag               | type     | default | description                                                                                             |
|--------------------|----------|---------|---------------------------------------------------------------------------------------------------------|
| multiomic          | boolean  | FALSE   | indicates whether to compute tumor deconvolution with both RNA expression and protein abundance         |
| abundanceFile      | filepath |         | path to tab-separated protein abundance table file, where rows are gene symbols and columns are samples |
| expressionFile     | filepath |         | path to tab-separated RNA expression table file, where rows are gene symbols and columns are samples    |
| signatureMatrix    | filepath |         | path to signature matrix table file, where rows are gene symbols and columns are cell types             |
| rowMeansImputation | boolean  | TRUE    | indicates whether to perform row means imputation for NA values in -omics files                         |


## Outputs
