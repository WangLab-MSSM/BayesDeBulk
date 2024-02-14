# BayesDeBulk

* [Introduction](#introduction)
* [Citation](#citation)
* [Installing BayesDeBulk package in R](#running-from-the-command-line)
* [Algorithm schematic](#algorithm-schematic)



## Introduction
Characterizing the tumor microenvironment is crucial in order to improve responsiveness to immunotherapy and develop new therapeutic strategies. The fraction of different cell-types in the tumor microenvironment can be estimated based on transcriptomic profiling of bulk RNA data via deconvolution algorithms. Some of these algorithms, defined as reference-based, rely on a reference signature containing gene expression for different cell-types. The limitation of these methods is that such signature is derived from gene expression of pure cells which might not be consistent to the transcriptomic profiling of different cells in solid tumors.sOn the other hand, reference-free methods, such as non-negative matrix factorization, require only a set of cell-specific markers to perform the deconvolution;showever, once estimated the labeling of different components might be problematic. To overcome these limitations, we propose BayesDeBulk - a new referencefree Bayesian method for bulk deconvolution. Given a list of markers expressedsin each cell-type (cell-specific markers), a repulsive prior is placed on the mean ofsgene expression in different cell-types to guarantee that cell-specific markers aresupregulated in a particular component. Contrary to existing reference-free methods, the labeling of different components is decided a priori through a repulsivesprior. On the other hand, the advantage over reference-based algorithms is thatsthe cell fractions as well as gene expression of different cells are estimated from thesdata, simultaneously. Given its flexibility, BayesDeBulk can be utilized to performs the deconvolution based on other data types such as proteomic profiling or thesintegration of both transcriptomic and proteomic profiling.

## Citation
For more information, please visit or cite the related preprint: 

[Petralia, F., Calinawan, A. P., Feng, S., Gosline, S., Pugliese, P., Ceccarelli, M., & Wang, P. (2021). BayesDeBulk: A Flexible Bayesian Algorithm for the Deconvolution of Bulk Tumor Data. doi.org/10.1101/2021.06.25.449763](https://www.biorxiv.org/content/10.1101/2021.06.25.449763v4)

## Installing BayesDeBulk package in R 
* Requires R >= 3.6
* Download BayesDeBulk_1.0.tar.gz file
* Install package from command line of R using

  install_packages("BayesDeBulk_1.0.tar.gz")


## Algorithm schematic
![alt text](./algorithm_schematic.png)
