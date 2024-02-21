\name{BayesDeBulk}
\alias{BayesDeBulk}
\title{BayesDeBulk to perform the deconvolution of bulk proteogenomics data.}

\usage{
BayesDeBulk(n.iter,burn.in,Y,markers)}
        
\arguments{
\item{n.iter} {Numeric value: number of Gibbs sampler iterations}

\item{burn.in} {Numeric value: number of iterations to be discarded (burn-in)}

\item{Y} {List object containing expression data for each omic, \code{Y=list(y_1,y_2)} where \code{y_j} is a \code{(p_j x n_j)} matrix with rows corresponding to genes and columns to samples. Genes and samples can vary across datasets. Missing values are not allowed. Rows of \code{Y[[j]]} does not need to be ordered in a specific way.}

\item{markers} {Two dimensional matrix contaning the list of markers differentiating each pair of cell-types. The first column contains the name of cell-type A while the second column of cell-type B. The third column includes the list of markers overexpressed in cell-type A compared to cell-type B. Funtion \code{LM22_markers} produces this matrix based on LM22 signature matrix.}

}

\value{
A list object containing: 1) \code{cell.fraction}: the \code{n x K} matrix of estimated cell-type fractions. Each element \code{(i,j)} contains the cell-type fraction of the \code{j} cell-type and the \code{i} patient. (2) \code{cell.expression}: the expression of genes/proteins in different cell types.
}

\description{
This function computes deconvolution analysis via BayesDeBulk.

}
\examples{

##---------------------------------------------------------------------------##
#-- example 1 (data-sets with same sample size)
##---------------------------------------------------------------------------##

# --- Generate data sets
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 rownames(data1)<-rownames(data2)<-genes.name
 
 colnames(data1)<-colnames(data2)<-paste("Sample",seq(1,n1))
 
 markers<-rbind(c("Cell 1","Cell 2","G1"),c("Cell 1","Cell 2","G2"),c("Cell 2","Cell 1","G3"),c("Cell 2","Cell 1","G4"))

 ## ---  Run BayesDeBulk
  out<-BayesDeBulk(n.iter=20,burn.in=10,Y=list(data1,data2),markers=markers)


##---------------------------------------------------------------------------##
#-- example 2 (data-sets with different sample sizes)
##---------------------------------------------------------------------------##

 # --- Generate data sets
 n1<-20               # sample size for each data sets
 n2<-10
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n2)       # generate data2
 rownames(data1)<-rownames(data2)<-genes.name
 
 colnames(data1)<-paste("Sample",seq(1,n1))
 colnames(data2)<-paste("Sample",seq(1,n2))

 markers<-rbind(c("Cell 1","Cell 2","G1"),c("Cell 1","Cell 2","G2"),c("Cell 2","Cell 1","G3"),c("Cell 2","Cell 1","G4"))
 
 ##---------------------------------------------------------------------------##
 ## ---  Run BayesDeBulk
 
  out<-BayesDeBulk(n.iter=20,burn.in=10,Y=list(data1,data2),markers=markers)

}
\references{
Petralia, F., Krek, A., Calinawan, A.P., Feng, S., Gosline, S., Pugliese, P., Ceccarelli, M. and Wang, P., 2021. BayesDeBulk: a flexible Bayesian algorithm for the deconvolution of bulk tumor data. bioRxiv, pp.2021-06.

}
