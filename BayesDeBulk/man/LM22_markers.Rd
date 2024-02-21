\name{LM22_markers}
\alias{LM22_markers}
\title{Compute cell-type markers from LM22 signature matrix.}

\usage{
LM22_markers(Y)}
        
\arguments{
\item{Y} {List object containing expression data for each omic, \code{Y=list(y_1,y_2)} where \code{y_j} is a \code{(p_j x n_j)} matrix with rows corresponding to genes and columns to samples. Genes and samples can vary across datasets. Rows of \code{Y[[j]]} does not need to be ordered in a specific way. Rownames need to be specified as gene names.}

}

\value{
Two dimensional matrix contaning the list of markers differentiating each pair of cell-types. The first column contains the name of cell-type A while the second column of cell-type B. The third column includes the list of markers overexpressed in cell-type A compared to cell-type B. This matrix is the input taken by \code{BayesDeBulk} function.

}

\description{
This function returns list of markers differentiating each pair of cells contained in the LM22 signature matrix from Newman et al (2015).

}
\examples{
 
 # --- Generate data sets
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-c("ABCB4","ABCB9","ACAP1","ACP5","ACHE")  # genes name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 rownames(data1)<-rownames(data2)<-genes.name
 
 colnames(data1)<-colnames(data2)<-paste("Sample",seq(1,n1))

 markers<-LM22_markers(list(data1,data2)) 

}
\references{
Petralia, F., Krek, A., Calinawan, A.P., Feng, S., Gosline, S., Pugliese, P., Ceccarelli, M. and Wang, P., 2021. BayesDeBulk: a flexible Bayesian algorithm for the deconvolution of bulk tumor data. bioRxiv, pp.2021-06.

Newman, A.M., Liu, C.L., Green, M.R., Gentles, A.J., Feng, W., Xu, Y., Hoang, C.D., Diehn, M. and Alizadeh, A.A., 2015. Robust enumeration of cell subsets from tissue expression profiles. Nature methods, 12(5), pp.453-457.
}
