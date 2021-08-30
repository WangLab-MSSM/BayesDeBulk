# Rscript main_docker_multiomic.R --multiomic=${multiomic} --abundanceFile=${abundanceFile} --expressionFile=${expressionFile} --signatureMatrix=${signatureMatrix} --rowMeansImputation=TRUE

source("./function_bulk_repulsive_github.R")

library("R.utils")
library(matrixStats)

args <- commandArgs(asValue=TRUE, excludeReserved=TRUE)
keys <- attachLocally(args)
multiomic=TRUE

# accept input file paths as command line args
df_f1 <- expressionFile # path to expression file
signature_matrix_f <- signatureMatrix # path to signature matrix

if (multiomic) {
  df_f2 <- abundanceFile # path to protein abundance file
}

if ('rowMeansImputation' %in% keys) { # boolean whether to run row means imputation, default TRUE
  rowMeansImputation <- eval(parse(text=rowMeansImputation))
  if (!is.logical(rowMeansImputation)) {
    rowMeansImputation <- TRUE
  }
} else {
  rowMeansImputation <- TRUE
}

# read input files into dataframes 
df1 <- read.csv(df_f1, sep = "\t", row.names = 1) # -- as input we can specify multiple files for multi-omic learning
if (multiomic==TRUE) df2 <- read.csv(df_f2, sep = "\t", row.names = 1) # -- if multiomic learning, load additional data file
ref <- read.csv(signature_matrix_f, sep = "\t", row.names = 1)

# if using rowMeansImputation, drop rows > 50% NA and impute
if (anyNA(df1) & rowMeansImputation) {
  print('Performing rowMeansImputation.')
  df1 <- df1[which(rowMeans(!is.na(df1)) > 0.5), ]
  k <- which(is.na(df1), arr.ind=TRUE)
  df1[k] <- rowMeans(df1, na.rm=TRUE)[k[,1]]
}

if (multiomic & anyNA(df2) & rowMeansImputation) {
  print('Performing rowMeansImputation.')
  df2 <- df2[which(rowMeans(!is.na(df2)) > 0.5), ]
  k <- which(is.na(df2), arr.ind=TRUE)
  df2[k] <- rowMeans(df2, na.rm=TRUE)[k[,1]]
}

# check that dataframe has no NA
if (anyNA(df1)) {
  stop("Missing values found in data 1. Consider setting rowMeansImputation=TRUE.")
}

if (multiomic & anyNA(df2)) {
  stop("Missing values found in data 1. Consider setting rowMeansImputation=TRUE.")
}

# slim down signature matrix and df based on intersecting genes
ref1 <- as.data.frame(ref[intersect(rownames(ref), rownames(df1)),])
df1 <- as.data.frame(df1[intersect(rownames(ref1), rownames(df1)),])

if (multiomic){
  ref2 <- as.data.frame(ref[intersect(rownames(ref), rownames(df2)),]) 
  df2 <- as.data.frame(df2[intersect(rownames(ref2), rownames(df2)),])
  
  df2<-df2[,match(colnames(df1),colnames(df2))] # match samples in df2 to that of df1
}

# create list.gene from signature matrix for data 1
cellTypeNames <- colnames(ref1)
marker.1 <- list()
taggedNames <- c()
for (i in seq(dim(ref1)[2])) {
  s = cellTypeNames[i]
  tempMarkers <- rownames(ref1[ref1[,s] > quantile(ref1[,s], prob=0.75),])
  markers <- c()
  if (length(args) > 2) {
    if (tolower(args[3]) == "pairwise") {
      tempTb = 2 * ref1[,cellTypeNames[cellTypeNames != s]] - ref1[,s]
      markers <- tempMarkers[tempMarkers %in% rownames(ref1[rowSums(tempTb < 0) == (length(cellTypeNames) - 1), ])]
    } else {
      markers <- tempMarkers[tempMarkers %in% rownames(ref1)[ref1[,s] > 2 * rowMedians(as.matrix(ref1))]]
    }
  } else {
    markers <- tempMarkers[tempMarkers %in% rownames(ref1)[ref1[,s] > 2 * rowMedians(as.matrix(ref1))]]
  }
  if (length(markers) == 0) {
    deOverMedians <- (ref1[tempMarkers,s] - rowMedians(as.matrix(ref1[tempMarkers,])))/rowMedians(as.matrix(ref1[tempMarkers,]))
    markers <- tempMarkers[deOverMedians > quantile(deOverMedians, prob=0.9)]
  }
  marker.1[[i]] <- markers
  taggedNames <- c(taggedNames, paste0(colnames(ref1)[i], "%", cellTypeNames[i], toString(i), "%", cellTypeNames[i], toString(i), ".txt"))
}
names(marker.1) <- taggedNames

list.gene.1 <- marker.1
data.1 <- df1

# create list.gene from signature matrix for data 2
if (multiomic==TRUE){
  cellTypeNames <- colnames(ref2)
  marker.2 <- list()
  taggedNames <- c()
  for (i in seq(dim(ref2)[2])) {
    s = cellTypeNames[i]
    tempMarkers <- rownames(ref2[ref2[,s] > quantile(ref2[,s], prob=0.75),])
    markers <- c()
    if (length(args) > 2) {
      if (tolower(args[3]) == "pairwise") {
        tempTb = 2 * ref2[,cellTypeNames[cellTypeNames != s]] - ref2[,s]
        markers <- tempMarkers[tempMarkers %in% rownames(ref2[rowSums(tempTb < 0) == (length(cellTypeNames) - 1), ])]
      } else {
        markers <- tempMarkers[tempMarkers %in% rownames(ref2)[ref2[,s] > 2 * rowMedians(as.matrix(ref2))]]
      }
    } else {
      markers <- tempMarkers[tempMarkers %in% rownames(ref2)[ref2[,s] > 2 * rowMedians(as.matrix(ref2))]]
    }
    if (length(markers) == 0) {
      deOverMedians <- (ref2[tempMarkers,s] - rowMedians(as.matrix(ref2[tempMarkers,])))/rowMedians(as.matrix(ref2[tempMarkers,]))
      markers <- tempMarkers[deOverMedians > quantile(deOverMedians, prob=0.9)]
    }
    marker.2[[i]] <- markers
    taggedNames <- c(taggedNames, paste0(colnames(ref2)[i], "%", cellTypeNames[i], toString(i), "%", cellTypeNames[i], toString(i), ".txt"))
  }
  
  names(marker.2) <- taggedNames
  list.gene.2 <- marker.2
  data.2 <- df2
}


# generate index.matrix, k.fix, and prior from list.gene
k.fix <- length(list.gene.1)

compute.index.matrix<-function(list.gene,data){
  k.fix<-length(list.gene)
  index.matrix<-c(NULL,NULL,NULL)
  for (s in 2:k.fix){
  for (k in 1:(s-1)) {
    index.s<-c(NULL,NULL,NULL)
    mg<-match(list.gene[[s]],list.gene[[k]])
    gene<-list.gene[[s]][is.na(mg)]
    if (length(gene)>0) index.s<-rbind(index.s,cbind(rep(s,length(gene)),rep(k,length(gene)),gene))
    
    mg<-match(list.gene[[k]],list.gene[[s]])
    gene<-list.gene[[k]][is.na(mg)]
    if (length(gene)>0) index.s<-rbind(index.s,cbind(rep(k,length(gene)),rep(s,length(gene)),gene))
    
    index.matrix<-rbind(index.matrix,index.s)
  }
  }
  # set gene symbols to indexes in df
  mg<-match(index.matrix[,3],rownames(data))
  index.matrix[,3]<-mg
  index.matrix<-index.matrix[!is.na(mg),]
  index.matrix<-apply(index.matrix,2,as.numeric)
  
return(index.matrix)
}

index.matrix.1<-compute.index.matrix(list.gene.1,data.1)
if (multiomic) index.matrix.2<-compute.index.matrix(list.gene.2,data.2)

# -- combine data
if (multiomic) {
  data<-rbind(data.1,data.2) 
  index.matrix.2[,3]<-index.matrix.2[,3]+dim(data.1)[1]
  index.matrix<-rbind(index.matrix.1,index.matrix.2) 
} else {
  index.matrix<-index.matrix.1
  data<-data.1
}

prior <- matrix(0,dim(data)[1],k.fix)

# set n.iter and burn.in
n.iter <- 10
burn.in <- 1

# run gibbs
gibbs<-gibbs.sampling(
  n.iter=n.iter,
  data,
  p=dim(data)[1],
  n=dim(data)[2],
  k.fix,
  index.matrix,
  burn.in,
  mean.prior=prior,
  sigma.prior=1
)

write.table(gibbs[[1]][[1]],"output_bayes_de_bulk.tsv",row.names=TRUE, col.names=NA, sep='\t', quote = FALSE)

# save gibbs object and signature matrix to rda file
save(gibbs, list.gene.1, index.matrix, file="output_bayes_de_bulk.Rdata")
