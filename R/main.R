# Rscript main.R --multiomic=${multiomic} --abundanceFile=${abundanceFile} --expressionFile=${expressionFile} --signatureMatrix=${signatureMatrix}

library("R.utils")
library(matrixStats)
library(base)

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

script_dir <- file.path(dirname(thisFile()))
source(paste(script_dir, '/function_bulk_repulsive_github.R', sep=''))

args <- commandArgs(asValue=TRUE, excludeReserved=TRUE)
keys <- attachLocally(args)

if (!('multiomic' %in% keys)) {
  multiomic = FALSE
  print('using default singleomic')
} else {
  multiomic <- eval(parse(text=multiomic))
}

if (multiomic == FALSE) {
  if ('expressionFile' %in% keys) {
    if ('abundanceFile' %in% keys) {
      warning('Two input files detected, but multiomic set to FALSE. Running single-omic with expression file.')
    }
    df_f1 <- expressionFile # path to expression file
  } else {
    df_f1 <- abundanceFile
  }
}

signature_matrix_f <- signatureMatrix # path to signature matrix

if (multiomic) {
  df_f1 <- expressionFile
  df_f2 <- abundanceFile
}

# read input files into dataframes 
df1 <- read.csv(df_f1, sep = ",", row.names = 1) # -- as input we can specify multiple files for multi-omic learning
if (multiomic==TRUE) df2 <- read.csv(df_f2, sep = ",", row.names = 1) # -- if multiomic learning, load additional data file

markers <- read.csv(signature_matrix_f, sep = ",",header=TRUE, row.names = 1) # -- genes x cell.type matrix, with (i,j) element equal to 1 if gene i is a marker of cell-type j and 0 otherwise

# -- transform ref in a list object
ref<-list()
for (s in 1:dim(markers)[2]) ref[[s]]<-rownames(markers)[markers[,s]==1]
cellTypeNames<-names(ref)<-colnames(markers)

# check that dataframe has no NA
if (anyNA(df1)) {
  stop("Missing values found in data 1. Missing values not allowed.")
}

if (multiomic) {
  if (anyNA(df2)) {
    stop("Missing values found in data 2. Missing values not allowed.")
  }
}

# slim down signature matrix and df based on intersecting genes
if (multiomic) {
  for (s in 1:length(ref)) ref[[s]] <- intersect(ref[[s]], unique(c(rownames(df1),rownames(df2))))
  df1 <- as.data.frame(df1[intersect(rownames(df1),unique(unlist(ref))),])
  df2 <- as.data.frame(df2[intersect(rownames(df2),unique(unlist(ref))),])
  df2<-df2[,match(colnames(df1),colnames(df2))] # match samples in df2 to that of df1
  
} else {
  for (s in 1:length(ref)) ref[[s]] <- ref[[s]][intersect(ref[[s]], rownames(df1))]
  df1 <- as.data.frame(df1[intersect(rownames(df1),unique(unlist(ref))),])
}  

data.1<-df1
if (multiomic) data.2<-df2

# generate index.matrix, k.fix, and prior from list.gene
k.fix <- length(ref)

compute.index.matrix<-function(list.gene,data){
  k.fix<-length(list.gene)
  index.matrix<-c(NULL,NULL,NULL)
  for (s in 1:k.fix){
  for (k in 1:k.fix) {
    if (s!=k){
    mg<-match(list.gene[[s]],list.gene[[k]])
    gene<-list.gene[[s]][is.na(mg)]
    if (length(gene)>0) index.matrix<-rbind(index.matrix,cbind(rep(s,length(gene)),rep(k,length(gene)),gene))
    
    }
    }
  }

  # set gene symbols to indexes in df
  mg<-match(index.matrix[,3],rownames(data))
  index.matrix[,3]<-mg
  index.matrix<-apply(index.matrix,2,as.numeric)
  
return(index.matrix)
}

index.matrix.1<-compute.index.matrix(ref,data.1)
if (multiomic) index.matrix.2<-compute.index.matrix(ref,data.2)

# -- combine data
if (multiomic) {
  data<-rbind(data.1,data.2) 
  index.matrix.2[,3]<-index.matrix.2[,3]+dim(data.1)[1]
  index.matrix<-rbind(index.matrix.1,index.matrix.2) 
  
  rownames(data)<-c(paste(rownames(data.1),"-RNA"),paste(rownames(data.2),"-Pro"))
} else {
  index.matrix<-index.matrix.1
  data<-data.1
}

prior <- matrix(0,dim(data)[1],k.fix)

# set n.iter and burn.in
n.iter <- 10000 
burn.in <- 1000

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

# -- derive point estimate for each cell type (average across MCMC iterations)
weights<-matrix(0,dim(gibbs[[1]][[1]]),length(gibbs[[1]]))
mu<-matrix(0,dim(gibbs[[2]][[1]]),length(gibbs[[1]]))
for (s in 1:length(gibbs[[1]])) {
  weights[,s]<-apply(gibbs[[1]][[s]],1,mean)
  mu[,s]<-apply(gibbs[[2]][[s]],1,mean)
}
weights<-apply(weights,1,function(x) (x/sum(x))) # -- normalize weights to sum to one
colnames(weights)<-colnames(mu)<-colnames(ref)

write.table(weights,"weights_bayesdebulk.tsv",row.names=TRUE, col.names=TRUE, sep='\t', quote = FALSE) # - matrix (samples x cell type), with (i,j) element being the fraction of the jth  cell type for sample i 
write.table(mu,"mean_parameter_bayesdebulk.tsv",row.names=TRUE, col.names=TRUE, sep='\t', quote = FALSE) # - matrix (genes x cell type), with (i,j) element being the mean of the i-the marker for the j-th cell type

# save gibbs object and signature matrix to rda file
save(gibbs,nu,weights, list.gene.1, index.matrix, file="output_bayes_de_bulk.Rdata")
