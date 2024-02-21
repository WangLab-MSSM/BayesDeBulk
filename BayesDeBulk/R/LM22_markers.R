
LM22_markers<-function(Y){
  
  data(LM22) # -- load LM22 signature matrix
  if (length(Y)>1) { genes<-unique(c(rownames(Y[[1]]),rownames(Y[[2]])))}  else { genes<-rownames(Y[[1]])}
  
  cell.type<-colnames(LM22)
  
  index.matrix<-NULL
  for (s in 1:dim(LM22)[2]){
    for (k in 1:dim(LM22)[2]){
      
      if (k!=s){
        i<-(LM22[,s]>1000 & (LM22[,s]>3*LM22[,k]) )
        marker.unique<-rownames(LM22)[i]
        marker.unique<-marker.unique[!is.na(match(marker.unique,genes))]
        
        if (length(marker.unique)>=1)  index.matrix<-rbind(index.matrix,cbind(rep(cell.type[s],each=length(marker.unique)),
                                                 rep(cell.type[k],each=length(marker.unique)),
                                                 marker.unique))
        
      }
    }
  }
  
  return(index.matrix)
  
}
