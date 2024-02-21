
gibbs.sampling.missing.samples<-function(n.iter,Y,p,n,K,index.matrix,burn.in,mean.prior,sigma.prior){

  # -- prior
  alpha.sigma<-3
  beta.sigma<-1 # -- prior inverse-gamma for variance parameter
  alpha<-rep(1,K)
  tau<-2; eta<-2
  
  mean.save<-NULL
  TP<-FP<-NULL
  
  # -- initial parameter
  mu<-sigma<-NULL
  mean.initial<-matrix(0,p,K)
  rho<-.5
  alpha<-rep(1,K)
  
  mean.initial[,1]<-rnorm(p,mean=mean.prior[,1],sd=sigma.prior); mean.initial[index.matrix[index.matrix[,1]==1,3],1]<-mean.initial[index.matrix[index.matrix[,1]==1,3],1]+1;
  mu[[1]]<-mean.initial[,1];  sigma[[1]]<-rinvgamma(p,1,1);
  for (k in 2:K) { # -- sample mu and sigma from prior
    
    if (sum((index.matrix[,2]==k & index.matrix[,1]<k) | (index.matrix[,1]==k &index.matrix[,2]< k))>0){
    # -- sample mean for unconstrained genes from normal distribution
    gene.constrain<-unique(index.matrix[(index.matrix[,2]==k & index.matrix[,1]<k) | (index.matrix[,1]==k &index.matrix[,2]< k) ,3])
    mean.initial[-gene.constrain,k]<-rnorm(p-length(gene.constrain),mean=mean.prior[-gene.constrain,k],sd=1)
    
    # -- sample mean for constrained genes upregulated in k-th component from truncated normal
    if (sum(index.matrix[,1]==k &index.matrix[,2]< k)!=0){
    gene.higher<-unique(index.matrix[index.matrix[,1]==k & index.matrix[,2]< k,3])
    lower.bound<-rep(NA,p)
    for (s in 1:(k-1)) {
       gene.higher.s<-index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3] 
       lower.bound[gene.higher.s]<-apply(cbind(lower.bound[gene.higher.s],mean.initial[gene.higher.s,s]),1,function(x) max(x[!is.na(x)]))
    }
    lower.bound<-lower.bound[gene.higher]
    bound<-(-tau/log(rho))^(1/eta); 
    bound<-bound+lower.bound 
    mean.initial[gene.higher,k]<-rtruncnorm(1, a=bound, b=Inf, mean =mean.prior[gene.higher,k], sd = sigma.prior)
    }
    
    if (sum((index.matrix[,2]==k & index.matrix[,1]<k))!=0) {
    gene.lower<-unique(index.matrix[index.matrix[,2]==k & index.matrix[,1]<k,3]) 
    upper.bound<-rep(NA,p)
    for (s in 1:(k-1)) { 
      gene.lower.s<-index.matrix[index.matrix[,1]==s & index.matrix[,2]==k,3] 
      upper.bound[gene.lower.s]<-apply(cbind(upper.bound[gene.lower.s],mean.initial[gene.lower.s,s]),1,function(x) min(x[!is.na(x)]))
    }
    upper.bound<-upper.bound[gene.lower]
    
    bound<-(-tau/log(rho))^(1/eta); 
    bound<- upper.bound-bound 
    mean.initial[gene.lower,k]<-rtruncnorm(1, a=-Inf, b=bound, mean =mean.prior[gene.lower,k], sd = sigma.prior)
    }
    } else {
      mean.initial[,k]<-rnorm(p,mean=0,sd=1)
    }
    mu[[k]]<-mean.initial[,k];
    }   # --- K dimensionali list, each element containing mean/sigma parameter of p genes
   sigma<-rinvgamma(p,1,1); 
  
  #######################################################################################################
  ## -- Run Gibbs Sampling
  pi.final<-NULL

   nu.mean<-mu.mean<-pi.mean<-vector("list", K)
  mean.final<-mean.initial
  
  nu<-rep(0.5,K) # -- probability of zero for each cell types
  gamma<-rinvgamma(K,3,1);  # -- variance of prior for pi's
  pi<-matrix(runif(n*k),n,K); 
  for (k in 1:K) pi[sample(n/2),k]<-0.0001
    
  for (j in 1:n.iter){ # -- for loop over iterations
    print(j)
    
    Z<-matrix(0,n,K)
    
    # --- sample pi's 
    for (j.n in 1:n){
      for (w in 1:(K)){

        # -- sample Zi for pi's spike and slab prior
        prob<-(nu[w]*dnorm(pi[j.n,w],0,sqrt(0.0001)))/(nu[w]*dnorm(pi[j.n,w],0,sqrt(0.0001))+(1-nu[w])*dnorm(pi[j.n,w],0,sqrt(gamma[w])))
        Z[j.n,w]<-rbinom(1,1,prob=prob)
        
        # -- smmple pi's
        index<-(!is.na(Y[,j.n])) # -- only genes with n missing
        T_w<-Y[index,j.n]
        for (jj in 1:K) {
          if (jj!=w) T_w<-T_w-pi[j.n,jj]*mu[[jj]][index]

        }
          
          if (Z[j.n,w]==0) sigma.V<-(gamma[w])
          if (Z[j.n,w]==1) sigma.V<-((0.0001))
          var.p<-(sum(mu[[w]][index]^2/sigma[index])+1/sigma.V)^(-1) 
          mean.p<-sum(T_w*mu[[w]][index]/sigma[index])*var.p
          
          # -- sample V from truncated normal distribution
          pi[j.n,w]<-rtruncnorm(1,a=0,b=1,mean=mean.p,sd=sqrt(var.p))

        }
    }
    
    for (w in 1:K){
      # -- update nu (hyperparameter of pi's prior)
      nu[w]<-rbeta(1,1+sum(Z[,w]==1),1+sum(Z[,w]==0))
    
    # -- update gamma (hyperparameter of pi's prior)
       gamma[w]<-1/rgamma(1,alpha.sigma+sum(Z[,w]==0)/2,beta.sigma+0.5*sum((pi[Z[,w]==0,w])^2))
    }
    # --- sample mean parameter
    
    for (k in 1:K) {
      
      total.mu<-0
      for (s in 1:K) { if (s!=k) total.mu<-total.mu+as.matrix(pi[,s]) %*% mu[[s]]}
      
      pi.mu<-NULL
      for (kkk in 1:length(mu[[k]])) pi.mu<-c(pi.mu,sum(pi[!is.na(Y[kkk,]),k]^2))
      
      var.k<-(1/sigma.prior+pi.mu/sigma)^(-1) 
      Y.residual<-as.matrix(Y - t(total.mu))
      Y.residual[is.na(Y.residual)]<-0 
      mu.k<-rowSums(sweep(Y.residual,2,pi[,k],"*"))/sigma 
      mean.p<-(var.k)*(mu.k+mean.prior[,k]/sigma.prior)
      var.p<-(var.k) 
      
      # -- sample mean for unconstrained genes from normal distribution
      if (sum(index.matrix[,2]==k | index.matrix[,1]==k)==0){
        mean.final[,k]<-rnorm(p,mean=0,sd=1)*sqrt((var.p))+mean.p
      }
      if (sum(index.matrix[,2]==k | index.matrix[,1]==k)!=0){
      gene.constrain<-unique(index.matrix[index.matrix[,2]==k | index.matrix[,1]==k ,3])
      mean.final[-gene.constrain,k]<-rnorm(p-length(gene.constrain),mean=0,sd=1)*sqrt((var.p[-gene.constrain]))+mean.p[-gene.constrain]

      # -- sample mean for constrained genes upregulated in k-th component from truncated normal
      gene.higher<-unique(index.matrix[index.matrix[,1]==k,3])
      mean.gene.h<-mean.p[gene.higher]; var.gene.h<-var.p[gene.higher]
      lower.bound<-rep(NA,p)
      for (s in 1:K) {
        if (s!=k & sum(index.matrix[,1]==k & index.matrix[,2]==s)>0) {
          gene.higher.s<-index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3]
          lower.bound[gene.higher.s]<-apply(cbind(lower.bound[gene.higher.s],mean.final[gene.higher.s,s]),1,function(x) max(x[!is.na(x)]))
        }
      }
      lower.bound<-lower.bound[gene.higher]
      bound<-(-tau/log(rho))^(1/eta); 
      if (bound<0) bound<-0
      bound<-bound+lower.bound 
      mean.final[gene.higher,k]<-rtruncnorm(1, a=bound, b=Inf, mean = mean.gene.h, sd = sqrt(var.gene.h))
      
      gene.lower<-unique(index.matrix[index.matrix[,2]==k,3]) 
      mean.gene.l<-mean.p[gene.lower]; var.gene.l<-var.p[gene.lower]
      upper.bound<-rep(NA,p)
      for (s in 1:K) { 
        if (s!=k & sum(index.matrix[,2]==k & index.matrix[,1]==s)>0) {
          gene.lower.s<-unique(index.matrix[index.matrix[,1]==s & index.matrix[,2]==k,3]) # -- select genes supposed to be up in s-th component compared to the k-th component
          upper.bound[gene.lower.s]<-apply(cbind(upper.bound[gene.lower.s],mean.final[gene.lower.s,s]),1,function(x) min(x[!is.na(x)]))
        }
      }
      upper.bound<-upper.bound[gene.lower]
      
      bound<-(-tau/log(rho))^(1/eta); 
      bound<- upper.bound-bound 
      mean.final[gene.lower,k]<-rtruncnorm(1, a=-Inf, b=bound, mean = mean.gene.l, sd = sqrt(var.gene.l))
      
      if (is.nan(sum(mean.final))) break
      
    }
      mu[[k]]<-mean.final[,k]
    }

  # -- Sample sigma from Full Conditional
   total.mu<-0
   for (k in 1:K) total.mu<-total.mu+as.matrix(pi[,k])%*%mu[[k]]
   total.mu<-t(total.mu)
   
   Y.residual<-(Y- total.mu); Y.residual[is.na(Y.residual)]<-0 
   beta.p<-(beta.sigma+0.5*rowSums((Y.residual)^2))
  
   sigma.final<-rep(0,p)      
   for (g in 1:p)  sigma.final[g]<- 1/rgamma(1,shape=alpha.sigma+sum(!is.na(Y.residual[g,]))/2,rate=beta.p[g])      
   sigma<-sigma.final
   if (is.na(sum(sigma))) break
   
# -- sample rho
      n.S<-NULL
      bound<-NULL
      for (k in 1:K)  {
        for (s in 1:K){
        if (s!=k) {
          gene<-unique(index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3])
          bound<-c(bound,min(exp(-tau*((mean.final[gene,k]-mean.final[gene,s]))^(-eta))))
        }
        }
      }  
      
        rho<-runif(1,min=0,max=min(bound))
        if (is.nan(rho)) break

        if (j>burn.in) { 
          for(k in 1:K) {pi.mean[[k]]<-cbind(pi.mean[[k]],pi[,k])
          mu.mean[[k]]<-cbind(mu.mean[[k]],mu[[k]])
          nu.mean[[k]]<-cbind(nu.mean[[k]],nu[k])
          }
        }
  }
  return(list(pi.mean,mu.mean,nu.mean))
}

gibbs.sampling<-function(n.iter,Y,p,n,K,index.matrix,burn.in,mean.prior,sigma.prior){
  
  # -- prior
  alpha.sigma<-3
  beta.sigma<-1 # -- prior inverse-gamma for variance parameter
  alpha<-rep(1,K)
  tau<-2; eta<-2
  
  mean.save<-NULL
  TP<-FP<-NULL
  
  # -- initial parameter
  mu<-sigma<-NULL
  mean.initial<-matrix(0,p,K)
  rho<-.5
  alpha<-rep(1,K)
  
  mean.initial[,1]<-rnorm(p,mean=mean.prior[,1],sd=sigma.prior); mean.initial[index.matrix[index.matrix[,1]==1,3],1]<-mean.initial[index.matrix[index.matrix[,1]==1,3],1]+1;
  mu[[1]]<-mean.initial[,1];  sigma[[1]]<-rinvgamma(p,1,1);
  for (k in 2:K) { # -- sample mu and sigma from prior
    
    if (sum((index.matrix[,2]==k & index.matrix[,1]<k) | (index.matrix[,1]==k &index.matrix[,2]< k))>0){
      # -- sample mean for unconstrained genes from normal distribution
      gene.constrain<-unique(index.matrix[(index.matrix[,2]==k & index.matrix[,1]<k) | (index.matrix[,1]==k &index.matrix[,2]< k) ,3])
      mean.initial[-gene.constrain,k]<-rnorm(p-length(gene.constrain),mean=mean.prior[-gene.constrain,k],sd=1)
      
      # -- sample mean for constrained genes upregulated in k-th component from truncated normal
      if (sum(index.matrix[,1]==k &index.matrix[,2]< k)!=0){
        gene.higher<-unique(index.matrix[index.matrix[,1]==k & index.matrix[,2]< k,3])
        lower.bound<-rep(NA,p)
        for (s in 1:(k-1)) {
          gene.higher.s<-index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3] 
          lower.bound[gene.higher.s]<-apply(cbind(lower.bound[gene.higher.s],mean.initial[gene.higher.s,s]),1,function(x) max(x[!is.na(x)]))
        }
        lower.bound<-lower.bound[gene.higher]
        bound<-(-tau/log(rho))^(1/eta); 
        bound<-bound+lower.bound 
        mean.initial[gene.higher,k]<-rtruncnorm(1, a=bound, b=Inf, mean =mean.prior[gene.higher,k], sd = sigma.prior)
      }
      
      if (sum((index.matrix[,2]==k & index.matrix[,1]<k))!=0) {
        gene.lower<-unique(index.matrix[index.matrix[,2]==k & index.matrix[,1]<k,3]) 
        upper.bound<-rep(NA,p)
        for (s in 1:(k-1)) { 
          gene.lower.s<-index.matrix[index.matrix[,1]==s & index.matrix[,2]==k,3] 
          upper.bound[gene.lower.s]<-apply(cbind(upper.bound[gene.lower.s],mean.initial[gene.lower.s,s]),1,function(x) min(x[!is.na(x)]))
        }
        upper.bound<-upper.bound[gene.lower]
        
        bound<-(-tau/log(rho))^(1/eta); 
        bound<- upper.bound-bound 
        mean.initial[gene.lower,k]<-rtruncnorm(1, a=-Inf, b=bound, mean =mean.prior[gene.lower,k], sd = sigma.prior)
      }
    } else {
      mean.initial[,k]<-rnorm(p,mean=0,sd=1)
    }
    mu[[k]]<-mean.initial[,k];
  }   # --- K dimensionali list, each element containing mean/sigma parameter of p genes
  sigma<-rinvgamma(p,1,1); 
  
  #######################################################################################################
  ## -- Run Gibbs sampling
  #######################################################################################################
  
  pi.final<-NULL

  nu.mean<-mu.mean<-pi.mean<-vector("list", K)
  mean.final<-mean.initial
  
  nu<-rep(0.5,K) # -- probability of zero for each cell types
  gamma<-rinvgamma(K,3,1);  # -- variance of prior for pi's
  pi<-matrix(runif(n*k),n,K); 
  for (k in 1:K) pi[sample(n/2),k]<-0.00001
  
  for (j in 1:n.iter){ # -- for loop over iterations
    print(j)
    
    Z<-matrix(0,n,K)
    
    # --- sample pi's 
    for (j.n in 1:n){
      for (w in 1:(K)){
        
        # -- sample Zi for pi's spike and slab prior
        prob<-(nu[w]*dnorm(pi[j.n,w],0,sqrt(0.00001)))/(nu[w]*dnorm(pi[j.n,w],0,sqrt(0.00001))+(1-nu[w])*dnorm(pi[j.n,w],0,sqrt(gamma[w])))
        Z[j.n,w]<-rbinom(1,1,prob=prob)
        
        # -- smmple pi's
        T_w<-Y[,j.n]
        for (jj in 1:K) {
          if (jj!=w) T_w<-T_w-pi[j.n,jj]*mu[[jj]]
          
        }
        
        if (Z[j.n,w]==0) sigma.V<-(gamma[w])
        if (Z[j.n,w]==1) sigma.V<-((0.00001))
        var.p<-(sum(mu[[w]]^2/sigma)+1/sigma.V)^(-1) 
        mean.p<-sum(T_w*mu[[w]]/sigma)*var.p
        
        # -- sample V from truncated normal distribution
        pi[j.n,w]<-rtruncnorm(1,a=0,b=1,mean=mean.p,sd=sqrt(var.p))
        
      }
    }
    
    for (w in 1:K){
      # -- update nu (hyperparameter of pi's prior)
      nu[w]<-rbeta(1,1+sum(Z[,w]==1),1+sum(Z[,w]==0))
      
      # -- update gamma (hyperparameter of pi's prior)
      gamma[w]<-1/rgamma(1,alpha.sigma+sum(Z[,w]==0)/2,beta.sigma+0.5*sum((pi[Z[,w]==0,w])^2))
    }

    # --- sample mean parameter
    
    for (k in 1:K) {
      
      total.mu<-0
      for (s in 1:K) { if (s!=k) total.mu<-total.mu+as.matrix(pi[,s]) %*% mu[[s]]}
      
      var.k<-(1/sigma.prior+sum(pi[,k]^2)/sigma)^(-1)
      mu.k<-rowSums(sweep(as.matrix(Y - t(total.mu)),2,pi[,k],"*"))/sigma 
      mean.p<-(var.k)*(mu.k+mean.prior[,k]/sigma.prior)
      var.p<-(var.k) 
      
      # -- sample mean for unconstrained genes from normal distribution
      if (sum(index.matrix[,2]==k | index.matrix[,1]==k)==0){
        mean.final[,k]<-rnorm(p,mean=0,sd=1)*sqrt((var.p))+mean.p
      }
      if (sum(index.matrix[,2]==k | index.matrix[,1]==k)!=0){
        gene.constrain<-unique(index.matrix[index.matrix[,2]==k | index.matrix[,1]==k ,3])
        mean.final[-gene.constrain,k]<-rnorm(p-length(gene.constrain),mean=0,sd=1)*sqrt((var.p[-gene.constrain]))+mean.p[-gene.constrain]
        
        # -- sample mean for constrained genes upregulated in k-th component from truncated normal
        gene.higher<-unique(index.matrix[index.matrix[,1]==k,3])
        mean.gene.h<-mean.p[gene.higher]; var.gene.h<-var.p[gene.higher]
        lower.bound<-rep(NA,p)
        for (s in 1:K) {
          if (s!=k & sum(index.matrix[,1]==k & index.matrix[,2]==s)>0) {
            gene.higher.s<-index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3] 
            lower.bound[gene.higher.s]<-apply(cbind(lower.bound[gene.higher.s],mean.final[gene.higher.s,s]),1,function(x) max(x[!is.na(x)]))
          }
        }
        lower.bound<-lower.bound[gene.higher]
        bound<-(-tau/log(rho))^(1/eta); 
        if (bound<0) bound<-0
        bound<-bound+lower.bound 
        mean.final[gene.higher,k]<-rtruncnorm(1, a=bound, b=Inf, mean = mean.gene.h, sd = sqrt(var.gene.h))
        
        gene.lower<-unique(index.matrix[index.matrix[,2]==k,3]) 
        mean.gene.l<-mean.p[gene.lower]; var.gene.l<-var.p[gene.lower]
        upper.bound<-rep(NA,p)
        for (s in 1:K) { 
          if (s!=k & sum(index.matrix[,2]==k & index.matrix[,1]==s)>0) {
            gene.lower.s<-unique(index.matrix[index.matrix[,1]==s & index.matrix[,2]==k,3])
            upper.bound[gene.lower.s]<-apply(cbind(upper.bound[gene.lower.s],mean.final[gene.lower.s,s]),1,function(x) min(x[!is.na(x)]))
          }
        }
        upper.bound<-upper.bound[gene.lower]
        
        bound<-(-tau/log(rho))^(1/eta); 
        bound<- upper.bound-bound 
        mean.final[gene.lower,k]<-rtruncnorm(1, a=-Inf, b=bound, mean = mean.gene.l, sd = sqrt(var.gene.l))
        
        if (is.nan(sum(mean.final))) break
        
      }
      mu[[k]]<-mean.final[,k]
    }
    
    # -- sample sigma from full conditional
    total.mu<-0
    for (k in 1:K) total.mu<-total.mu+as.matrix(pi[,k])%*%mu[[k]]
    total.mu<-t(total.mu)
    
    alpha.p<-(alpha.sigma+n/2); 
    beta.p<-(beta.sigma+0.5*rowSums((Y- total.mu)^2))
    
    sigma.final<-rep(0,p)      
    for (g in 1:p)  sigma.final[g]<- 1/rgamma(1,shape=alpha.p,rate=beta.p[g])      
    sigma<-sigma.final
    if (is.na(sum(sigma))) break
    
    # -- sample rho
    n.S<-NULL
    bound<-NULL
    for (k in 1:K)  {
      for (s in 1:K){
        if (s!=k) {
          gene<-unique(index.matrix[index.matrix[,1]==k & index.matrix[,2]==s,3])
          bound<-c(bound,min(exp(-tau*(mean.final[gene,k]-mean.final[gene,s])^(-eta))))
        }
      }
    }  
    
    rho<-runif(1,min=0,max=min(bound))
    if (is.nan(rho)) break
    
    
     if (j>burn.in) { 
      for(k in 1:K) {pi.mean[[k]]<-cbind(pi.mean[[k]],pi[,k])
      mu.mean[[k]]<-cbind(mu.mean[[k]],mu[[k]])
      }
    }
  }
  return(list(pi.mean,mu.mean))
}

BayesDeBulk<-function(n.iter,burn.in,Y,markers){
  
  if (length(Y)==2) multiomic=TRUE else multiomic=FALSE
    
  df1 <- Y[[1]]
  if (multiomic == TRUE) df2 <- Y[[2]]
  
  ref<-markers
  cellTypeNames<-unique(c(markers[,1],markers[,2]))
  
  mg<-match(ref[,1],cellTypeNames)
  ref[,1]<-mg
  
  mg<-match(ref[,2],cellTypeNames)
  ref[,2]<-mg
  
  if (max(df1[!is.na(df1)])>500) {
    df1<-log(df1+1)
    print("Log normalizing data 1.")
  }

  if (max(df2[!is.na(df2)])>500) {
    df1<-log(df2+1)
    print("Log normalizing data 2.")
  }
 
  # slim down signature matrix and data based on markers with prior information
  if (multiomic) {
    ref <- ref[!is.na(match(ref[,3], unique(c(rownames(df1),rownames(df2))))),]
    df1 <- as.data.frame(df1[intersect(rownames(df1),unique((ref[,3]))),])
    df2 <- as.data.frame(df2[intersect(rownames(df2),unique((ref[,3]))),])
    
    mg<-match(colnames(df1),colnames(df2)) # -- extract samples shared between two data
    df2.new<-df2[,mg[!is.na(mg)]] 
    df1.new<-df1[,!is.na(mg)] 

    sampleid<-colnames(df1)[!is.na(mg)]
    
    missing=FALSE
    
    mg<-match(colnames(df1),colnames(df2))
    if (sum(is.na(mg))>0){
      df1.new<-cbind(df1.new,df1[,is.na(mg)]) 
      df2.new<-cbind(df2.new,matrix(NA,dim(df2.new)[1],sum(is.na(mg))))
      missing=TRUE
      sampleid<-c(sampleid,colnames(df1)[is.na(mg)])
}
    mg<-match(colnames(df2),colnames(df1))
    if (sum(is.na(mg))>0){
    df2.new<-cbind(df2.new,df2[,is.na(mg)]) 
    df1.new<-cbind(df1.new,matrix(NA,dim(df1.new)[1],sum(is.na(mg))))
    missing=TRUE
    sampleid<-c(sampleid,colnames(df2)[is.na(mg)])
    }
    
    df1<-df1.new
    df2<-df2.new
    
    colnames(df1)<-colnames(df2)<-sampleid

      } else {
        ref <- ref[!is.na(match(ref[,3], unique(rownames(df1)))),]
        df1 <- as.data.frame(df1[intersect(rownames(df1),unique((ref[,3]))),])
        sampleid<-colnames(df1)
  }  
  
  data.1<-df1
  if (multiomic) data.2<-df2
  
  k.fix <- length(cellTypeNames)
  
  compute.index.matrix<-function(index.matrix,data){

    mg<-match(index.matrix[,3],rownames(data))
    index.matrix[,3]<-mg
    index.matrix<-index.matrix[!is.na(index.matrix[,3]),]
    
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
    
    rownames(data)<-c(paste(rownames(data.1),"-data1"),paste(rownames(data.2),"-data2"))
  } else {
    index.matrix<-index.matrix.1
    data<-data.1
  }
  
  prior.levels <- matrix(0,dim(data)[1],k.fix)
  
  # - Z-score each marker
  data<-t(apply(data,1,function(x) (x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))
  
  # Run Gibbs Sampling
if (missing==FALSE)  
  gibbs<-gibbs.sampling(
    n.iter=n.iter,
    data,
    p=dim(data)[1],
    n=dim(data)[2],
    k.fix,
    index.matrix,
    burn.in,
    mean.prior=prior.levels,
    sigma.prior=1
  )
  
  if (missing==TRUE)  
    gibbs<-gibbs.sampling.missing.samples(
      n.iter=n.iter,
      data,
      p=dim(data)[1],
      n=dim(data)[2],
      k.fix,
      index.matrix,
      burn.in,
      mean.prior=prior.levels,
      sigma.prior=1
    )
  
  # -- derive point estimate for each cell type (average across MCMC iterations)
  weights<-matrix(0,dim(gibbs[[1]][[1]]),length(gibbs[[1]]))
  mu<-matrix(0,dim(gibbs[[2]][[1]]),length(gibbs[[2]]))
  for (s in 1:length(gibbs[[1]])) {
    weights[,s]<-apply(gibbs[[1]][[s]],1,mean)
    mu[,s]<-apply(gibbs[[2]][[s]],1,mean)
  }
  weights<-t(apply(weights,1,function(x) (x/sum(x)))) # -- normalize weights to sum to one
  colnames(weights)<-colnames(mu)<-cellTypeNames
  rownames(weights)<-sampleid
  
  out<-NULL
  out$cell.fraction<-weights

  out$cell.expression<-mu
  rownames(out$cell.expression)<-rownames(data)
  
  return(out)
}










