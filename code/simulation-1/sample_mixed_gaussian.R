################################################################################
# Generate mixed type data from a mixture of multivariate Gaussian distributions
# code from https://github.com/Alsanchez13/Clustering-approaches-for-mixed-type-data-A-comparative-study/blob/master/Code
################################################################################

mixt_data<-function(cont_var=NULL,cont_cov=NULL, over_cont=NULL, ord_var=NULL,ord_lev=NULL, over_ord=NULL, bin_var=NULL,over_bin=NULL,nom_var=NULL, nom_lev=NULL, over_nom=NULL, clusters=NULL, obser=NULL,symm=TRUE)
{
  library(Matrix)
  library(mvtnorm)
  #data check
  
  if(is.null(clusters)){stop("PLease specify the total number of clusters.")}
  if(is.null(obser)){stop("PLease specify the total number of observations.")}
  
  if(!is.null(clusters)&&!is.null(obser)&&(is.null(ord_var)&&is.null(bin_var)&&is.null(cont_var)&&is.null(nom_var))){stop("Please specify the data type.")}
  # partition cluster sizes 
  
  if(symm==TRUE)
  {clus_dat_part<-floor(rep((obser/clusters),clusters))
  clus_dat_part[clusters]<-obser-sum(clus_dat_part[-clusters])
  }
  
  if(symm!=TRUE)
  {
    prob_part<-runif((clusters-1),.1,(1/(clusters)))
    prob_part[clusters]<-1-sum(prob_part)
    clus_dat_part=as.vector(rmultinom(1,obser, prob=prob_part))
  }
  
  tot_data<-vector() 
  #generate continuous data
  if(!(is.null(cont_var)))
  {
    if(is.null(over_cont))
    {
      over_cont<-rep(0,times=cont_var)
    }
    else{if(length(over_cont)!=cont_var){stop("Please make sure that the overlap parameters are of the same length as the continous variables")}}
    
    #generate covariance matrix
    
    cov<-diag(1,nrow=cont_var)
    
    data<-matrix(0,ncol=cont_var, nrow=1)
    mean_clus<-matrix(0,ncol=cont_var, nrow=clusters)
    mean_int<-seq(from=0, to=10, length.out = cont_var)
    cont_data<-matrix(0, nrow=1,ncol=cont_var)
    for(m in 1:clusters)
    {
      mean_clus[m,]<-mean_int
      temp<-mvtnorm::rmvnorm(clus_dat_part[m], mean=mean_int, sigma=cov)
      mean_int<-(mean_int+5)+(-5*over_cont)
      cont_data<-rbind(cont_data, temp)
    }
    
    cont_data<-cont_data[-1,]
    tot_data<-c(tot_data,as.vector(cont_data))
  }
  
  #generate ord _data
  if(!is.null(ord_var))
  {
    if(is.null(over_ord))
    {
      over_ord<-rep(0,times=ord_var)
    }
    else{if(length(over_ord)!=ord_var){stop("Please make sure that the overlap parameters are of the same length as the ordinal variables")}}
    
    mean_int_ord<-seq(from=0, to=10, length.out = ord_var)
    ord_data<-matrix(0, nrow=1,ncol=ord_var)
    cov_ord<-diag(1,nrow=ord_var)
    for(m in 1:clusters)
    {
      #mean_clus[m,]<-mean_int
      temp<-abs(mvtnorm::rmvnorm(clus_dat_part[m], mean=mean_int_ord, sigma=cov_ord))
      mean_int_ord<-(mean_int_ord+7)+(-7*over_ord)
      ord_data<-rbind(ord_data, temp)
    }
    
    ord_data<-ord_data[-1,]
    ord_data<-as.matrix(ord_data)
    mean_int_ord<-seq(from=0, to=10, length.out = ord_var)
    ordinal_data<-vector()
    for(k in 1:ord_var)
    {
      
      int_lin<-seq(from=mean_int_ord[k],to=(mean_int_ord[k]+7*(clusters-1)), length.out=clusters*2)
      
      ord_lab<-rep(0,nrow(ord_data))
      
      for(l in 1:(length(int_lin)-1))
      {
        for(i in 1:nrow(as.matrix(ord_data)))
        {
          
          if(ord_data[i,k]>= int_lin[l] && ord_data[i,k]<int_lin[l+1])
          {
            ord_lab[i]=l
          }
          if(ord_data[i,k]>=int_lin[length(int_lin)])
          {
            ord_lab[i]=length(int_lin)
          } 
        }
      }
      
      ordinal_data<-c(ordinal_data,ord_lab)
    }
    
    tot_data<-c(tot_data,ordinal_data)
  }
  
  #generate binary Data 
  if(!(is.null(bin_var)))
  {
    if(is.null(over_bin))
    {
      over_bin<-rep(0,times=bin_var)
    }
    else{if(length(over_bin)!=bin_var){stop("Please make sure that the overlap parameters are of the same length as the binary variables")}}
    
    bin_data<-vector()
    for(k in 1:bin_var)
    {
      cov_bin<-diag((.5), nrow=2)
      mean_bin<-rep(0,times=2)
      bin_data_temp<-vector()
      count<-1
      for(m in 1:clusters)
      {
        bin_lab<-rep(0, times=clus_dat_part[m])
        mean_bin[count]<-c(2.6)
        mean_bin[-count]<-(2.6*over_bin[k])
        data_bin<-mvtnorm::rmvnorm(clus_dat_part[m], mean=mean_bin, sigma=cov_bin)
        
        for(i in 1:clus_dat_part[m])
        {
          bin_lab[i]<-which(data_bin[i,]==max(data_bin[i,]))-1
        }
        bin_data_temp<-c(bin_data_temp, bin_lab)
        mean_bin<-rep(0,times=2)
        count<-count+1
        if(count>2)
        {count<-1}
      }
      bin_data<-c(bin_data, bin_data_temp)
      
    }
    tot_data<-c(tot_data,bin_data)
    
  }
  
  #generate nominal Data 
  if(!(is.null(nom_var)))
  {
    if(is.null(over_nom))
    {
      over_nom<-rep(0,times=nom_var)
    }
    else{if(length(over_nom)!=nom_var){stop("Please make sure that the overlap parameters are of the same length as the nominal variables")}}
    
    nom_data<-vector()
    for(k in 1:nom_var)
    {
      cov_ord<-diag((.5), nrow=2*clusters)
      mean_ord<-rep(0,times=2*clusters)
      nom_data_temp<-vector()
      for(m in 1:clusters)
      {
        nom_lab<-rep(0, times=clus_dat_part[m])
        mean_ord[(2*m-1):(2*m)]<-c(2.5,2.5)
        mean_ord[-((2*m-1):(2*m))]<-(2.5*over_nom[k])
        data_nom<-mvtnorm::rmvnorm(clus_dat_part[m], mean=mean_ord, sigma=cov_ord)
        
        for(i in 1:clus_dat_part[m])
        {
          if(max(data_nom[i,])<0)
          {
            nom_lab[i]=1
          }
          else{
            nom_lab[i]<-which(data_nom[i,]==max(data_nom[i,]))+1
          }}
        nom_data_temp<-c(nom_data_temp, nom_lab)
        mean_ord<-rep(0,times=2*clusters)
      }
      nom_data<-c(nom_data, nom_data_temp)
      
    }
    tot_data<-c(tot_data,nom_data)
    
  }
  
  tot_data_mat<-matrix(tot_data,nrow=obser,byrow=F)
  result<-list(data=tot_data_mat)
  if(!(is.null(cont_var)))
  {
    result<-list(cont_means=mean_clus, data=tot_data_mat)
  }
  
  return(result)
}
