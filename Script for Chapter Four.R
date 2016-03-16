# ------------------------------------------------------------------------------
# Script 2: Chapter Four
# Chapter Title: Development of Robust Mean Squared Error Estimation of ELL Poverty Estimates  
# Description: Simulation Work based on Parametric Bootstrap procedure 
# considering Homoskedastic (HM) level-specific errors
# -------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
# -------------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------------
library(nlme)
library(mvtnorm)
library(foreach)
library(doMC)
registerDoMC(40)
# -------------------------------------------------------------------------------
# Functions 
# -------------------------------------------------------------------------------
tapply.order<-function(value,area,FUN,target.order){
  # value: The values
  # area: target area
  # target.order: The order u wish to see the outcome ....
  
  raw.output<-tapply(value,area,FUN)
  data<-data.frame(key=names(raw.output), value=raw.output)
  ordered.value<-data[match(target.order, data$key),]$value
  return(ordered.value)
}
# -------------------------------------------------------------------------------#
addTrans <- function(color,trans) {
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
# -------------------------------------------------------------------------------#
Var.Com.MM.2<-function(level.2,level.1,res.2,res.1){
  
  # Homoskedasticity at both levels
  # level.2: ID number of level 2
  # level.1: ID number of level 1
  # res.2: Cluster level residuals (average of marginal residuals of respective cluster)
  # res.1: HH level residuals (Marginal residuals of HHs / OLS residuals)
  
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.2<-as.vector(table(level.2))
  C<-length(n.2)
  ID.C<-unique(level.2)
  
  n<-length(level.1)
  
  n0.bar.2<-sum(n.2^2)/sum(n.2)
  s1.e<-sum((res.1-mean(res.1))^2)/(n-1)
  s2.e<-sum(n.2*(res.2-mean(res.1))^2)/(C-1)
  sigma2.2<-max(((n-1)*(C-1))/((n-n0.bar.2)*(n-C))*(-s1.e+s2.e),0)
  sigma2.1<-max((n-1)/(n-C)*s1.e-(C-1)/(n-C)*s2.e,0)
  result<-list(sigma2.1=sigma2.1,sigma2.2=sigma2.2)
  return(result)
  
}
# -------------------------------------------------------------------------------#
Var.Com.MM.3<-function(level.3,level.2,level.1,res.3,res.2,res.1){
  
  # Homoskedasticity at all three levels
  # level.3: ID number of level 3
  # level.2: ID number of level 2
  # level.1: ID number of level 1
  # res.2: Cluster level residuals (average of marginal residuals of respective cluster)
  # res.1: HH level residuals (Marginal residuals of HHs / OLS residuals)
  
  level.3<-as.numeric(factor(level.3))
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.3<-as.vector(table(level.3))
  n.2<-as.vector(table(level.2))
  n<-length(level.1)
  
  D<-length(n.3)
  C<-length(n.2)
  ID.D<-unique(level.3)
  ID.C<-unique(level.2)
  
  n0.bar.2<-sum(n.2^2)/sum(n.2)
  n0.bar.3<-sum(n.3^2)/sum(n.3)
  
  # Area wise calculation
  
  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  for (d in 1:D){
    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  }
  
  #  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  #  for (d in unique(level.3)){
  #    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
  #    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  #  }
  
  
  sum.n.0i<-sum(n.0i)
  
  s1.e<-sum((res.1-mean(res.1))^2)/(n-1)
  s2.e<-sum(n.2*(res.2-mean(res.1))^2)/(C-1)
  s3.e<-sum(n.3*(res.3-mean(res.1))^2)/(D-1)
  
  sigma2.1<-(n-1)*s1.e/(n-C)-(C-1)*s2.e/(n-C)
  sigma2.2<-(-(n-1)*(C-D)*s1.e+(C-1)*(n-D)*s2.e-(D-1)*(n-C)*s3.e)/((n-sum.n.0i)*(n-C))
  sigma2.3<-((n-1)*((C-1)*(sum.n.0i-n0.bar.2)-(D-1)*(n-n0.bar.2))*s1.e+
               (C-1)*((D-1)*(n-n0.bar.2)-(n-1)*(sum.n.0i-n0.bar.2))*s2.e+
               (D-1)*(n-n0.bar.2)*(n-C)*s3.e)/((n-n0.bar.3)*(n-sum.n.0i)*(n-C))
  
  
  result<-list(sigma2.1=sigma2.1,sigma2.2=sigma2.2,sigma2.3=sigma2.3)
  return(result)
  
}
# -------------------------------------------------------------------------------#
GLS.EST<-function(n.c,sigma2.ch,sigma2.c,x.matrix,y.s){
  # n.c is size of clusters
  # sigma2.ch is a vector (or constant) of variance components at HH level
  # sigma2.c is variance components at cluster level
  # Output: beta estimates with thier variance covariance matrix
  
  library(Matrix)
  n<-sum(n.c)
  number.cluster.s<-length(n.c)
  
  if (length(sigma2.ch)==1){
    sampleList <- list()
    for (i in 1:number.cluster.s) sampleList[[i]]<-diag(rep(sigma2.ch,n.c[i]))+sigma2.c*matrix(1,n.c[i],n.c[i])
  }
  
  if (length(sigma2.ch)>1){
    sampleList <- list()
    j<-1
    for (i in 1:number.cluster.s) {
      sampleList[[i]]<-diag(sigma2.ch[j:(j+n.c[i]-1)])+sigma2.c*matrix(1,n.c[i],n.c[i])
      j<-j+n.c[i]
    }
  }
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls))
}
# -------------------------------------------------------------------------------#
GLS.EST.2L<-function(level.2,level.1,sigma2.2,sigma2.1,x.matrix,y.s){
  # level.3: Area Level ID
  # level.2: Cluster Level ID
  # level.1: HH Level ID
  # sigma2.3, sigma2.2, sigma2.1 level specific VC (Homoskedastic here)
  # Output: beta estimates with thier variance covariance matrix
  
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  library(Matrix)
  n.2<-as.vector(table(level.2))
  n<-sum(n.2)
  number.cluster.s<-length(n.2)
  
  if (length(sigma2.1)==1){
    sampleList <- list()
    for (i in 1:number.cluster.s) sampleList[[i]]<-diag(rep(sigma2.1,n.2[i]))+sigma2.2*matrix(1,n.2[i],n.2[i])
  }
  
  
  if (length(sigma2.1)>1){
    sampleList <- list()
    j<-1
    for (i in 1:number.cluster.s) {
      sampleList[[i]]<-diag(sigma2.1[j:(j+n.2[i]-1)])+sigma2.2*matrix(1,n.2[i],n.2[i])
      j<-j+n.2[i]
    }
  }
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls),V=V)
}
# -------------------------------------------------------------------------------#
GLS.EST.3L<-function(level.3,level.2,level.1,sigma2.3,sigma2.2,sigma2.1,x.matrix,y.s){
  # level.3: Area Level ID
  # level.2: Cluster Level ID
  # level.1: HH Level ID
  # sigma2.3, sigma2.2, sigma2.1 level specific VC (Homoskedastic here)
  # Output: beta estimates with thier variance covariance matrix
  
  library(Matrix)
  
  level.3<-as.numeric(factor(level.3))
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.3<-as.vector(table(level.3))
  n.2<-as.vector(table(level.2))
  D<-length(n.3)
  
  sampleList <- list()
  
  for (d in 1:D){
    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
    n.i<-sum(n.i.j)
    
    if (length(sigma2.1)==1) {
      clusterlist<-list()
      for (i in 1:length(n.i.j)) clusterlist[[i]]<-matrix(1,n.i.j[i],n.i.j[i])
      
      sampleList[[d]]<-diag(rep(sigma2.1,n.i))+bdiag(clusterlist)*sigma2.2+sigma2.3*matrix(1,n.i,n.i)
    }
    
    if (length(sigma2.1)>1){
      clusterlist<-list()
      for (i in 1:length(n.i.j)) clusterlist[[i]]<-matrix(1,n.i.j[i],n.i.j[i])
      
      sampleList[[d]]<-diag(sigma2.1[level.1[level.3==d]])+bdiag(clusterlist)*sigma2.2+sigma2.3*matrix(1,n.i,n.i)
    }
    
  }
  
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls),V=V)
}
# -------------------------------------------------------------------------------#
rsquared.lmm.mom<-function(beta,Model.Matrix,est.sigma2){
  # Get design matrix of fixed effects from model
  Fmat <- Model.Matrix
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(t(beta)%*%t(Model.Matrix)))
  # Get variance of random effects by extracting variance components
  VarRand <- sum(as.numeric(unlist(est.sigma2)[-1]), na.rm=T)
  VarResid <- as.numeric(unlist(est.sigma2)[1])
  R.squared.M<-VarF/(VarF+VarRand+VarResid)
  R.squared.C<-(VarF+VarRand)/(VarF+VarRand+VarResid)
  list(VarF=VarF,VarRand=VarRand,VarResid=VarResid,R.squared.M=R.squared.M,R.squared.C=R.squared.C)
}
# -------------------------------------------------------------------------------#
Var.Com.MM.2.H<-function(level.2,level.1,res.2,con.res.1){
  # Heteroskedasticity at 1st level
  # Homoskedasticity at 2nd level
  
  # Homoskedasticity at both levels
  # level.2: ID number of level 2
  # level.1: ID number of level 1
  # res.2: Cluster level residuals (average of marginal residuals of respective cluster)
  # con.res.1: HH level conditional residuals (Marginal residuals of HHs - Cluster level residuals)
  
  
  n.2<-as.vector(table(level.2))
  n.1<-length(level.1)
  
  C<-length(n.2)
  ID.C<-unique(level.2)
  
  wc<-n.2/sum(n.2)
  tau.2.c<-tapply(con.res.1,level.2,var)/n.2
  sigma2.2<-max((sum(wc*(res.2-mean(res.2))^2)-sum(wc*(1-wc)*tau.2.c))/sum(wc*(1-wc)),0)
  
  
  result<-list(sigma2.2=sigma2.2)
  return(result)
  
}
# -------------------------------------------------------------------------------#
Var.Com.MM.3.H<-function(level.3,level.2,level.1,res.3,res.2,con.res.1){
  # Heteroskedasticity at 1st level
  # Homoskedasticity at 2nd and 3rd level
  # level.3  : ID of level 3
  # level.2  : ID of level 2
  # level.1  : ID of level 1
  # ID will be in ascending orders
  #   res.3  : Area level residuals (average of marginal residuals of respective area)
  #   res.2  : Cluster level residuals (average of marginal residuals of respective cluster)
  # con.res.1: HH level conditional residuals (Marginal residuals of HHs - Cluster level residuals - Area level residuals)
  
  level.3<-as.numeric(factor(level.3))
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.3<-as.vector(table(level.3))
  n.2<-as.vector(table(level.2))
  n.1<-length(level.1)
  
  D<-length(n.3)
  C<-length(n.2)
  ID.D<-unique(level.3)
  ID.C<-unique(level.2)
  
  # Number of cluster per area
  n.2.3<-NULL
  for (d in 1: D) n.2.3<-c(n.2.3,length(unique(level.2[level.3==d])))
  ID.C.D<-rep(ID.D,n.2.3)
  
  
  #  n0.bar.2<-sum(n.2^2)/sum(n.2)
  #  n0.bar.3<-sum(n.3^2)/sum(n.3)
  
  # Area wise calculation
  
  #  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  #  for (d in 1:D){
  #    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
  #    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  #  }
  # sum.n.0i<-sum(n.0i)
  
  
  
  #  s1.e<-sum((res.1-mean(res.1))^2)/(n-1)
  #  s2.e<-sum(n.2*(res.2-mean(res.1))^2)/(C-1)
  #  s3.e<-sum(n.3*(res.3-mean(res.1))^2)/(D-1)
  
  tau.2.c<-tapply(con.res.1,level.2,var)/n.2
  tau.2.d<-tapply(n.2^2*tau.2.c,ID.C.D,sum)/n.3^2
  wc<-n.2/sum(n.2)
  wd<-n.3/sum(n.3)
  
  
  sigma2.2<-max( ( sum(wc*(res.2-mean(res.2))^2) - sum(wd*(res.3-mean(res.3))^2) - sum(wc*(1-wc)*tau.2.c) + sum(wd*(1-wd)*tau.2.d) ) /
                   (sum(wc) - sum (1/wd*tapply(wc^2,ID.C.D,sum)))  , 0)
  
  sigma2.3<-max(
    ( sum( (1/wd-1)*tapply(wc^2,ID.C.D,sum) ) * sum(wc*(1-wc)*tau.2.c)
      - sum(wc*(1-wc)) * sum(wd*(1-wd)*tau.2.d)
      - sum( (1/wd-1)*tapply(wc^2,ID.C.D,sum) ) * sum(wc*(res.2-mean(res.2))^2)
      + sum(wc*(1-wc)) * sum(wd*(res.3-mean(res.3))^2) )  /
      ( sum(wd*(1-wd)) * (sum(wc) - sum (1/wd*tapply(wc^2,ID.C.D,sum))) )
  )
  
  
  result<-list(sigma2.2=sigma2.2,sigma2.3=sigma2.3)
  return(result)
  
}
# -------------------------------------------------------------------------------#
FGT.alpha<-function(y,z,alpha){ 
  # Function for FGT indicators # z: Poverty line
  
  if (length(z)==1){
    t.z=ifelse(y<z,1,0)
    t.z.alpha=t.z*((rep(1,length(y))-y/z)^alpha)
    povert=sum(t.z.alpha)/length(y)
  }
  
  if (length(z)>1){
    
    povert<-rep(0,length(z))
    
    for (i in 1:length(z)){
      t.z=ifelse(y<z[i],1,0)
      t.z.alpha=t.z*((rep(1,length(y))-y/z[i])^alpha)
      povert[i]=sum(t.z.alpha)/length(y)}
  }
  
  povert
}
# -------------------------------------------------------------------------------#
Population<-function(Area,Cluster.Area,HH.Cluster,Mu,Sigma,X,Model=c("Normal","Log-Normal")){
  
  # This function creats population for
  # Area=Number of Areas
  # Cluster.Area=Vector(Number of clusters per area)
  # HH.Cluster=Fixed Number of HHs per cluster
  # Mu=Population Mean
  # Sigma=Variance of random errors of Two/Three-level Model
  
  No.Area<-Area
  No.Cluster<-sum(Cluster.Area)
  N.Area<-Cluster.Area*HH.Cluster
  N.Cluster<-rep(HH.Cluster,Cluster.Area)
  N=sum(N.Area)
  
  if (is.null(X))  X.design<-cbind(rep(1,N))
  if (!is.null(X)) X.design<-cbind(rep(1,N),X)
  
  if(length(Sigma)==2){
    e1=rnorm(N,0,sqrt(Sigma[1]))
    e2=rnorm(No.Cluster,0,sqrt(Sigma[2]))
    e.2<-rep(e2,N.Cluster)
    y.ijk<-Mu+e1+e.2
    if (Model=="Normal") y.ijk<-X.design%*%Mu+e1+e.2
    if (Model=="Log-Normal") y.ijk<-(X.design%*%Mu+e1+e.2)
  }
  
  if(length(Sigma)==3){
    e1=rnorm(N,0,sqrt(Sigma[1]))
    e2=rnorm(No.Cluster,0,sqrt(Sigma[2]))
    e3=rnorm(No.Area,0,sqrt(Sigma[3]))
    e.3<-rep(e3,N.Area)
    e.2<-rep(e2,N.Cluster)
    if (Model=="Normal") y.ijk<-X.design%*%Mu+e1+e.2+e.3
    if (Model=="Log-Normal") y.ijk<-(X.design%*%Mu+e1+e.2+e.3)
  }
  
  # ID of Area, Cluster & HHs
  
  ID.Area<-c(1:No.Area)
  ID.Cluster<-c(1:No.Cluster)
  ID.HH<-c(1:N)
  
  # Construction of Population Data
  
  ID.D=rep(c(1:Area),N.Area)
  ID.EA.D=rep(ID.Cluster,N.Cluster)
  if (Model=="Normal") pop.data<-data.frame(ID.D,ID.EA.D,ID.HH,y.ijk)
  if (Model=="Log-Normal") pop.data<-data.frame(ID.D,ID.EA.D,ID.HH,y.ijk,X)
  # Main population census
  return(pop.data)
}
# -------------------------------------------------------------------------------#
Sample<-function(Pop.Data,Cluster.Area.s,HH.Cluster.s){
  # This function darw sample in two-stage for every area 
  # Pop.Data<-data.frame(Area ID, Cluster ID, HH ID, Y.ijk)
  # Cluster.Area.s<-Vector(Number of Sampled cluster per area)
  # HH.Cluster.s<-No. of HHs per Cluster
  
  # Selection of cluster area wise
  ID.Cluster.s<-NULL
  for (d in unique(Pop.Data[,1])){
    ID.Cluster.s<-c(ID.Cluster.s,sample(unique(Pop.Data[,2][Pop.Data[,1]==d]),Cluster.Area.s[d]))
  }
  
  # Selection of individuals cluster wise
  ID.HH.s<-NULL
  for (c in ID.Cluster.s){
    ID.HH.s<-c(ID.HH.s,sample(unique(Pop.Data[,3][Pop.Data[,2]==c]),HH.Cluster.s))
  }
  
  # Selection of the sampling units from the population
  
  data.s<-Pop.Data[ID.HH.s,] ; dim(data.s)
  row.names(data.s) <- NULL 
  
  return(data.s)
}
# -------------------------------------------------------------------------------#
# Estimation of Variance component under a three-level null model
# -------------------------------------------------------------------------------#
mme<-function(Sample.Data){
  
  # Function for Estimating Variance Component: Method of Moments --------------------------------------------
  
  # Parameters
  n.3<-tapply(Sample.Data[,3],Sample.Data[,1],length)
  n.2<-tapply(Sample.Data[,3],Sample.Data[,2],length)
  n<-sum(n.3)
  D<-length(unique(Sample.Data[,1]))
  C<-length(unique(Sample.Data[,2]))
  n0.bar.2<-sum(n.2^2)/n
  n0.bar.3<-sum(n.3^2)/n
  
  # Area wise calculation 
  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  for (d in 1:D){
    n.i.j<-tapply(Sample.Data[,3][Sample.Data[,1]==d],Sample.Data[,2][Sample.Data[,1]==d],length)
    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  }
  sum.n.0i<-sum(n.0i)
  
  # Means
  mean.s<-mean(Sample.Data[,4])
  cluster.mean.s<-tapply(Sample.Data[,4],Sample.Data[,2],mean)
  area.mean.s<-tapply(Sample.Data[,4],Sample.Data[,1],mean)
  
  # Variances
  s1.3<-var(Sample.Data[,4])
  s2.3<-sum(as.vector(n.2)*(as.vector(cluster.mean.s)-mean.s)^2)/(C-1)
  s3.3<-sum(as.vector(n.3)*(as.vector(area.mean.s)-mean.s)^2)/(D-1)
  
  # variance component under three-level model
  lemda.1.3<-(n-1)*s1.3/(n-C)-(C-1)*s2.3/(n-C)
  lemda.2.3<-(-(n-1)*(C-D)*s1.3+(C-1)*(n-D)*s2.3-(D-1)*(n-C)*s3.3)/((n-sum.n.0i)*(n-C))
  lemda.3.3<-((n-1)*((C-1)*(sum.n.0i-n0.bar.2)-(D-1)*(n-n0.bar.2))*s1.3+
                (C-1)*((D-1)*(n-n0.bar.2)-(n-1)*(sum.n.0i-n0.bar.2))*s2.3+
                (D-1)*(n-n0.bar.2)*(n-C)*s3.3)/((n-n0.bar.3)*(n-sum.n.0i)*(n-C))
  lemda.3.3[lemda.3.3<0]=0 
  est.lemda.3<-c(lemda.1.3,lemda.2.3,lemda.3.3)
  
  # variance component under two-level model
  lemda.1.2<-(n-1)/(n-C)*s1.3-(C-1)/(n-C)*s2.3
  lemda.2.2<-(n-1)*(C-1)/((n-n0.bar.2)*(n-C))*(-s1.3+s2.3)
  est.lemda.2<-c(lemda.1.2,lemda.2.2)
  
  result<-list(lemda.3=est.lemda.3,lemda.2=est.lemda.2,size.cluster=n.2,size.area=n.3,area=D,cluster=C,n.0i=n.0i)
  # names(result)<-c("lemda.3","lemda.2","size.cluster","size.area","area","cluster","n.0i")
  return(result)
}
# -------------------------------------------------------------------------------#
Var.Com.MM.2<-function(level.2,level.1,res.2,res.1){
  
  # Homoskedasticity at both levels
  # level.2: ID number of level 2
  # level.1: ID number of level 1
  # res.2: Cluster level residuals (average of marginal residuals of respective cluster)
  # res.1: HH level residuals (Marginal residuals of HHs / OLS residuals)
  
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.2<-as.vector(table(level.2))
  C<-length(n.2)
  ID.C<-unique(level.2)
  
  n<-length(level.1)
  
  n0.bar.2<-sum(n.2^2)/sum(n.2)
  s1.e<-sum((res.1-mean(res.1))^2)/(n-1)
  s2.e<-sum(n.2*(res.2-mean(res.1))^2)/(C-1)
  sigma2.2<-max(((n-1)*(C-1))/((n-n0.bar.2)*(n-C))*(-s1.e+s2.e),0)
  sigma2.1<-max((n-1)/(n-C)*s1.e-(C-1)/(n-C)*s2.e,0)
  result<-list(sigma2.1=sigma2.1,sigma2.2=sigma2.2)
  return(result)
  
}
# -------------------------------------------------------------------------------#
Var.Com.MM.3<-function(level.3,level.2,level.1,res.3,res.2,res.1){
  
  # Homoskedasticity at all three levels
  # level.3: ID number of level 3
  # level.2: ID number of level 2
  # level.1: ID number of level 1
  # res.2: Cluster level residuals (average of marginal residuals of respective cluster)
  # res.1: HH level residuals (Marginal residuals of HHs / OLS residuals)
  
  level.3<-as.numeric(factor(level.3))
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.3<-as.vector(table(level.3))
  n.2<-as.vector(table(level.2))
  n<-length(level.1)
  
  D<-length(n.3)
  C<-length(n.2)
  ID.D<-unique(level.3)
  ID.C<-unique(level.2)
  
  n0.bar.2<-sum(n.2^2)/sum(n.2)
  n0.bar.3<-sum(n.3^2)/sum(n.3)
  
  # Area wise calculation 
  
  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  for (d in 1:D){
    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  }
  
  #  n.0i<-NULL # sum(n.2^2)/sum(n.2)
  #  for (d in unique(level.3)){
  #    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
  #    n.0i<-c(n.0i,sum(n.i.j^2)/sum(n.i.j))
  #  }  
  
  
  sum.n.0i<-sum(n.0i)
  
  s1.e<-sum((res.1-mean(res.1))^2)/(n-1)
  s2.e<-sum(n.2*(res.2-mean(res.1))^2)/(C-1)
  s3.e<-sum(n.3*(res.3-mean(res.1))^2)/(D-1)
  
  sigma2.1<-(n-1)*s1.e/(n-C)-(C-1)*s2.e/(n-C)
  sigma2.2<-(-(n-1)*(C-D)*s1.e+(C-1)*(n-D)*s2.e-(D-1)*(n-C)*s3.e)/((n-sum.n.0i)*(n-C))
  sigma2.3<-((n-1)*((C-1)*(sum.n.0i-n0.bar.2)-(D-1)*(n-n0.bar.2))*s1.e+
               (C-1)*((D-1)*(n-n0.bar.2)-(n-1)*(sum.n.0i-n0.bar.2))*s2.e+
               (D-1)*(n-n0.bar.2)*(n-C)*s3.e)/((n-n0.bar.3)*(n-sum.n.0i)*(n-C))
  
  
  result<-list(sigma2.1=sigma2.1,sigma2.2=sigma2.2,sigma2.3=sigma2.3)
  return(result)
  
}
# -------------------------------------------------------------------------------#
GLS.EST.2L<-function(level.2,level.1,sigma2.2,sigma2.1,x.matrix,y.s){
  # level.3: Area Level ID
  # level.2: Cluster Level ID
  # level.1: HH Level ID
  # sigma2.3, sigma2.2, sigma2.1 level specific VC (Homoskedastic here)
  # Output: beta estimates with thier variance covariance matrix
  
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  library(Matrix)
  n.2<-as.vector(table(level.2))
  n<-sum(n.2)
  number.cluster.s<-length(n.2)
  
  if (length(sigma2.1)==1){
    sampleList <- list()
    for (i in 1:number.cluster.s) sampleList[[i]]<-diag(rep(sigma2.1,n.2[i]))+sigma2.2*matrix(1,n.2[i],n.2[i])
  }
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls)) 
}
# -------------------------------------------------------------------------------#
GLS.EST.2L<-function(level.2,level.1,sigma2.2,sigma2.1,x.matrix,y.s){
  # level.3: Area Level ID
  # level.2: Cluster Level ID
  # level.1: HH Level ID
  # sigma2.3, sigma2.2, sigma2.1 level specific VC (Homoskedastic here)
  # Output: beta estimates with thier variance covariance matrix
  
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  library(Matrix)
  n.2<-as.vector(table(level.2))
  n<-sum(n.2)
  number.cluster.s<-length(n.2)
  
  if (length(sigma2.1)==1){
    sampleList <- list()
    for (i in 1:number.cluster.s) sampleList[[i]]<-diag(rep(sigma2.1,n.2[i]))+sigma2.2*matrix(1,n.2[i],n.2[i])
  }
  
  
  if (length(sigma2.1)>1){
    sampleList <- list()
    j<-1
    for (i in 1:number.cluster.s) {
      sampleList[[i]]<-diag(sigma2.1[j:(j+n.2[i]-1)])+sigma2.2*matrix(1,n.2[i],n.2[i])   
      j<-j+n.2[i]
    }
  }
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls),V=V) 
}
# -------------------------------------------------------------------------------#
GLS.EST.3L<-function(level.3,level.2,level.1,sigma2.3,sigma2.2,sigma2.1,x.matrix,y.s){
  # level.3: Area Level ID
  # level.2: Cluster Level ID
  # level.1: HH Level ID
  # sigma2.3, sigma2.2, sigma2.1 level specific VC (Homoskedastic here)
  # Output: beta estimates with thier variance covariance matrix
  
  library(Matrix)
  
  level.3<-as.numeric(factor(level.3))
  level.2<-as.numeric(factor(level.2))
  level.1<-as.numeric(factor(level.1))
  
  n.3<-as.vector(table(level.3))
  n.2<-as.vector(table(level.2))
  D<-length(n.3)
  
  sampleList <- list()
  
  for (d in 1:D){ 
    n.i.j<-tapply(level.1[level.3==d],level.2[level.3==d],length)
    n.i<-sum(n.i.j)
    
    if (length(sigma2.1)==1) {
      clusterlist<-list()
      for (i in 1:length(n.i.j)) clusterlist[[i]]<-matrix(1,n.i.j[i],n.i.j[i])
      
      sampleList[[d]]<-diag(rep(sigma2.1,n.i))+bdiag(clusterlist)*sigma2.2+sigma2.3*matrix(1,n.i,n.i)
    }
    
    if (length(sigma2.1)>1){
      clusterlist<-list()
      for (i in 1:length(n.i.j)) clusterlist[[i]]<-matrix(1,n.i.j[i],n.i.j[i])
      
      sampleList[[d]]<-diag(sigma2.1[level.1[level.3==d]])+bdiag(clusterlist)*sigma2.2+sigma2.3*matrix(1,n.i,n.i)   
    }
    
  }
  
  
  V<-bdiag(sampleList)
  inv.V<-solve(V,sparse=TRUE)
  xtx<-(t(x.matrix)%*%inv.V%*%(x.matrix))
  xty<-(t(x.matrix)%*%inv.V%*%(y.s))
  beta.gls <- solve(xtx)%*%xty
  vcov.beta.gls<- solve(t(x.matrix)%*%(x.matrix)) %*% (t(x.matrix)%*%V%*%(x.matrix)) %*% solve(t(x.matrix)%*%(x.matrix))
  list(beta.gls=as.matrix(beta.gls),vcov.beta.gls=as.matrix(vcov.beta.gls),V=V) 
}
# -------------------------------------------------------------------------------#
ELL.PB.HM<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (length(t)==1){
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
  }
  
  if (length(t)>1){
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.Mean<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  F11<-tapply(y.l,ID.D,mean)
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.FGT<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t){
  
  # This function is for estimating FGT estimates
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # y.l is in original scale
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  if (! is.null(X.U)) z.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
  
  y.l<-exp(z.l) # y.l is in original scale
  
  if (length(t)==1){
    F00<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))
    F22<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))
  }
  
  if (length(t)>1){
    F00<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F22<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F00=F00,F11=F11,F22=F22)
  
}
# -------------------------------------------------------------------------------#
Cons.ELL.PB.HM<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(D,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.d)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.d)+eps.l
  
  if (length(t)==1){
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
  }
  
  if (length(t)>1){
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
Cons.ELL.PB.HM.Mean<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(D,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.d)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.d)+eps.l
  
  F11<-tapply(y.l,ID.D,mean)
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
Cons.ELL.PB.HM.FGT<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t){
  
  # This function is for estimating FGT estimates
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # y.l is in original scale
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  eta.l<-rnorm(D,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.d)+eps.l
  if (! is.null(X.U)) z.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.d)+eps.l # z.l is in logarithm scale
  
  y.l<-exp(z.l) # y.l is in original scale
  
  if (length(t)==1){
    F00<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))
    F22<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))
  }
  
  if (length(t)>1){
    F00<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F22<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F00=F00,F11=F11,F22=F22)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.3L<- function(beta,var.beta,var.com.1,var.com.2,var.com.3,ID.D,ID.C,X.U,t){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  u.l<-rnorm(D,0,sqrt(var.com.3))
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  if (length(t)==1){
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
  }
  
  if (length(t)>1){
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.3L.Mean<- function(beta,var.beta,var.com.1,var.com.2,var.com.3,ID.D,ID.C,X.U){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  u.l<-rnorm(D,0,sqrt(var.com.3))
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  F11<-tapply(y.l,ID.D,mean)
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.3L.FGT<- function(beta,var.beta,var.com.1,var.com.2,var.com.3,ID.D,ID.C,X.U,t){
  
  # This function is for estimating FGT estimates
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # y.l is in original scale
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  u.l<-rnorm(D,0,sqrt(var.com.3))
  eta.l<-rnorm(C,0,sqrt(var.com.2))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) z.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
  
  y.l<-exp(z.l) # y.l is in original scale
  
  if (length(t)==1){
    F00<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))
    F22<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))
  }
  
  if (length(t)>1){
    F00<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F22<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F00=F00,F11=F11,F22=F22)
  
}
# -------------------------------------------------------------------------------#
MELL.PB.HM.2L<- function(beta,var.beta,var.com.1,var.com.2,NoClusterBlock,ID.D,ID.C,X.U,t){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  # var.com.2: Should be vector of length - number of Block
  # NoClusterBlock: Should be a vector of length - number of Block
  
  NoBlock<-length(var.com.2)
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  adj.eta.l<-list() 
  
  for(i in 1:NoBlock) adj.eta.l[[i]]<-rnorm(NoClusterBlock[i],0,sqrt(var.com.2[i]))
  
  eta.l<-c(unlist(adj.eta.l))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (length(t)==1){
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
  }
  
  if (length(t)>1){
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F11=F11)
  
}
# -------------------------------------------------------------------------------#
MELL.PB.HM.2L.Mean<- function(beta,var.beta,var.com.1,var.com.2,NoClusterBlock,ID.D,ID.C,X.U){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  # var.com.2: Should be vector of length - number of Block
  # NoClusterBlock: Should be a vector of length - number of Block
  
  NoBlock<-length(var.com.2)
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  adj.eta.l<-list() 
  
  for(i in 1:NoBlock) adj.eta.l[[i]]<-rnorm(NoClusterBlock[i],0,sqrt(var.com.2[i]))
  
  eta.l<-c(unlist(adj.eta.l))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) y.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) y.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  F11<-tapply(y.l,ID.D,mean)
  
  list(F11=F11)
  
  
}
# -------------------------------------------------------------------------------#
MELL.PB.HM.2L.FGT<- function(beta,var.beta,var.com.1,var.com.2,NoClusterBlock,ID.D,ID.C,X.U,t){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  # var.com.2: Should be vector of length - number of Block
  # NoClusterBlock: Should be a vector of length - number of Block
  
  NoBlock<-length(var.com.2)
  
  N<-length(ID.D)
  
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
  
  adj.eta.l<-list() 
  
  for(i in 1:NoBlock) adj.eta.l[[i]]<-rnorm(NoClusterBlock[i],0,sqrt(var.com.2[i]))
  
  eta.l<-c(unlist(adj.eta.l))
  eps.l<-rnorm(N,0,sqrt(var.com.1))
  
  if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  if (! is.null(X.U)) z.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
  
  y.l<-exp(z.l) # y.l is in original scale
  
  
  if (length(t)==1){
    F00<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))
    F11<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))
    F22<-tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))
  }
  
  if (length(t)>1){
    F00<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,0))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F11<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,1))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
    F22<-array(matrix(simplify2array(tapply(y.l,ID.D,function (x) FGT.alpha(x,t,2))),nrow=length(unique(ID.D)),ncol=length(t),byrow=TRUE),
               dim=c(1,length(unique(ID.D)),length(t)))
  }
  
  list(F00=F00,F11=F11,F22=F22)
  
}
# -------------------------------------------------------------------------------#
MCS.MELL<-function(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal","Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("MM","REML"),
   Parameter=c("Mean","DF","FGT"),ELL.Method=c("ELL.2L","Opt.ELL.2L","Cons.ELL.2L","MELL1.2L","MELL2.2L","MELL3.2L","ELL.3L"),quant,NoSim,No.Boot){
  # MCS.MELL : Monte Carlo Simulation 
  # MELL : Modified ELL
  # t: the targeted Quantile (Single value or Vector)
  # Block.ID: Strata/Block ID for the area
  
  # Parallel Object Creation (Inner loop)   ==============================================================
  
  # f2 is for single Parameter (here Mean and DF for single quantile)
  f2 <- function(obj1,obj2) {
    z <- list(
      F11.True=rbind(obj1$F11.True,obj2$F11.True)
    )
  }
  
  # f2.array is for multiple Parameter (here DF for multiple quantile)
  f2.array <- function(obj1,obj2) {
    z <- list(
      F11.True=abind(obj1$F11.True,obj2$F11.True,along=1)
    )
  }
  
  # f2.FGT is for FGT at single quantile
  f2.FGT <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      F11.FGT0=rbind(obj1$F11.FGT0,obj2$F11.FGT0),
      F11.FGT1=rbind(obj1$F11.FGT1,obj2$F11.FGT1),
      F11.FGT2=rbind(obj1$F11.FGT2,obj2$F11.FGT2)
    )
  }
  
  # f2.array.FGT is for FGT at multiple quantiles
  f2.array.FGT <- function(obj1,obj2) {
    # f2.array is for multiple quantile
    z <- list(
      F11.FGT0=abind(obj1$F11.FGT0,obj2$F11.FGT0,along=1),
      F11.FGT1=abind(obj1$F11.FGT1,obj2$F11.FGT1,along=1),
      F11.FGT2=abind(obj1$F11.FGT2,obj2$F11.FGT2,along=1)
    )
  }
  
  # Estimated variance components: 2-L & 3-L
  est.lemda.3<-matrix(0,NoSim,3) ; est.lemda.2<-matrix(0,NoSim,2) 
  
  # Population Area Mean: True, ELL, Optimistic ELL, Conservative ELL, Adjusted ELL  =====================
  
  
  if (Parameter=="Mean"){
    
    Mean.True<-matrix(0,NoSim,Area);
    Mean.ELL<-matrix(0,NoSim,Area); 
    Mean.ELL.MSE<-matrix(0,NoSim,Area); 
    Mean.CR.I<-matrix(0,NoSim,Area)
    
    # Mean.Opt.ELL<-matrix(0,NoSim,Area); Mean.Con.ELL<-matrix(0,NoSim,Area); Mean.Adj.ELL.1<-matrix(0,NoSim,Area); Mean.Adj.ELL.2<-matrix(0,NoSim,Area); Mean.Adj.ELL.3<-matrix(0,NoSim,Area)
    # Mean.Opt.ELL.var<-matrix(0,NoSim,Area); Mean.Con.ELL.var<-matrix(0,NoSim,Area); Mean.Adj.ELL.1.var<-matrix(0,NoSim,Area); Mean.Adj.ELL.2.var<-matrix(0,NoSim,Area); Mean.Adj.ELL.3.var<-matrix(0,NoSim,Area)
    
  }
  
  # Population Area-specific DF & FGT: True, ELL, Optimistic ELL, Conservative ELL, Adjusted ELL  ========
  
  No.Q=length(quant)
  
  if (Parameter=="DF"){
    
    if (No.Q==1) {
      
      DF.True<-array(0,dim=c(NoSim,No.Area))
      DF.ELL<-array(0,dim=c(NoSim,No.Area))
      DF.ELL.MSE<-array(0,dim=c(NoSim,No.Area))
      DF.CR.I<-array(0,dim=c(NoSim,No.Area))    
      
      # Estimate of Population Parameter: FGT 
      # DF.Opt.ELL<-array(0,dim=c(NoSim,No.Area)) # DF.Con.ELL<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.1<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.2<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.3<-array(0,dim=c(NoSim,No.Area))    
      
      # Estimated Variance  of Population Parameter: FGT 
      # DF.Opt.ELL.var<-array(0,dim=c(NoSim,No.Area)) # DF.Con.ELL.var<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.1.var<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.2.var<-array(0,dim=c(NoSim,No.Area)) # DF.Adj.ELL.3.var<-array(0,dim=c(NoSim,No.Area))    
      
    }
    
    if (No.Q>1) {
      
      # True Population Parameter 
      
      DF.True<-array(0,dim=c(NoSim,No.Area,No.Q))
      DF.ELL<-array(0,dim=c(NoSim,No.Area,No.Q))
      DF.ELL.MSE<-array(0,dim=c(NoSim,No.Area,No.Q)) 
      DF.CR.I<-array(0,dim=c(NoSim,No.Area,No.Q))    
      
      # Estimate of Population Parameter: FGT 
      # DF.ELL<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Opt.ELL<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Con.ELL<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Adj.ELL.1<-array(0,dim=c(NoSim,No.Area,No.Q))
      # DF.Adj.ELL.2<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Adj.ELL.3<-array(0,dim=c(NoSim,No.Area,No.Q))    
      
      # Estimated Variance  of Population Parameter: FGT 
      # DF.ELL.var<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Opt.ELL.var<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Con.ELL.var<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Adj.ELL.1.var<-array(0,dim=c(NoSim,No.Area,No.Q))
      # DF.Adj.ELL.2.var<-array(0,dim=c(NoSim,No.Area,No.Q)) # DF.Adj.ELL.3.var<-array(0,dim=c(NoSim,No.Area,No.Q))    
      
    }
    
  }
  
  if (Parameter=="FGT"){
    
    if (No.Q==1) {
      
      # True Population Parameter 
      
      F0.True<-array(0,dim=c(NoSim,No.Area))     
      F1.True<-array(0,dim=c(NoSim,No.Area))     
      F2.True<-array(0,dim=c(NoSim,No.Area))      
      
      F0.F11<-array(0,dim=c(NoSim,No.Area))       # ELL estimates FGT 0
      F1.F11<-array(0,dim=c(NoSim,No.Area))       # ELL estimates FGT 1
      F2.F11<-array(0,dim=c(NoSim,No.Area))       # ELL estimates FGT 2
      
      F0.F11.MSE<-array(0,dim=c(NoSim,No.Area))   # MSE of estimated FGT 0
      F1.F11.MSE<-array(0,dim=c(NoSim,No.Area))   # MSE of estimated FGT 1
      F2.F11.MSE<-array(0,dim=c(NoSim,No.Area))   # MSE of estimated FGT 2
      
      F0.CR.I<-array(0,dim=c(NoSim,No.Area))      # Coverage Indicator of 95% CI  
      F1.CR.I<-array(0,dim=c(NoSim,No.Area))      # Coverage Indicator of 95% CI  
      F2.CR.I<-array(0,dim=c(NoSim,No.Area))      # Coverage Indicator of 95% CI  
      
    }
    
    if (No.Q>1) {
      
      # True Population Parameter 
      
      F0.True<-array(0,dim=c(NoSim,No.Area,No.Q))      # True estimates: FGT 0
      F1.True<-array(0,dim=c(NoSim,No.Area,No.Q))      # True estimates: FGT 1
      F2.True<-array(0,dim=c(NoSim,No.Area,No.Q))      # True estimates: FGT 2
      
      F0.F11<-array(0,dim=c(NoSim,No.Area,No.Q))       # ELL estimates FGT 0
      F1.F11<-array(0,dim=c(NoSim,No.Area,No.Q))       # ELL estimates FGT 1
      F2.F11<-array(0,dim=c(NoSim,No.Area,No.Q))       # ELL estimates FGT 2
      
      F0.F11.MSE<-array(0,dim=c(NoSim,No.Area,No.Q))   # MSE of estimated FGT 0
      F1.F11.MSE<-array(0,dim=c(NoSim,No.Area,No.Q))   # MSE of estimated FGT 1
      F2.F11.MSE<-array(0,dim=c(NoSim,No.Area,No.Q))   # MSE of estimated FGT 2
      
      F0.CR.I<-array(0,dim=c(NoSim,No.Area,No.Q))      # Coverage Indicator of 95% CI  
      F1.CR.I<-array(0,dim=c(NoSim,No.Area,No.Q))      # Coverage Indicator of 95% CI  
      F2.CR.I<-array(0,dim=c(NoSim,No.Area,No.Q))      # Coverage Indicator of 95% CI  
      
    }
    
  }
  
  # Store the significant test results: Variance component
  
  test.results<-rep(0,NoSim)
  
  # Number of block/Strata to split the whole population -------------------------------------------
  
  NoBlock<-length(unique(Block.ID)) 
  N<-sum(Cluster.Area*HH.Cluster)
  # Simulation start from Here ===========================================================================
  
  for(s in 1:NoSim){
    
    cat(date(),"Iteration number",s,"starting","\n",fill=T)
    
    if (Model==c("Normal")) X=NULL
    
    if (Model==c("Log-Normal")) {
      Sigma.Mu<-c(0.5,0.75)
      Sigma2<-matrix(c(1.50,0.10,0.10,0.95),2,2,byrow=TRUE)
      X<-rmvnorm(N,Sigma.Mu,Sigma2)
    }
    
    # Population and Sample generation -------------------------------------------------------------------------
    #Pop.Data<-Population(Area,Cluster.Area,HH.Cluster,Mu,Sigma) 
    if (Model==c("Normal")) Pop.Data<-Population(Area,Cluster.Area,HH.Cluster,Mu,Sigma,X=NULL,Model=c("Normal"))
    if (Model==c("Log-Normal")) Pop.Data<-Population(Area,Cluster.Area,HH.Cluster,Mu,Sigma,X,Model=c("Log-Normal"))
    
    Sample.Data<-Sample(Pop.Data,Cluster.Area.s,HH.Cluster.s)
    if (Model==c("Normal")) z<-quantile(Pop.Data$y.ijk,quant)
    if (Model==c("Log-Normal") & Parameter==c("FGT")) z<-quantile(exp(Pop.Data$y.ijk),quant)
    # Population & Sample parameters -------------------------------------------------------------------
    No.Cluster<-length(unique(Pop.Data$ID.EA.D))
    No.Cluster.s<-length(unique(Sample.Data$ID.EA.D))
    N.Area<-tapply(Pop.Data$y.ijk,Pop.Data$ID.D,length)
    N.Cluster<-tapply(Pop.Data$y.ijk,Pop.Data$ID.EA.D,length)
    N<-sum(N.Area)
    n<-length(Sample.Data$y.ijk)
    n.i<-tapply(Sample.Data$ID.HH,Sample.Data$ID.D,length)
    n.ij<-tapply(Sample.Data$ID.HH,Sample.Data$ID.EA.D,length)
    
    if (is.null(X)) {
      lm1<-lm(y.ijk~1,Sample.Data)
      X.ols<-model.matrix(lm1)
    }
    if (!is.null(X)) {
      lm1<-lm(y.ijk~X1+X2,Sample.Data)
      X.ols<-model.matrix(lm1)
    }
    
    # Estimation of VC, regression parametrs and their var-cov matrix: MM ------------------------------    
    
    if (Var.Com=="MM") {
      # Needs to be change ================#
      
      #    fit.mm<-mme(Sample.Data) # Only for mean model
      #    est.lemda.3[s,]<-fit.mm$lemda.3
      #    est.lemda.2[s,]<-fit.mm$lemda.2
      
      Sample.Data$u.ch<-resid(lm1)
      eta<-as.vector(tapply(Sample.Data$u.ch,Sample.Data$ID.D,mean))
      u<-as.vector(tapply(Sample.Data$u.ch,Sample.Data$ID.EA.D,mean))
      Sample.Data$eps<-Sample.Data$u.ch-rep(eta,n.i)-rep(u,n.ij) # No application latter
      Sample.Data$y.hat<-fitted.values(lm1)
      
      est.sigma2.2<-Var.Com.MM.2(Sample.Data$ID.EA.D,Sample.Data$ID.HH,u,Sample.Data$u.ch)
      est.sigma2.3<-Var.Com.MM.3(Sample.Data$ID.D,Sample.Data$ID.EA.D,Sample.Data$ID.HH,eta,u,Sample.Data$u.ch)
      
      est.lemda.2[s,]<-c(est.sigma2.2$sigma2.1,est.sigma2.2$sigma2.2)
      est.lemda.3[s,]<-c(est.sigma2.3$sigma2.1,est.sigma2.3$sigma2.2,est.sigma2.3$sigma2.3)
      
      gls.lm2<-GLS.EST.2L(Sample.Data$ID.EA.D,Sample.Data$ID.HH,est.sigma2.2$sigma2.2,est.sigma2.2$sigma2.1,x.matrix=X.ols,Sample.Data$y.ijk)
      beta.gls.2<-gls.lm2$beta.gls
      var.beta.gls.2<-gls.lm2$vcov.beta.gls
      
      # I have to modify this
      gls.lm3<-GLS.EST.3L(Sample.Data$ID.D,Sample.Data$ID.EA.D,Sample.Data$ID.HH,est.sigma2.3$sigma2.3,est.sigma2.3$sigma2.2,est.sigma2.3$sigma2.1,x.matrix=X.ols,Sample.Data$y.ijk)
      beta.gls.3<-gls.lm3$beta.gls
      var.beta.gls.3<-gls.lm3$vcov.beta.gls
      
      # Need to check whether
    }
    
    # Estimation of VC, regression parametrs and their var-cov matrix: REML ------------------------------
    
    if (Var.Com=="REML") {
      
      if (Model=="Normal") fit.lme.2<-lme(fixed=y.ijk~1, random=~1|ID.EA.D,data=Sample.Data)
      if (Model=="Log-Normal") fit.lme.2<-lme(fixed=y.ijk~X1+X2, random=~1|ID.EA.D,data=Sample.Data)
      
      lemda.2.lme<-nlme::VarCorr(fit.lme.2)
      est.lemda.2[s,]<-as.numeric(c(lemda.2.lme[2,1],lemda.2.lme[1,1]))
      beta.gls.2<-fixed.effects(fit.lme.2)
      var.beta.gls.2<-summary(fit.lme.2)$varFix
      
      if (Model=="Normal") fit.lme.3<-lme(fixed=y.ijk~1, random=~1|ID.D/ID.EA.D,data=Sample.Data)
      if (Model=="Log-Normal") fit.lme.3<-lme(fixed=y.ijk~X1+X2, random=~1|ID.D/ID.EA.D,data=Sample.Data)
      
      lemda.3.lme<-nlme::VarCorr(fit.lme.3)
      est.lemda.3[s,]<-as.numeric(c(lemda.3.lme[5,1],lemda.3.lme[4,1],lemda.3.lme[2,1]))
      beta.gls.3<-fixed.effects(fit.lme.3)
      var.beta.gls.3<-summary(fit.lme.3)$varFix
      test.results[s]<-anova.lme(fit.lme.2,fit.lme.3)$ 'p-value'[2]
    }
    
    # names(summary(fit.lme.3)) # to see the objests under summary(lm.object)
    
    # Block ID for every HHs
    
    ID.Block.U<-rep(Block.ID,N.Area)
    
    # Adjustment of the 2nd variance component under 2-L  ------------------------------------------------
    
    Ki.s<-matrix(0,Area,1) ## Adjustment factor by area from sample structure
    Ki.P<-matrix(0,Area,1) ## Adjustment factor by area from population structure
    
    for (d in 1:Area){
      nij<-tapply(Sample.Data$ID.HH[Sample.Data$ID.D==d],Sample.Data$ID.EA.D[Sample.Data$ID.D==d],length)
      Ki.s[d,]<-sum(nij)^2/sum(nij^2) ## Adjustment factor by area
      Nij<-tapply(Pop.Data$ID.HH[Pop.Data$ID.D==d],Pop.Data$ID.EA.D[Pop.Data$ID.D==d],length)
      Ki.P[d,]<-sum(Nij)^2/sum(Nij^2) ## Adjustment factor by area from population structure
    }
    
    # Blockwise parameter ------------------------------------------------------------------------
    
    Ki.Block<-list(); for(i in 1:NoBlock) Ki.Block[[i]]<-Ki.P[unique(Pop.Data$ID.D[ID.Block.U==i])]
    Cluster.B<-list(); for(i in 1:NoBlock) Cluster.B[[i]]<-unique(Pop.Data$ID.EA.D[ID.Block.U==i])
    No.Cluster.B<-list(); for(i in 1:NoBlock) No.Cluster.B[[i]]<-length(unique(Pop.Data$ID.EA.D[ID.Block.U==i]))
    N.B<-list(); for(i in 1:NoBlock) N.B[[i]]<-length(Pop.Data$ID.EA.D[ID.Block.U==i])
    
    # Blockwise calculations by manually ---------------------------------------------------------
    
    # for(i in 1:NoBlock) assign(paste("Ki.Block", i, sep = "."),Ki.P[unique(Pop.Data$ID.D[ID.Block.U==i])])
    # Ki.Block.1<-Ki.P[01:15] ; Ki.Block.2<-Ki.P[16:30]; Ki.Block.3<-Ki.P[31:45]; Ki.Block.4<-Ki.P[46:60]; Ki.Block.5<-Ki.P[61:75]
    # for(i in 1:NoBlock) assign(paste("Cluster.B", i, sep = ""),unique(Pop.Data$ID.EA.D[ID.Block.U==i]))
    # Cluster.B1<-unique(Pop.Data$ID.EA.D[Pop.Data$ID.D<=15]); Cluster.B2<-unique(Pop.Data$ID.EA.D[Pop.Data$ID.D>=16 & Pop.Data$ID.D<=30]); Cluster.B3<-unique(Pop.Data$ID.EA.D[Pop.Data$ID.D>=31 & Pop.Data$ID.D<=45])
    # Cluster.B4<-unique(Pop.Data$ID.EA.D[Pop.Data$ID.D>=46 & Pop.Data$ID.D<=60]); Cluster.B5<-unique(Pop.Data$ID.EA.D[Pop.Data$ID.D>=61 & Pop.Data$ID.D<=75])
    
    # for(i in 1:NoBlock) assign(paste("No.Cluster.B", i, sep = ""),length(unique(Pop.Data$ID.EA.D[ID.Block.U==i])))
    # No.Cluster.B1<-length(Cluster.B1); No.Cluster.B2<-length(Cluster.B2); No.Cluster.B3<-length(Cluster.B3); No.Cluster.B4<-length(Cluster.B4); No.Cluster.B5<-length(Cluster.B5)
    
    # for(i in 1:NoBlock) assign(paste("N.B", i, sep = ""),length(Pop.Data$ID.EA.D[ID.Block.U==i]))
    # N.B1<-length(Pop.Data[,2][Pop.Data$ID.D<=15]); N.B2<-length(Pop.Data[,2][Pop.Data$ID.D>=16 & Pop.Data$ID.D<=30]); N.B3<-length(Pop.Data[,2][Pop.Data$ID.D>=31 & Pop.Data$ID.D<=45])
    # N.B4<-length(Pop.Data[,2][Pop.Data$ID.D>=46 & Pop.Data$ID.D<=60]); N.B5<-length(Pop.Data[,2][Pop.Data$ID.D>=61 & Pop.Data$ID.D<=75])
    
    # All together in a single function: Cluster ID by Block; No. of Cluster by block; No. of HH by block
    
    # for(i in 1:NoBlock){
    #   assign(paste("Ki.Block", i, sep = "."),Ki.P[unique(Pop.Data$ID.D[ID.Block.U==i])])              # Area-specific weight by block
    #   assign(paste("Cluster.B", i, sep = ""),unique(Pop.Data$ID.EA.D[ID.Block.U==i]))                 # Block specific cluster
    #   assign(paste("No.Cluster.B", i, sep = ""),length(unique(Pop.Data$ID.EA.D[ID.Block.U==i])))      # Block specific number of cluster  
    #   assign(paste("N.B", i, sep = ""),length(Pop.Data$ID.EA.D[ID.Block.U==i]))                       # Block size
    # }
    
    
    # Adjustment factors --------------------------------------------------------------------------
    
    # 1st adjustment factor
    K1<-(est.lemda.3[s,2]+mean(Ki.s)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    adj.lemda2.2.K1<-K1*est.lemda.2[s,][2]
    
    # 2nd adjustment factor
    K2<-(est.lemda.3[s,2]+mean(Ki.P)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    adj.lemda2.2.K2<-K2*est.lemda.2[s,][2]
    
    # 3rd adjustment factor: Block wise
    K3<-list(); for(i in 1:NoBlock) K3[[i]]<-(est.lemda.3[s,2]+mean(Ki.Block[[i]])*est.lemda.3[s,3])/est.lemda.2[s,][2]
    adj.lemda2.2.K3.i<-list(); for(i in 1:NoBlock) adj.lemda2.2.K3.i[[i]]<-K3[[i]]*est.lemda.2[s,][2]
    
    # Manual calculation of the Weights ========================================================================================
    # 3rd adjustment factor: Block wise
    # K3.Block.1<-(est.lemda.3[s,2]+mean(Ki.Block.1)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    # adj.lemda2.2.K3.1<-K3.Block.1*est.lemda.2[s,][2]
    
    # K3.Block.2<-(est.lemda.3[s,2]+mean(Ki.Block.2)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    # adj.lemda2.2.K3.2<-K3.Block.2*est.lemda.2[s,][2]
    
    # K3.Block.3<-(est.lemda.3[s,2]+mean(Ki.Block.3)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    # adj.lemda2.2.K3.3<-K3.Block.3*est.lemda.2[s,][2]
    
    # K3.Block.4<-(est.lemda.3[s,2]+mean(Ki.Block.4)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    # adj.lemda2.2.K3.4<-K3.Block.4*est.lemda.2[s,][2]
    
    # K3.Block.5<-(est.lemda.3[s,2]+mean(Ki.Block.5)*est.lemda.3[s,3])/est.lemda.2[s,][2]
    # adj.lemda2.2.K3.5<-K3.Block.5*est.lemda.2[s,][2]
    
    # adj.lemda2.2.K3<-mean(adj.lemda2.2.K3.1,adj.lemda2.2.K3.2,adj.lemda2.2.K3.3,adj.lemda2.2.K3.4,adj.lemda2.2.K3.5)
    
    # Application of the ELL and MELL Functions ================================================================
    
    if (Parameter=="Mean") {
      
      if (ELL.Method=="ELL.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      if (ELL.Method=="Opt.ELL.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True) }
      }
      
      if (ELL.Method=="Cons.ELL.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-Cons.ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      if (ELL.Method=="MELL1.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      if (ELL.Method=="MELL2.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      if (ELL.Method=="MELL3.2L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-MELL.PB.HM.2L.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      if (ELL.Method=="ELL.3L") {
        r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
          F11.True<-ELL.PB.HM.3L.Mean(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X)$F11
          list(F11.True=F11.True)}
      }
      
      #  ELL.2L<-ELL.PB.HM.Mean(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  Opt.ELL.2L<-ELL.PB.HM.Mean(Opt.beta.gls.2,Opt.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  Cons.ELL.2L<-Cons.ELL.PB.HM.Mean(Cons.beta.gls.2,Cons.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  MELL1.2L<-ELL.PB.HM.Mean(beta.gls.k1,var.beta.gls.k1,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  MELL2.2L<-ELL.PB.HM.Mean(beta.gls.k2,var.beta.gls.k2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  MELL3.2L<-MELL.PB.HM.2L.Mean(beta.gls.k3,var.beta.gls.k3,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      #  ELL.3L<-ELL.PB.HM.3L.Mean(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL)$F11
      
      Mean.True[s,]<-tapply(Pop.Data$y.ijk,Pop.Data$ID.D,mean)
      
      Mean.ELL[s,]<-colMeans(r.FGT$F11.True)
      Mean.ELL.MSE[s,]<-apply(r.FGT$F11.True,2,sd) 
      F11.Q.02.5<-apply(r.FGT$F11.True,2,function(x) quantile(x,0.025,na.rm=TRUE))
      F11.Q.97.5<-apply(r.FGT$F11.True,2,function(x) quantile(x,0.975,na.rm=TRUE))
      Mean.CR.I[s,]<-(Mean.True[s,]>=F11.Q.02.5 & Mean.True[s,]<=F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix  
      
    }
    
    if (Parameter=="DF") {
      
      if (No.Q==1) {
        
        if (ELL.Method=="ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="Opt.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True) }
        }
        
        if (ELL.Method=="Cons.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-Cons.ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL1.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL2.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL3.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-MELL.PB.HM.2L(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="ELL.3L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
            F11.True<-ELL.PB.HM.3L(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        DF.True[s,]<-tapply(Pop.Data$y.ijk,Pop.Data$ID.D,function(x) FGT.alpha(x,z,0))
        
        DF.ELL[s,]<-colMeans(r.FGT$F11.True)
        DF.ELL.MSE[s,]<-apply(r.FGT$F11.True,2,sd) 
        F11.Q.02.5<-apply(r.FGT$F11.True,2,function(x) quantile(x,0.025,na.rm=TRUE))
        F11.Q.97.5<-apply(r.FGT$F11.True,2,function(x) quantile(x,0.975,na.rm=TRUE))
        DF.CR.I[s,]<-(DF.True[s,]>=F11.Q.02.5 & DF.True[s,]<=F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix  
      }
      
      #  ELL.2L<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      #  Opt.ELL.2L<-ELL.PB.HM(Opt.beta.gls.2,Opt.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      #  Con.ELL.2L<-Cons.ELL.PB.HM(Cons.beta.gls.2,Cons.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      #  MELL1.2L<-ELL.PB.HM(beta.gls.k1,var.beta.gls.k1,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      #  MELL2.2L<-ELL.PB.HM(beta.gls.k2,var.beta.gls.k2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      #  MELL3.2L<-MELL.PB.HM.2L(beta.gls.k3,var.beta.gls.k3,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)$F11
      #  ELL.3L<-ELL.PB.HM.3L(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U==NULL,z)
      
      if (No.Q>1) {
        
        if (ELL.Method=="ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="Opt.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True) }
        }
        
        if (ELL.Method=="Cons.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-Cons.ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL1.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL2.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-ELL.PB.HM(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="MELL3.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-MELL.PB.HM.2L(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        if (ELL.Method=="ELL.3L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array) %dopar% {
            F11.True<-ELL.PB.HM.3L(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)$F11
            list(F11.True=F11.True)}
        }
        
        DF.True[s,,]<-array(matrix(simplify2array(tapply(Pop.Data$y.ijk,Pop.Data$ID.D,function (x) FGT.alpha(x,z,0))),nrow=No.Area,ncol=length(z),byrow=TRUE),dim=c(1,No.Area,length(z))) # True estimates: FGT 0
        
        DF.ELL[s,,]<-colMeans(r.FGT$F11.True,dims=1)
        DF.ELL.MSE[s,,]<-apply(r.FGT$F11.True,c(2,3),sd)
        F11.Q.02.5<-apply(r.FGT$F11.True,c(2,3),function(x) quantile(x,0.025,na.rm=TRUE))
        F11.Q.97.5<-apply(r.FGT$F11.True,c(2,3),function(x) quantile(x,0.975,na.rm=TRUE))
        DF.CR.I[s,,]<-(DF.True[s,,]>=F11.Q.02.5 & DF.True[s,,]<=F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
      }
      
    }
    
    if (Parameter=="FGT") {
      
      if (No.Q==1) {
        
        if (ELL.Method=="ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="Opt.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="Cons.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-Cons.ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL1.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL2.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL3.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-MELL.PB.HM.2L.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="ELL.3L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.FGT) %dopar% {
            F11.True<-ELL.PB.HM.3L.FGT(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        F0.True[s,]<-tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function(x) FGT.alpha(x,z,0))
        F1.True[s,]<-tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function(x) FGT.alpha(x,z,1))
        F2.True[s,]<-tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function(x) FGT.alpha(x,z,2))
        
        F0.F11[s,]<-colMeans(r.FGT$F11.FGT0) ; F0.F11.MSE[s,]<-apply(r.FGT$F11.FGT0,2,sd)
        F0.F11.Q.02.5<-apply(r.FGT$F11.FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
        F0.F11.Q.97.5<-apply(r.FGT$F11.FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
        F0.CR.I[s,]<-(F0.True[s,]>=F0.F11.Q.02.5 & F0.True[s,]<=F0.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
        
        F1.F11[s,]<-colMeans(r.FGT$F11.FGT1) ; F1.F11.MSE[s,]<-apply(r.FGT$F11.FGT1,2,sd)
        F1.F11.Q.02.5<-apply(r.FGT$F11.FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
        F1.F11.Q.97.5<-apply(r.FGT$F11.FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
        F1.CR.I[s,]<-(F1.True[s,]>=F1.F11.Q.02.5 & F1.True[s,]<=F1.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
        
        F2.F11[s,]<-colMeans(r.FGT$F11.FGT2) ; F2.F11.MSE[s,]<-apply(r.FGT$F11.FGT2,2,sd)
        F2.F11.Q.02.5<-apply(r.FGT$F11.FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
        F2.F11.Q.97.5<-apply(r.FGT$F11.FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
        F2.CR.I[s,]<-(F2.True[s,]>=F2.F11.Q.02.5 & F2.True[s,]<=F2.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
      }
      
      #  ELL.2L<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      #  Opt.ELL.2L<-ELL.PB.HM.FGT(Opt.beta.gls.2,Opt.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      #  Con.ELL.2L<-Cons.ELL.PB.HM.FGT(Cons.beta.gls.2,Cons.var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      #  MELL1.2L<-ELL.PB.HM.FGT(beta.gls.k1,var.beta.gls.k1,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      #  MELL2.2L<-ELL.PB.HM.FGT(beta.gls.k2,var.beta.gls.k2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      #  MELL3.2L<-MELL.PB.HM.2L.FGT(beta.gls.k3,var.beta.gls.k3,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)$F11
      #  ELL.3L<-ELL.PB.HM.3L.FGT(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X==NULL,z)
      
      if (No.Q>1) {
        
        if (ELL.Method=="ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=est.lemda.2[s,2],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="Opt.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="Cons.ELL.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-Cons.ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=(est.lemda.3[s,2]+est.lemda.3[s,3]),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL1.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K1,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL2.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-ELL.PB.HM.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=adj.lemda2.2.K2,Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="MELL3.2L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-MELL.PB.HM.2L.FGT(beta.gls.2,var.beta.gls.2,var.com.1=est.lemda.2[s,1],var.com.2=c(unlist(adj.lemda2.2.K3.i)),NoClusterBlock=c(unlist(No.Cluster.B)),Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        if (ELL.Method=="ELL.3L") {
          r.FGT <- foreach(icount(No.Boot), .combine=f2.array.FGT) %dopar% {
            F11.True<-ELL.PB.HM.3L.FGT(beta.gls.3,var.beta.gls.3,var.com.1=est.lemda.3[s,1],var.com.2=est.lemda.3[s,2],var.com.3=est.lemda.3[s,3],Pop.Data$ID.D,Pop.Data$ID.EA.D,X.U=X,z)
            F11.FGT0<-F11.True$F00; F11.FGT1<-F11.True$F11; F11.FGT2<-F11.True$F22
            F11.True<-NULL
            list(F11.FGT0=F11.FGT0,F11.FGT1=F11.FGT1,F11.FGT2=F11.FGT2)
          }
        }
        
        F0.True[s,,]<-array(matrix(simplify2array(tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function (x) FGT.alpha(x,z,0))),nrow=No.Area,ncol=length(z),byrow=TRUE),dim=c(1,No.Area,length(z))) # True estimates: FGT 0
        F1.True[s,,]<-array(matrix(simplify2array(tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function (x) FGT.alpha(x,z,1))),nrow=No.Area,ncol=length(z),byrow=TRUE),dim=c(1,No.Area,length(z))) # True estimates: FGT 1
        F2.True[s,,]<-array(matrix(simplify2array(tapply(exp(Pop.Data$y.ijk),Pop.Data$ID.D,function (x) FGT.alpha(x,z,2))),nrow=No.Area,ncol=length(z),byrow=TRUE),dim=c(1,No.Area,length(z))) # True estimates: FGT 2
        
        F0.F11[s,,]<-colMeans(r.FGT$F11.FGT0,dims=1) ; F0.F11.MSE[s,,]<-apply(r.FGT$F11.FGT0,c(2,3),sd)        
        F0.F11.Q.02.5<-apply(r.FGT$F11.FGT0,c(2,3),function(x) quantile(x,0.025,na.rm=TRUE))
        F0.F11.Q.97.5<-apply(r.FGT$F11.FGT0,c(2,3),function(x) quantile(x,0.975,na.rm=TRUE))
        F0.CR.I[s,,]<-(F0.True[s,,]>=F0.F11.Q.02.5 & F0.True[s,,]<=F0.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
        
        F1.F11[s,,]<-colMeans(r.FGT$F11.FGT1,dims=1) ; F1.F11.MSE[s,,]<-apply(r.FGT$F11.FGT1,c(2,3),sd)
        F1.F11.Q.02.5<-apply(r.FGT$F11.FGT1,c(2,3),function(x) quantile(x,0.025,na.rm=TRUE))
        F1.F11.Q.97.5<-apply(r.FGT$F11.FGT1,c(2,3),function(x) quantile(x,0.975,na.rm=TRUE))
        F1.CR.I[s,,]<-(F1.True[s,,]>=F1.F11.Q.02.5 & F1.True[s,,]<=F1.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
        
        F2.F11[s,,]<-colMeans(r.FGT$F11.FGT2,dims=1) ; F2.F11.MSE[s,,]<-apply(r.FGT$F11.FGT2,c(2,3),sd)
        F2.F11.Q.02.5<-apply(r.FGT$F11.FGT2,c(2,3),function(x) quantile(x,0.025,na.rm=TRUE))
        F2.F11.Q.97.5<-apply(r.FGT$F11.FGT2,c(2,3),function(x) quantile(x,0.975,na.rm=TRUE))
        F2.CR.I[s,,]<-(F2.True[s,,]>=F2.F11.Q.02.5 & F2.True[s,,]<=F2.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
        
      }
      
    }
    
  } # Simulation End
  
  if (Parameter==c("Mean")) Results<-list(Sigma2.HM.3l=est.lemda.3,Sigma2.HM.2l=est.lemda.2,test.results=test.results,Mean.True=Mean.True, Mean.ELL=Mean.ELL, Mean.ELL.MSE=Mean.ELL.MSE,Mean.CR.I=Mean.CR.I)
  
  if (Parameter==c("DF"))  Results<-list(Sigma2.HM.3l=est.lemda.3,Sigma2.HM.2l=est.lemda.2,test.results=test.results,DF.True=DF.True, DF.ELL=DF.ELL, DF.ELL.MSE=DF.ELL.MSE,DF.CR.I=DF.CR.I)
  
  if (Parameter==c("FGT")) Results<-list(Sigma2.HM.3l=est.lemda.3,Sigma2.HM.2l=est.lemda.2,test.results=test.results,F0.True=F0.True,F1.True=F1.True,F2.True=F2.True,
                                         F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.CR.I=F0.CR.I,F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.CR.I=F1.CR.I,F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.CR.I=F2.CR.I)
  
  return(Results )
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.2L.FGT.Estimator<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t,HH.Size,No.Boot){
  
  # Validation
  #beta<-beta.gls.2;var.beta<-var.beta.gls.2;var.com.1<-est.sigma2.2$sigma2.1;var.com.2<-est.sigma2.2$sigma2.2;
  #ID.D<-Census.data.fitted$ID.UPZ;ID.C<-Census.data.fitted$psu;
  #X.U<-x.matrix.U; t<-Census.data.fitted$lpovln;HH.Size<-Census.data.fitted$HH.size.U
  
  # This function is for estimating FGT estimates under ELL HM
  # Considering HH Size, No HH size put HH.Size=NULL
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # ID.D, ID.C,X.U,t,HH.Size all correspond to HH-level information
  # All variables are ordered according to ID.D
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  N.d<-as.vector(table(ID.D))
  C<-length(unique(ID.C))
  
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    eta.l<-rnorm(C,0,sqrt(var.com.2))
    eps.l<-rnorm(N,0,sqrt(var.com.1))
    
    if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
    if (! is.null(X.U)) z.l<-as.matrix(cbind(rep(1,N),X.U))%*%t(beta.l)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
      
    }
    
    if (! is.null(HH.Size)) {
      
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(Area.ID=unique(ID.D),
       F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
}
# -------------------------------------------------------------------------------#
ELL.PB.HM.3L.FGT.Estimator<- function(beta,var.beta,var.com.1,var.com.2,var.com.3,ID.D,ID.C,X.U,t,HH.Size,No.Boot){
  
  # This function is for estimating FGT estimates
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # y.l is in original scale
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    
    u.l<-rnorm(D,0,sqrt(var.com.3))
    eta.l<-rnorm(C,0,sqrt(var.com.2))
    eps.l<-rnorm(N,0,sqrt(var.com.1))
    
    if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
    if (! is.null(X.U)) z.l<-as.matrix(cbind(rep(1,N),X.U))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
      
    }
    
    if (! is.null(HH.Size)) {
      
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
}
# -------------------------------------------------------------------------------#
MELL.PB.HM.2L.FGT.Estimator<- function(beta,var.beta,var.com.1,Mod.var.com.2,NoClusterBlock,ID.D,ID.C,X.U,t,HH.Size,No.Boot){
  
  # This function is for estimating Distribution Function
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector
  # Mod.var.com.2: Should be vector of length - number of Block
  # NoClusterBlock: Should be a vector of length - number of Block
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  NoBlock<-length(Mod.var.com.2)
  
  N<-length(ID.D)
  
  #  N.c<-as.vector(table(ID.C))
  #  C<-length(unique(ID.C))
  
  N.c<-as.vector(table(ID.C))[unique(ID.C)]
  C<-length(unique(ID.C))
  
  
  N.d<-as.vector(table(ID.D))[unique(ID.D)]
  D<-length(unique(ID.D))
  
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    
    adj.eta.l<-list()
    
    for(i in 1:NoBlock) adj.eta.l[[i]]<-rnorm(NoClusterBlock[i],0,sqrt(Mod.var.com.2[i]))
    
    eta.l<-c(unlist(adj.eta.l))
    eps.l<-rnorm(N,0,sqrt(var.com.1))
    
    if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
    if (! is.null(X.U)) z.l<-as.matrix(cbind(rep(1,N),X.U))%*%t(beta.l)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
      
    }
    
    if (! is.null(HH.Size)) {
      
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
  
}
# -------------------------------------------------------------------------------#
MELL.NPB.HM.2L.FGT.Estimator<- function(beta,var.beta,eta.s,eps.s,var.com.1,var.com.2,Mod.var.com.2,NoClusterBlock,ID.C.s,ID.D,ID.C,X.U,t,HH.Size,No.Boot,residual.1=c("unconditional","conditional")){
  
  # Non-parametric Bootstrap procedure (ELL, 2002,2003)
  # This function is for estimating Poverty Indicators
  # Heteroskedastic error variances based on logistic function
  # t is here a value but need to be a vector
  # Conditional: Draw individual level effects for a population cluster from the sample cluster whose random effect is randomly drawn for the population cluster
  # Unconditional: Draw without restriction of drawing cluster random effects
  # eta.s should be ordered in ID.C.s
  # eps.s and sig2.eps.s should be ordered in ID.C.s
  # Sample and population data sets should be ordered by area and clusters
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  N<-length(ID.D)
  N.d<-as.vector(table(ID.D))[unique(ID.D)]
  N.c<-as.vector(table(ID.C))[unique(ID.C)]
  C<-length(unique(ID.C))
  # standardized the residuals using the individual specific residuals
  
  # eps.std<-eps.s/sqrt(sig2.eps.s)-mean(eps.s/sqrt(sig2.eps.s))
  
  scale.2<-sqrt(var.com.2/(sum(eta.s^2)/(length(eta.s)-1)))
  eta.s.scl<-eta.s*scale.2
  
  scale.1<-sqrt(var.com.1/(sum(eps.s^2)/(length(eps.s)-1)))
  eps.s.scl<-eps.s*scale.1
  
  # M2.eta.s<-u/sqrt(est.sigma2.2.HT)*sqrt(M2.est.sigma2.2)
  
  NoBlock<-length(Mod.var.com.2)
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    
    if (residual.1=="unconditional") {
      
      adj.eta.l<-list()
      
      #for(i in 1:NoBlock) adj.eta.l[[i]]<-sample(eta.s/sqrt(var.com.2)*sqrt(Mod.var.com.2[i]),NoClusterBlock[i],replace=TRUE)
      for(i in 1:NoBlock) adj.eta.l[[i]]<-sample(eta.s.scl/sqrt(var.com.2)*sqrt(Mod.var.com.2[i]),NoClusterBlock[i],replace=TRUE)
      
      eta.l<-c(unlist(adj.eta.l))
      eps.l<-sample(eps.s.scl,N,replace=TRUE)
    }
    
    if (residual.1=="conditional") {
      
      adj.eta.l<-list()
      
      ID.C.sampling<-sample(c(1:length(unique(ID.C.s))),C,replace=TRUE)
      eta.l.origin<-eta.s.scl[ID.C.sampling]
      
      for(i in 1:NoBlock) adj.eta.l[[i]]<-sample(eta.l.origin/sqrt(var.com.2)*sqrt(Mod.var.com.2[i]),NoClusterBlock[i],replace=TRUE)
      eta.l<-c(unlist(adj.eta.l))
      
      ID.C.sampled<-unique(ID.C.s)[ID.C.sampling]
      
      eps.l<-NULL
      for (c in 1:C){
        eps.U<-sample(eps.s.scl[ID.C.s==ID.C.sampled[c]],N.c[c],replace=TRUE)
        eps.l<-c(eps.l,eps.U)
      }
      
    }
    
    z.l<-cbind(rep(1,N),X.U)%*%t(beta.l)+rep(eta.l,N.c)+eps.l
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
    }
    
    if (! is.null(HH.Size)) {
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
}
# -------------------------------------------------------------------------------
simulation.Results.Robust.Work<-function(True.Value.object,Estimate.object,Estimate.object.MSE){
  
  No.Area<-dim(True.Value.object)[2]
  
  True.Area.Mean<-apply(True.Value.object,2,mean)
  True.Area.Mean.MSE<-apply(True.Value.object,2,var)
  
  Estimate.Mean<-apply(Estimate.object,2,mean)
  Estimate.Mean.MSE<-apply(Estimate.object.MSE,2,mean)
  
  RB.Estimate.Mean<-(Estimate.Mean-True.Area.Mean)/True.Area.Mean
  RB.Estimate.Mean.MSE<-(Estimate.Mean.MSE-True.Area.Mean.MSE)/True.Area.Mean.MSE
  
  Average.RB.Estimate.Mean<-mean(RB.Estimate.Mean)
  Average.RB.Estimate.Mean.MSE<-mean(RB.Estimate.Mean.MSE)
  
  ARB.Estimate.Mean<-abs(Estimate.Mean-True.Area.Mean)/True.Area.Mean
  ARB.Estimate.Mean.MSE<-abs(Estimate.Mean.MSE-True.Area.Mean.MSE)/True.Area.Mean.MSE
  
  Average.ARB.Estimate.Mean<-mean(ARB.Estimate.Mean)
  Average.ARB.Estimate.Mean.MSE<-mean(ARB.Estimate.Mean.MSE)
  
  RMSE.Estimate.Mean.MSE<-sqrt(colMeans((Estimate.object.MSE-matrix(True.Area.Mean.MSE,nrow=dim(Estimate.object.MSE)[1],ncol=dim(Estimate.object.MSE)[2],byrow=TRUE))^2))
  RRMSE.Estimate.Mean.MSE<-sqrt(colMeans((Estimate.object.MSE-matrix(True.Area.Mean.MSE,nrow=dim(Estimate.object.MSE)[1],ncol=dim(Estimate.object.MSE)[2],byrow=TRUE))^2)/True.Area.Mean.MSE)
  
  Average.RMSE.Estimate.Mean.MSE<-mean(RMSE.Estimate.Mean.MSE)
  Average.RRMSE.Estimate.Mean.MSE<-mean(RRMSE.Estimate.Mean.MSE)
  
  IND<-(True.Value.object>=(Estimate.object-1.96*sqrt(Estimate.object.MSE)) & True.Value.object<=(Estimate.object+1.96*sqrt(Estimate.object.MSE)))*1
  
  # IND<-ifelse(True.Value.object>=(Estimate.object-1.96*sqrt(Estimate.object.MSE)) & True.Value.object<=(Estimate.object+1.96*sqrt(Estimate.object.MSE)),1,0)
  
  
  CR<-apply(IND,2,mean)
  Average.CR<-mean(CR)
  
  CIW<-apply((Estimate.object+1.96*sqrt(Estimate.object.MSE))-(Estimate.object-1.96*sqrt(Estimate.object.MSE)),2,mean)
  Average.CIW<-mean(CIW)
  
  list(True.Area.Mean.MSE=True.Area.Mean.MSE,Estimate.Mean.MSE=Estimate.Mean.MSE,
       RB.Estimate.Mean=RB.Estimate.Mean,Average.RB.Estimate.Mean=Average.RB.Estimate.Mean,
       RB.Estimate.Mean.MSE=RB.Estimate.Mean.MSE,Average.RB.Estimate.Mean.MSE=Average.RB.Estimate.Mean.MSE,
       ARB.Estimate.Mean=ARB.Estimate.Mean,Average.ARB.Estimate.Mean=Average.ARB.Estimate.Mean,
       ARB.Estimate.Mean.MSE=ARB.Estimate.Mean.MSE,Average.ARB.Estimate.Mean.MSE=Average.ARB.Estimate.Mean.MSE,
       RMSE.Estimate.Mean.MSE=RMSE.Estimate.Mean.MSE,Average.RMSE.Estimate.Mean.MSE=Average.RMSE.Estimate.Mean.MSE,RRMSE.Estimate.Mean.MSE=RRMSE.Estimate.Mean.MSE,
       Average.RRMSE.Estimate.Mean.MSE=Average.RRMSE.Estimate.Mean.MSE,CR=CR,Average.CR=Average.CR,CIW=CIW,Average.CIW=Average.CIW)
  
}
# -------------------------------------------------------------------------------
# Simulation Study
# -------------------------------------------------------------------------------
# Population 1 : Three-Level mean model with Normal distributed random errors
# -------------------------------------------------------------------------------
Area=75
Cluster.Area=rep(c(15:29),each=5)
HH.Cluster=c(rep(c(96:100),3),rep(c(101:105),3),rep(c(106:110),3),rep(c(111:115),3),rep(c(116:120),3))
Mu=20
Sigma<-c(0.80,0.15,0.05)
No.Area<-Area
# -------------------------------------------------------------------------------#
Cluster.Area.s<-rep(c(2:4),25)
HH.Cluster.s<-10
D=c(1:Area)  
quant<-c(0.10,0.25,0.50,0.75,0.90)
Block.ID<-rep(1:5,each=15)
# -------------------------------------------------------------------------------
# Simulation for Population 1: Three-level Normal Population - Area specific Mean
# -------------------------------------------------------------------------------
set.seed(1000)
Sce1.Mean.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("Mean"),ELL.Method=c("ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.Opt.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                               Parameter=c("Mean"),ELL.Method=c("Opt.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.Cons.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                                Parameter=c("Mean"),ELL.Method=c("Cons.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.MELL1.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL1.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.MELL2.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL2.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.MELL3.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL3.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.Mean.ELL.3L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("Mean"),ELL.Method=c("ELL.3L"),quant,NoSim=5,No.Boot=5)
# -------------------------------------------------------------------------------
# Simulation for Population 1: Three-level Normal Population -  Area-specific DF
# -------------------------------------------------------------------------------
set.seed(1000)
Sce1.DF.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                         Parameter=c("DF"),ELL.Method=c("ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.Opt.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("DF"),ELL.Method=c("Opt.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.Cons.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                              Parameter=c("DF"),ELL.Method=c("Cons.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.MELL1.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL1.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.MELL2.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL2.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.MELL3.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL3.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce1.DF.ELL.3L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                         Parameter=c("DF"),ELL.Method=c("ELL.3L"),quant,NoSim=5,No.Boot=5)
# -------------------------------------------------------------------------------
# Population 2: Two-Level mean model with Normally distributed Random Errors
# -------------------------------------------------------------------------------
Area=75
Cluster.Area=rep(c(15:29),each=5)
HH.Cluster=c(rep(c(96:100),3),rep(c(101:105),3),rep(c(106:110),3),rep(c(111:115),3),rep(c(116:120),3))
Mu=20
Sigma<-c(0.80,0.20)
No.Area<-Area
# -------------------------------------------------------------------------------#
Cluster.Area.s<-rep(c(2:4),25)
HH.Cluster.s<-10
D=c(1:Area)  
quant<-c(0.10,0.25,0.50,0.75,0.90)
# -------------------------------------------------------------------------------
# Simulation for Population 2: Two-level Normal Population - NULL Model - Area specific Mean
# -------------------------------------------------------------------------------
set.seed(1000)
Sce2.Mean.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("Mean"),ELL.Method=c("ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.Opt.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                               Parameter=c("Mean"),ELL.Method=c("Opt.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.Cons.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                                Parameter=c("Mean"),ELL.Method=c("Cons.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.MELL1.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL1.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.MELL2.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL2.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.MELL3.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("Mean"),ELL.Method=c("MELL3.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.Mean.ELL.3L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("Mean"),ELL.Method=c("ELL.3L"),quant,NoSim=5,No.Boot=5)
# -------------------------------------------------------------------------------
# Simulation for Population 2: Two-level Normal Population - NULL MOdel - Area-specific DF
# -------------------------------------------------------------------------------
set.seed(1000)
Sce2.DF.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                         Parameter=c("DF"),ELL.Method=c("ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.Opt.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                             Parameter=c("DF"),ELL.Method=c("Opt.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.Cons.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                              Parameter=c("DF"),ELL.Method=c("Cons.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.MELL1.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL1.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.MELL2.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL2.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.MELL3.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                           Parameter=c("DF"),ELL.Method=c("MELL3.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce2.DF.ELL.3L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                         Parameter=c("DF"),ELL.Method=c("ELL.3L"),quant,NoSim=5,No.Boot=5)
# -------------------------------------------------------------------------------
# Population 3 : Three Level full model with log-Normal Population 
# -------------------------------------------------------------------------------
Area=75
Cluster.Area=rep(c(15:29),each=5)
HH.Cluster=c(rep(c(96:100),3),rep(c(101:105),3),rep(c(106:110),3),rep(c(111:115),3),rep(c(116:120),3))
Mu<-c(6,0.5,-0.55)
Sigma<-c(0.20,0.035,0.015)
#Sigma<-c(0.50,0.20,0.05)
N<-sum(Cluster.Area*HH.Cluster)
No.Area<-Area
# Generate Design Matrix (Fixed)
Sigma.Mu<-c(0.5,0.75)
Sigma2<-matrix(c(1.50,0.10,0.10,0.95),2,2,byrow=TRUE)
X<-rmvnorm(N,Sigma.Mu,Sigma2)
# -------------------------------------------------------------------------------#
Cluster.Area.s<-rep(c(2:4),25)
HH.Cluster.s<-10
D=c(1:Area)  
quant<-c(0.10,0.25)
# -------------------------------------------------------------------------------
# Simulation for Population 3: Three-level log-normal Population - Area-specific FGT
# -------------------------------------------------------------------------------
set.seed(1000)
Sce3.FGT.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                          Parameter=c("FGT"),ELL.Method=c("ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.Opt.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                              Parameter=c("FGT"),ELL.Method=c("Opt.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.Cons.ELL.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                               Parameter=c("FGT"),ELL.Method=c("Cons.ELL.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.MELL1.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                            Parameter=c("FGT"),ELL.Method=c("MELL1.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.MELL2.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                            Parameter=c("FGT"),ELL.Method=c("MELL2.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.MELL3.2L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                            Parameter=c("FGT"),ELL.Method=c("MELL3.2L"),quant,NoSim=5,No.Boot=5)
set.seed(1000)
Sce3.FGT.ELL.3L<-MCS.MELL(Area,Cluster.Area,HH.Cluster,Mu,Sigma,Model=c("Log-Normal"),Cluster.Area.s,HH.Cluster.s,Block.ID,Var.Com=c("REML"),
                          Parameter=c("FGT"),ELL.Method=c("ELL.3L"),quant,NoSim=5,No.Boot=5)
# -------------------------------------------------------------------------------
