FLCM.select<-function(y,mydata,nbasis1=7,nbasis2=15,cvl=100)
{
  GrandY<-y
  n=length(unique(mydata$id))
  T<-list()
  
  for(i in 1:n)
  {
    ind<-which(mydata$id==i)
    T[[i]]<-mydata$time[ind]
  }
  
  name<-names(mydata)
  name<-name[-c(1,2)]
  p<- ncol(mydata)-2
  UT<-c()

  for( i in 1:n)
  {
    UT<-union(UT,T[[i]])
  }
  UT<-sort(UT)

  library(fda)
  
  knots = seq(min(UT),max(UT),l=(nbasis2-2))
  norder = 4
  nbasis = length(knots) + norder - 2
  dayrng = c(min(UT),max(UT))
  bbasis = create.bspline.basis(dayrng,nbasis,norder,knots)
  
  #bbasis$nbasis    # number of basis functions
  #bbasis$rangeval   # basis range
  #plot(bbasis)
  
  bbasisMat<-list()

  for (i in 1:n){
    bbasisMat[[i]] = eval.basis(T[[i]],bbasis)
  }

  #dim(bbasisMat[[1]])
  in.mat = inprod(bbasis,bbasis) #penalty matrix for function
  #par(mfrow=c(1,1))
  #image(in.mat)
  in.derivmat=inprod(bbasis,bbasis,int2Lfd(2),int2Lfd(2)) #penalty matrix for 2nd deriv.
  #par(mfrow=c(1,1))
  #image(in.derivmat)
  
  basismat<-bbasisMat
  for (i in 1:n){
    colnames(basismat[[i]])<-NULL}
  
  XMatstar<-list()
  for(i in 1:n){
    XMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
    for (j in 1:p){
      for(l in 1:length(T[[i]])){
        indtemp<-which(mydata$id==i)
        datatemp<-mydata[indtemp,]
        datatempj<-datatemp[,2+j]
        XMatstar[[i]][j,l,]=datatempj[l]*basismat[[i]][l,]
      }
    }
  }
  knotsint = seq(min(UT),max(UT),l=(nbasis1-2))
  norderint = 4
  nbasisint = length(knotsint) + norder - 2
  dayrngint = c(min(UT),max(UT))
  bbasisint = create.bspline.basis(dayrngint,nbasisint,norderint,knotsint)
  bbasisMatint<-list()
  for (i in 1:n)
  {bbasisMatint[[i]] = eval.basis(T[[i]],bbasisint)}
  B<-bbasisMatint
  GrandB<-NULL
  for(i in 1:n)
  {GrandB<-rbind(GrandB,B[[i]])
  }
  group<-c()
  d<-nbasis-nbasisint
  for (j in 1:(p+1))
  {if(j==1)
  {group[((nbasisint*(j-1))+1):(nbasisint*j)]<-0}
    if(j>1) 
    {group[((nbasis*(j-1))+1-d):((nbasis*j)-d)]<-j-1}
  }
  ###getting prewihtning matrix
  ZMatstar<-list()
  for(i in 1:n)
  {ZMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
  for (j in 1:p){
    for(l in 1:length(T[[i]])){
      ZMatstar[[i]][j,l,]=XMatstar[[i]][j,l,]
    }
  }
  }
  Zimatlist<-list()
  for (i in 1:n)
  {Zimatlist[[i]]<-matrix(0,length(T[[i]]),nbasis*p)
  for(l in 1:length(T[[i]]))
  {for(j in 1 :p){
    Zimatlist[[i]][l,((nbasis*(j-1))+1):(nbasis*j)]<-ZMatstar[[i]][j,l,]
  }
  }
  }
  temp<-NULL
  for(i in 1:n)
  {temp<-rbind(temp, Zimatlist[[i]])
  }
  GrandZ<-temp
  GrandZ1<-cbind(GrandB,GrandZ)
  library(mgcv)
  fit <- bam(GrandY ~ -1 + GrandZ1,method = "REML")
  res <- fit$residuals

  lenT <- c()
  for (i in 1:n)
  {
    lenT[i] <- length(T[[i]])
  }
  .id = rep(1:n, lenT)
  .index <- unlist(T[1:n])
  resdata <- data.frame(.id , .index)
  resdata$.value <- res
  #Using FPCA to get estimate of the covariance matrix
  fpca1 <-  fpca.sc(ydata = resdata, center = TRUE, nbasis = 15, pve = 0.999, var = TRUE)
  timeU <- as.numeric(attr(fpca1$Yhat, "dimnames")[[2]])
  m <- length(timeU)
  if (length(fpca1$evalues) > 1)
  {
    sigmat <- fpca1$efunctions %*% diag(fpca1$evalues) %*% t(fpca1$efunctions) + fpca1$sigma2 *diag(m)
  }
  if (length(fpca1$evalues) == 1)
  {
    sigmat <- fpca1$evalues * fpca1$efunctions %*% t(fpca1$efunctions) + fpca1$sigma2 *diag(m)
  }
  indtime <- list()
  for (i in 1:n)
  {
    indtime[[i]] <- vector(mode = "logical", length = length(T[[i]]))
    for (l in 1:length(T[[i]])) {
      indtime[[i]][l] = which(round(timeU, 4) == round(T[[i]][l], 4))
    }
  }
  subsigmat <- list()
  for (i in 1:n)
  {
    subsigmat[[i]] <- sigmat[indtime[[i]], indtime[[i]]]
  }

  Linv <- list()
  for (i in 1:n)
  {
    R <- chol(subsigmat[[i]])
    L <- t(R)
    Linv[[i]] <- forwardsolve(L, diag(dim(L)[1]))
  }
  Linvhalf <- bdiag(Linv)


  #do crossvalidation for phi
  library(parallel)
  ncores <-detectCores()-1
  cl <- makeCluster(ncores)
  cvphi<-function(phi){
    #p=18
    #phi=10 ##How to choose phi later
    Kphi<-in.mat+phi*in.derivmat
    R<-chol(Kphi)
    L<-t(R)
    M<-solve(R)
    ZMatstar<-list()
    for(i in 1:n)
    {ZMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
    for (j in 1:p){
      for(l in 1:length(T[[i]])){
        ZMatstar[[i]][j,l,]=XMatstar[[i]][j,l,]%*%M
      }
    }
    }
    Zimatlist<-list()
    for (i in 1:n)
    {Zimatlist[[i]]<-matrix(0,length(T[[i]]),nbasis*p)
    for(l in 1:length(T[[i]]))
    {for(j in 1 :p){
      Zimatlist[[i]][l,((nbasis*(j-1))+1):(nbasis*j)]<-ZMatstar[[i]][j,l,]
    }
    }
    }
    GrandZ<-NULL
    for(i in 1:n)
    {GrandZ<-rbind(GrandZ, Zimatlist[[i]])
    }
    Grandzz<-cbind(GrandB,GrandZ)
    Grandzz <- as.matrix(Linvhalf %*% Grandzz)
    GrandY<- as.matrix(Linvhalf %*% GrandY)
    Groupvar<-as.factor(group)
    library(grpreg)
    fit2 <- grpreg(Grandzz,GrandY,Groupvar,penalty="grSCAD")
    fit3 <- grpreg(Grandzz,GrandY,Groupvar,penalty="grMCP")
    cvfit2<-select(fit2,crit="EBIC")
    cvfit3<-select(fit3,crit="EBIC")
    gamma2<-cvfit2$beta
    gamma3<-cvfit3$beta
    cvescad<-min(cvfit2$IC)
    cvemcp<-min(cvfit3$IC)
    if(!is.finite(cvescad))
    {cvescad<-999999}
    if(!is.finite(cvemcp))
    {cvemcp<-999999}
    lambdascad<-cvfit2$lambda
    lambdamcp<-cvfit3$lambda
    cvresult<-list(cvescad,cvemcp,lambdascad,lambdamcp,gamma2,gamma3,M)
    names(cvresult)<-c("cvescad","cvemcp","lambdascad","lambdamcp","gamma2","gamma3","M")
    return(cvresult)
  }
  clusterExport(cl,c())
  clusterEvalQ(cl, library("fda"))
  clusterEvalQ(cl, library("grpreg"))
  clusterExport(cl, list("in.mat","in.derivmat", "XMatstar", "T", "nbasis","n","p","GrandY","GrandB","bbasisMat","group","Linvhalf"), envir=environment())
  phiseq<-10^seq(-2,15,l=cvl)
  cvresultfinal<- parLapply(cl,phiseq,cvphi)
  stopCluster(cl)
  len<-length(phiseq)
  cves<-c()
  cvem<-c()
  for (k in 1:len)  
  {
    cves<-c(cves,cvresultfinal[[k]]$cvescad)
    cvem<-c(cvem,cvresultfinal[[k]]$cvemcp) 
  }
  indscad<-which.min((cves))
  indmcp<-which.min(unlist(cvem))
  gamma2<-cvresultfinal[[indscad]]$gamma2
  gamma3<-cvresultfinal[[indmcp]]$gamma3
  Mscad<-cvresultfinal[[indscad]]$M
  Mmcp<-cvresultfinal[[indmcp]]$M
  indsel2<-which(gamma2[-1]!=0)
  Groupsel2<-group[indsel2]
  grind2<-unique(Groupsel2)
  varselected2<-grind2
  indsel3<-which(gamma3[-1]!=0)
  Groupsel3<-group[indsel3]
  grind3<-unique(Groupsel3)
  varselected3<-grind3
  namescad<-name[varselected2[-1]]
  namemcp<-name[varselected3[-1]]
  result<-list(namescad,namemcp)
  #################return following objects too in order get estimate####################################
  ###result<-list(gamma2,gamma3,Mscad,Mmcp,bbasisMat,bbasis,bbasisint,bbasisMatint)
  names(result)<-c("scad","mcp")
  return(result)
}