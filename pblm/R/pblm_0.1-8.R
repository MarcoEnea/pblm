#da crreggere penalty.vlasso cn diag(lam1[-1]) che si blcca
require(Matrix)
require(splines)
require(MASS)

Mdiag <- function(x){ 
  if (is.null(x)) M <- NULL else
  {
    x <- x[!sapply(x,is.null)]
    dimlist <- sapply(x,dim)
    a <- apply(dimlist,1,cumsum) 
    dimMat <- rowSums(dimlist)
    M <- array(0,dimMat)    
    indexdim <- rbind(c(0,0),a)    
    for (i in 1:length(x))
      M[(indexdim[i,1]+1):indexdim[i+1,1],
        (indexdim[i,2]+1):indexdim[i+1,2]] <- x[[i]]
  } 
  M 
}

Mdiag2 <- function(x){ 
  dimlist <- sapply(x,dim)
  a <- apply(dimlist,1,cumsum) 
  dimMat <- rowSums(dimlist)  
  M <- array(0,dimMat)
  indexdim <- rbind(c(0,0),a)    
  for (i in 1:length(x))
    M[(indexdim[i,1]+1):indexdim[i+1,1],
      (indexdim[i,2]+1):indexdim[i+1,2]] <- x[[i]] 
  M 
}

PORF <- function(par){
  par[1:2] <- 1/(1+exp(-par[1:2]))
  par[3] <- exp(par[3])
  p3 <- par[3]-1
  if (par[3]!=1) {
    a <- 1+sum(par[-3])*p3
    b <- -4*p3*prod(par)
    p11 <- (a-sqrt(a*a+b))/(2*p3)
  }  else p11 <- prod(par[-3])
  p11
}

MPORF <- function(a,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type,...){
  #        a: vector of parameters etas with the first eta[1]=0
  #    ncat1: number of categories of the first response
  #    ncat2: number of categories of the second response
  #  ncat1_1: ncat1-1
  #  ncat2_1: ncat2-1
  #     ncat: ncat1*ncat2
  #   ncat_1: ncat-1
  #   ncat12: ncat1 + ncat2
  # ncat12_1: ncat12-1
  # ncat_1_2: ncat1_1*ncat2_1
  if (type=="ss") a[,2:(ncat1_1+ncat2_1+1)] <- -1*a[,2:(ncat1_1+ncat2_1+1)]
  P <- A <- matrix(0,ncat1,ncat2)
  k <- 0
  for (i in 1:(ncat1_1)){
    for (j in 1:(ncat2_1)){
      k <- k + 1
      i1 <- i + 1
      j1 <- j + 1
      A[i1,j1] <- PORF(a[c(i1, ncat1+j, ncat12_1 + k)])
      P[i,j] <- A[i,j] - A[i,j1] - A[i1,j] + A[i1,j1]
    }
  }
  P[-ncat1,ncat2]<-1/(1+exp(-a[2:ncat1])) - if (ncat1==2) P[-ncat1,-ncat2] else 
                     (rowSums(P[-ncat1,-ncat2])+c(0,1/(1+exp(-a[2:(ncat1_1)]))))
  P[ncat1,-ncat2]<-1/(1+exp(-a[(ncat1+1):(ncat12_1)]))-if (ncat2==2) P[-ncat1,-ncat2] else 
                     (colSums(P[-ncat1,-ncat2])+c(0,1/(1+exp(-a[(ncat1+1):(ncat12-2)]))))

  onesP <- 1 - sum(P) 
  P[ncat1,ncat2] <- pmax(0,onesP)  
  as.vector(t(P))
}



mPORF2 <- function(par){
  par[, 1:2] <- 1 / (1 + exp( - par[, 1:2]))
  par[, 3] <- exp(par[, 3])
  correl <- par[, 3] - 1
  or <- (correl != 0)
  p11 <- vector(, nrow(par))
  a <- 1 + rowSums(par[or, - 3]) * correl[or]
  b <- - 4 * correl[or] * par[or, 1] * par[or, 2] * par[or, 3] 
  p11[or] <- (a - sqrt( a * a + b)) / (2 * correl[or])
  if (any(!or)) p11[!or] <- par[!or,1] * par[!or,2]
  p11
}


mMPORF2 <- function(a,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type,...){
  if (type=="ss") a[,2:(ncat1_1+ncat2_1+1)] <- -1*a[,2:(ncat1_1+ncat2_1+1)]
  P <- A <- array(0,c(ncat1,ncat2,nrow(a)))
  k <- 0
  for (i in 1:(ncat1_1)){
    for (j in 1:(ncat2_1)){
      k <- k + 1
      i1 <- i + 1
      j1 <- j + 1
      A[i1,j1,] <- mPORF2(a[,c(i1, ncat1+j, ncat12_1 + k)])
      P[i,j,] <- A[i,j,] - A[i,j1,] - A[i1,j,] + A[i1,j1,]
    }
  }
  
  P[-ncat1,ncat2,] <- t(1/(1+exp(-a[,2:ncat1]))) - if (ncat1==2) P[-ncat1,-ncat2,] else 
    (apply( P[-ncat1,-ncat2,], 3, rowSums) + t(cbind(0,1/(1+exp(-a[,2:(ncat1_1)]))) ))
  
  P[ncat1,-ncat2,] <- t(1/(1+exp(-a[,(ncat1+1):(ncat12_1)]))) - if (ncat2==2) P[-ncat1,-ncat2,] else 
    (apply( P[-ncat1,-ncat2,],3,colSums) + rbind(0,t(1/(1+exp(-a[,(ncat1+1):(ncat12-2)])))))
  
  onesP <- 1 - apply(P,3,sum) 
  P[ncat1,ncat2,] <- pmax(0,onesP)
  matrix(aperm(P,c(3,2,1)),nrow(a),ncat)
}


dplackett.pblm <- function(a,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type,...){
  if (nrow(a)>1) 
    mMPORF2(a,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type) else
    MPORF(a,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type)
}



find.pi3 <- function(neweta,new.nu,M,C,Ct,acc2=1e-06,maxit2=100,tol2){
  C_Cdqrls<-getNativeSymbolInfo('Cdqrls', PACKAGE=getLoadedDLLs()$stats)
  iter2 <- 0
  while (any(tol2>acc2)&&(iter2<maxit2)){
    iter2 <- iter2+1
    nu <- new.nu
    enu <- exp(nu)
    Menu <- drop(M%*%enu)
    der1 <- crossprod(C/Menu,M)
    H <- -der1%*%diag(enu)
    U <- H%*%nu + Ct%*%log(Menu)-neweta
    new.nu <- .Call(C_Cdqrls,H,U, 1e-11, FALSE)$coefficients
    tol2 <- abs(new.nu-nu)
  }
  p <- exp(new.nu)
  attr(p,"der1") <- der1
  return(p)
}

find.score <- function(C,M,prob) crossprod((C/as.vector(M%*%prob)),M)

find.score2 <- function(C,M,prob){
 Nc <- ncol(C)
 Nr <- nrow(C)
 Nrp <- if (is.null(nrow(prob))) 1 else nrow(prob)
 Mp <- tcrossprod(prob,M) 
 GRAD0 <- aperm(array(apply(Mp,1,function(x)C/x),c(Nr,Nc,Nrp)),c(2,1,3))
 GRAD <- array(apply(GRAD0,3,function(x)x%*%M),c(Nc,Nc,Nrp))
 attributes(GRAD)$"Mp" <- Mp
 GRAD
}

make.marginals <- function(ncat1,ncat2,type="gg"){
  totalsum <- matrix(1,1,ncat1*ncat2)
  marginal.logit1 <- function(i,ncat1,ncat2,type){
    switch(substr(type,1,1),
           "g" = rbind(kronecker(((1:ncat1)<=i),rep(1,ncat2)),kronecker(((1:ncat1)>i),rep(1,ncat2))),
           "c" = rbind(kronecker(((1:ncat1)<=i),rep(1,ncat2)),kronecker(((1:ncat1)==(i+1)),rep(1,ncat2))),
           "l" = rbind(kronecker(((1:ncat1)==i),rep(1,ncat2)),kronecker(((1:ncat1)==(i+1)),rep(1,ncat2))),         
           "r" = rbind(kronecker(((1:ncat1)==i),rep(1,ncat2)),kronecker(((1:ncat1)>i),rep(1,ncat2))),
           "b" = rbind(kronecker(((1:ncat1)==1),rep(1,ncat2)),kronecker(((1:ncat1)==(i+1)),rep(1,ncat2))),
           "s" = rbind(kronecker(((1:ncat1)==ncat1),rep(1,ncat2)),kronecker(((1:ncat1)==i),rep(1,ncat2))))  
  }
  marginal.logit2 <- function(i,ncat1,ncat2,type){
    switch(substr(type,2,2),
           "g" = rbind(rep(((1:ncat1)<=i)*1,ncat2),rep(((1:ncat1)>i)*1,ncat2)),
           "c" = rbind(rep(((1:ncat1)<=i)*1,ncat2),rep(((1:ncat1)==(i+1))*1,ncat2)),
           "r" = rbind(rep(((1:ncat1)==i)*1,ncat2),rep(((1:ncat1)>i)*1,ncat2)),         
           "l" = rbind(rep(((1:ncat1)==i)*1,ncat2),rep(((1:ncat1)==(i+1))*1,ncat2)),
           "b" = rbind(rep(((1:ncat1)==1)*1,ncat2),rep(((1:ncat1)==(i+1))*1,ncat2)),
           "s" = rbind(rep(((1:ncat1)==ncat1)*1,ncat2),rep(((1:ncat1)==i)*1,ncat2)))
  }
  Mrow <- Mcol <-NULL
  for (i in 1:(ncat1-1)) Mrow <- rbind(Mrow,marginal.logit1(i,ncat1,ncat2,type))
  for (i in 1:(ncat2-1)) Mcol <- rbind(Mcol,marginal.logit2(i,ncat2,ncat1,type))
  M <- as.matrix(rbind(totalsum,Mrow,Mcol))
  k<-0
  P<-Q<-R<-S<-matrix(0,nrow=ncat1,ncol=ncat2)
  P1<-Q1<-R1<-S1<-matrix(0,nrow=(ncat1-1)*(ncat2-1),ncol=ncat1*ncat2)
  
  marginal.ORs <- function(r,s,u,v,type,ncat1,ncat2){
    switch(type,
           "gg" = c((u<=r)&(v<=s),(u<=r)&(v>s),(u>r)&(v<=s),(u>r)&(v>s)),
           "gl" = c((u<=r)&(v==s),(u<=r)&(v==(s+1)),(u>r)&(v==s),(u>r)&(v==(s+1))),
           "gc" = c((u<=r)&(v==s),(u<=r)&(v>s),(u>r)&(v==s),(u>r)&(v>s)),
           "gr" = c((u<=r)&(v<=s),(u<=r)&(v==(s+1)),(u>r)&(v<=s),(u>r)&(v==(s+1))),              
           "gb" = c((u<=r)&(v==1),(u<=r)&(v==(s+1)),(u>r)&(v==1),(u>r)&(v==(s+1))),    
           "gs" = c((u<=r)&(v==ncat2),(u<=r)&(v==s),(u>r)&(v==ncat2),(u>r)&(v==s)),    
           
           "ll" = c((u==r)&(v==s),(u==r)&(v==(s+1)),(u==(r+1))&(v==s),(u==(r+1))&(v==(s+1))),
           "lg" = c((u==r)&(v<=s),(u==r)&(v>s),(u==(r+1))&(v<=s),(u==(r+1))&(v>s)),
           "lc" = c((u==r)&(v==s),(u==r)&(v>s),(u==(r+1))&(v==s),(u==(r+1))&(v>s)),
           "lr" = c((u==r)&(v<=s),(u==r)&(v==(s+1)),(u==(r+1))&(v<=s),(u==(r+1))&(v==(s+1))),
           "lb" = c((u==r)&(v==1),(u==r)&(v==(s+1)),(u==(r+1))&(v==1),(u==(r+1))&(v==(s+1))),
           "ls" = c((u==r)&(v==ncat2),(u==r)&(v==s),(u==(r+1))&(v==ncat2),(u==(r+1))&(v==s)),
           
           "cc" = c((u==r)&(v==s),(u==r)&(v>s),(u>r)&(v==s),(u>r)&(v>s)),
           "cg" = c((u==r)&(v<=s),(u==r)&(v>s),(u>r)&(v<=s),(u>r)&(v>s)),
           "cl" = c((u==r)&(v==s),(u==r)&(v==(s+1)),(u>r)&(v==s),(u>r)&(v==(s+1))),
           "cr" = c((u==r)&(v<=s),(u==r)&(v==(s+1)),(u>r)&(v<=s),(u>r)&(v==(s+1))),
           "cb" = c((u==r)&(v==1),(u==r)&(v==(s+1)),(u>r)&(v==1),(u>r)&(v==(s+1))),
           "cs" = c((u==r)&(v==ncat2),(u==r)&(v==s),(u>r)&(v==ncat2),(u>r)&(v==s)),
           
           "rr" = c((u<=r)&(v<=s),(u<=r)&(v==(s+1)),(u==(r+1))&(v<=s),(u==(r+1))&(v==(s+1))),
           "rg" = c((u<=r)&(v<=s),(u<=r)&(v>s),(u==(r+1))&(v<=s),(u==(r+1))&(v>s)),
           "rl" = c((u<=r)&(v==s),(u<=r)&(v==(s+1)),(u==(r+1))&(v==s),(u==(r+1))&(v==(s+1))),
           "rc" = c((u<=r)&(v==s),(u<=r)&(v>s),(u==(r+1))&(v==s),(u==(r+1))&(v>s)),
           "rb" = c((u<=r)&(v==1),(u<=r)&(v==(s+1)),(u==(r+1))&(v==1),(u==(r+1))&(v==(s+1))),
           "rs" = c((u<=r)&(v==ncat2),(u<=r)&(v==s),(u==(r+1))&(v==ncat2),(u==(r+1))&(v==s)),
           
           "bb" = c((u==1)&(v==1),(u==1)&(v==(s+1)),(u==(r+1))&(v==1),(u==(r+1))&(v==(s+1))),
           "br" = c((u==1)&(v<=s),(u==1)&(v==(s+1)),(u==(r+1))&(v<=s),(u==(r+1))&(v==(s+1))),    
           "bl" = c((u==1)&(v==s),(u==1)&(v==(s+1)),(u==(r+1))&(v==s),(u==(r+1))&(v==(s+1))), 
           "bc" = c((u==1)&(v==s),  (u==1)&(v>s),   (u==(r+1))&(v==s),(u==(r+1))&(v>s)),
           "bg" = c((u==1)&(v<=s),  (u==1)&(v>s),   (u==(r+1))&(v<=s),(u==(r+1))&(v>s)),
           "bs" = c((u==1)&(v==ncat2),(u==1)&(v==s),(u==(r+1))&(v==ncat2),(u==(r+1))&(v==s)),
           
           "ss" = c((u==ncat1)&(v==ncat2),(u==ncat1)&(v==s),(u==r)&(v==ncat2),(u==r)&(v==s)),    
           "sb" = c((u==ncat1)&(v==1),(u==ncat1)&(v==(s+1)),(u==r)&(v==1),(u==r)&(v==(s+1))),
           "sr" = c((u==ncat1)&(v<=s),(u==ncat1)&(v==(s+1)),(u==r)&(v<=s),(u==r)&(v==(s+1))),    
           "sl" = c((u==ncat1)&(v==s),(u==ncat1)&(v==(s+1)),(u==r)&(v==s),(u==r)&(v==(s+1))), 
           "sc" = c((u==ncat1)&(v==s),  (u==ncat1)&(v>s),   (u==r)&(v==s),(u==r)&(v>s)),
           "sg" = c((u==ncat1)&(v<=s),  (u==ncat1)&(v>s),   (u==r)&(v<=s),(u==r)&(v>s))     )
  }          
  for (r in 1:(ncat1-1)){
    for (s in 1:(ncat2-1)){
      for (u in 1:ncat1){
        for (v in 1:ncat2){
          mm <- 1*marginal.ORs(r,s,u,v,type,ncat1,ncat2)
          P[u,v] <- mm[1]
          Q[u,v] <- mm[2]
          R[u,v] <- mm[3]
          S[u,v] <- mm[4]   
        }
      }   
      k<-k+1
      P1[k,] <- as.vector(t(P))
      Q1[k,] <- as.vector(t(Q))
      R1[k,] <- as.vector(t(R))
      S1[k,] <- as.vector(t(S))
    }
  }
  rbind(M,P1,Q1,R1,S1)
} 

make.contrasts <- function(ncat1,ncat2){
  Mtot <- matrix(1,1,1)
  Crow <- kronecker(diag(ncat1-1),c(1,-1))
  Ccol <- kronecker(diag(ncat2-1),c(1,-1)) 
  Id <- diag(1,(ncat1-1)*(ncat2-1))
  Mdiag2(list(Mtot,Crow,Ccol,rbind(Id,-Id,-Id,Id))) 
}



prop <-function(x,xx,nc,z){
  propo[[1]] <- matrix(0,nrow=nc,ncol=sum(xx))
  cxx <- cumsum(xx)
  aa <- 0
  for (i in 1:length(x)){
    propo[[1]][,(1+aa):cxx[i]] <- if (x[i]) t(rep(z[i],nc)) else diag(z[i],nc)  
    aa<-cxx[i] 
  } 
  propo[[1]]
}  


make.RC <- function(RC.fo,ncat1_1,ncat2_1,contrasts){
  dd <- expand.grid("Col"=1:ncat2_1,"Row"=1:ncat1_1)
  dd$Row <- as.factor(dd$Row)
  dd$Col <- as.factor(dd$Col)
  model.matrix(RC.fo,dd,contrasts=contrasts) 
}

prop.RC <-function(x,xx,nc,z,RC.fo,ncat1,ncat2,contrasts){
  if (is.null(RC.fo)) propo[[1]] <- prop(x,xx,nc,z) else  {
    RC <- make.RC(RC.fo,ncat1-1,ncat2-1,contrasts)
    if (length(x)>1){
      xx <- xx[-1]
      z <- z[-1]
      x <- x[-1]
      propo[[1]] <- matrix(0,nrow=nc,ncol=sum(xx))
      aa <- 0
      cxx <- cumsum(xx)
      nom <- c()
      for (i in 1:length(x)){
        propo[[1]][,(1+aa):cxx[i]] <- if (x[i]) t(rep(z[i],nc)) else diag(z[i],nc) 
        nom<-c(nom,names(z[i]))
        aa<-cxx[i] 
      } 
      colnames(propo[[1]])<-nom
    } else propo[[1]] <- NULL  
    propo[[1]]<-cbind(RC,propo[[1]])
  }
  propo[[1]]
}  



penalty.INEQUALITY <- function(lamv1,lamv2,ncat1,ncat2,eta,Z,D3,m2,type){
  if ((ncat1==2)&(ncat2==2)) matrix(0,ncol(Z),ncol(Z)) else 
  {            
    lamv1 <- ifelse(!(substr(type,1,1)=="g"),0,lamv1)
    lamv2 <- ifelse(!(substr(type,1,1)=="g"),0,lamv2)
    lam <- rep(c(rep(lamv1,ncat1-2),rep(lamv2,ncat2-2)),m2)
    DE <- D3%*%as.vector(t(eta))
    II <- crossprod(sqrt(lam)*I(DE<=0),(1+D3%*%Z))
    crossprod(II) 
  }   
}

DIFF<- function(ncat1,ncat2,m2){
  bin <- (ncat1==2)&(ncat2==2)
  D1q <- if (ncat1==2) NULL else diff(diag(1,ncat1-1))
  D2q <- if (ncat2==2) NULL else diff(diag(1,ncat2-1))
  D12 <- if (bin) NULL else Mdiag(list(D1q,D2q))
  if (bin) NULL else kronecker(diag(1,m2),D12)
}


penalization <- function(nc,d,pnl.type){
  switch(pnl.type,   
         "ARC1"  =  diff(diag(1,nc),d=d),
         "ARC2"  =  diff(diag(1,nc),d=d),
         "ridge" =  diag(1,nc),   
         "equal" =  diag(1,nc),
         "lasso" =  diag(1,nc),
         "lassoV" = diff(diag(1,nc),d=d)
  )                                      
}     

penalty.TERM <- function(lam1,lam2,lam3,lam4,propo,xx1,xx2,xx12,
                         totcol, ncat1,ncat2,d1,d2,d3,d4,pnl.type,new.par,iter){
  if (is.null(lam1)&is.null(lam2)&is.null(lam3)&is.null(lam4)) 
    matrix(0,totcol,totcol) else {       
      switch(pnl.type,
             "ARC1" = penalty.ARC1(lam1,lam2,lam3,propo,xx1,xx2,xx12,
                                   totcol,ncat1,ncat2,d1,d2,d3,pnl.type),  
             "ARC2" = penalty.ARC2(lam1,lam2,lam3,lam4,propo,xx1,xx2,xx12,
                                   totcol,ncat1,ncat2,d1,d2,d3,d4,pnl.type),                    
             "ridge" = penalty.RIDGE(lam1,lam2,lam3,propo, xx1,xx2,xx12,
                                     totcol,ncat1,ncat2,pnl.type),
             "lasso" = penalty.LASSO(lam1,lam2,lam3,propo, xx1,xx2,xx12,
                                     totcol, ncat1,ncat2,pnl.type,new.par,iter),
             "lassoV" = penalty.VLASSO(lam1,lam2,lam3,propo, xx1,xx2,xx12,
                                     totcol, ncat1,ncat2,pnl.type,new.par,iter),
             "equal" = penalty.equality(lam1,xx1,ncat1,d1,propo[[1]],pnl.type,xx12))                                          
    }      
}  

# penalty.LASSO<- function(lam1,lam2,lam3,propo, xx1,xx2,xx12,
#                          totcol, ncat1,ncat2,pnl.type,new.par,iter){
#   
#   if (all(is.null(lam1),is.null(lam2),is.null(lam3))) 
#     matrix(0,totcol,totcol) else {    
#       if (is.null(lam1))  {
#         P1 <- matrix(0,sum(xx1),sum(xx1)) 
#         D1q<-matrix(0,ncat1-2,ncat1-1)
#       } else {
#         p1 <- propo[[1]][-1]
#         D1 <- NULL
#         D1q <- penalization((ncat1-1),NULL,pnl.type)
#         index<-rep(1,ncat1-1)   
#         if ((length(xx1)>1)&(sum(!p1))) {
#           for (i in 2:length(xx1)) index<-c(index,rep(propo[[1]][i],xx1[i]))
#           npar <- new.par[1:sum(xx1)]  
#           D1 <- kronecker(diag(sqrt(lam1[-1]),sum(p1==0)),D1q)
#           D1 <- D1*ginv(diag(sqrt(abs(npar[!index])),sum(!index)))
#         }  
#         D1q <- D1q*ginv(diag(sqrt(abs(npar[1:(ncat1-1)])),ncat1-1))
#         
#         INT <- if (!propo[[1]][1]) lam1[1]*crossprod(D1q) else as.matrix(0)  
#         PEN_PROP <- if (any(p1)) matrix(0,sum(p1),sum(p1)) else NULL  
#         PEN_NOT_PROP <- if (any(!p1)) crossprod(D1) else NULL
#         P1 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))  
#         
#       }
#      
#       if (is.null(lam2))  {
#         P2 <- matrix(0,sum(xx2),sum(xx2)) 
#         D2q <- matrix(0,ncat2-2,ncat2-1)
#       } else {
#         p2 <- propo[[2]][-1]
#         D2 <- NULL
#         D2q <- penalization((ncat2-1),NULL,pnl.type)
#         index <- rep(TRUE,ncat2-1)            
#         npar <- new.par[(1+sum(xx1)):sum(xx1,xx2)]
#         if ((length(xx2)>1)&(sum(!p2))) {
#           for (i in 2:length(xx2)) index<-c(index,rep(propo[[2]][i],xx2[i]))
#           D2 <- kronecker(diag(sqrt(lam2[-1]),sum(!p2)),D2q)
#           D2 <- D2*ginv(diag(sqrt(abs(npar[!index])),sum(!index) ))
#         }  
#         D2q <- D2q*ginv(diag(sqrt(abs(npar[1:(ncat2-1)])),ncat2-1))
#         
#         INT <- if (!propo[[2]][1]) lam2[1]*crossprod(D2q) else as.matrix(0)
#         PEN_PROP <- if (any(p2)) matrix(0,sum(p2),sum(p2)) else NULL  
#         PEN_NOT_PROP <- if (any(!p2)) crossprod(D2) else NULL
#         P2 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))      
#       }
# 
#       if (is.null(lam3)) P12 <- matrix(0,sum(xx12),sum(xx12)) else {
#         p3 <- propo[[3]][-1]
#         D12 <- NULL
#         ncat12_1 <- (ncat1-1)*(ncat2-1)
#         D12q <- penalization(ncat12_1,NULL,pnl.type)     
#         npar <- new.par[(1+sum(xx1,xx2)):length(new.par)]   
#         index <- rep(1,((ncat1-1)*(ncat2-1)))
#         D12q <- D12q*ginv(diag(sqrt(abs(npar[1:ncat12_1])),ncat12_1))                    
#         if ((length(xx12)>1)&(sum(!p3)))  {
#           for (i in 2:length(xx12)) index<-c(index,rep(propo[[3]][i],xx12[i]))            
#           D12 <- kronecker(diag(sqrt(lam3[-1]),sum(!p3)),
#                            D12q)*ginv(diag(sqrt(abs(npar[!index])),sum(!index)))
#         }                
#         D12qr<-kronecker(diag(1,(ncat1-1)),D1q)
#         D12qc<-kronecker(D2q,diag(1,(ncat2-1)))
#         INT <- if (!propo[[3]][1])  lam3[1]*crossprod(D12q) else as.matrix(0)
#         PEN_PROP <- if (any(p3)) matrix(0,sum(p3),sum(p3)) else NULL 
#         PEN_NOT_PROP <- if (any(!p3)) crossprod(D12)  else NULL
#         P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))        
#       }              
#       Mdiag(list(P1,P2,P12))
#     }  
# }

penalty.LASSO <- function(lam1,lam2,lam3,propo, xx1,xx2,xx12,
                         totcol, ncat1,ncat2,pnl.type,new.par,iter){
  
  if (all(is.null(lam1),is.null(lam2),is.null(lam3))) 
    matrix(0,totcol,totcol) else {    
      if (is.null(lam1))  {
        P1 <- matrix(0,sum(xx1),sum(xx1)) 
        D1q<-matrix(0,ncat1-2,ncat1-1)
      } else {
        p1 <- propo[[1]][-1]
        D1 <- NULL
        D1q <- penalization((ncat1-1),NULL,pnl.type)
        index<-rep(1,ncat1-1)   
        if ((length(xx1)>1)&(sum(!p1))) {
          for (i in 2:length(xx1)) index<-c(index,rep(propo[[1]][i],xx1[i]))
          npar <- new.par[1:sum(xx1)]  
          D1 <- kronecker(diag(lam1[-1],sum(p1==0)),D1q)
          #browser()
         # D1 <- D1*ginv(diag(sqrt(abs(npar[!index])),sum(!index)))
          D1 <- D1*ginv(diag(abs(npar[!index]),sum(!index)))
        }  
        D1q <- D1q*ginv(diag(abs(npar[1:(ncat1-1)]),ncat1-1))
        
        INT <- if (!propo[[1]][1]) lam1[1]*D1q else as.matrix(0)  
        PEN_PROP <- if (any(p1)) matrix(0,sum(p1),sum(p1)) else NULL  
        PEN_NOT_PROP <- if (any(!p1)) D1 else NULL
        P1 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))  
        
      }
      
      if (is.null(lam2))  {
        P2 <- matrix(0,sum(xx2),sum(xx2)) 
        D2q <- matrix(0,ncat2-2,ncat2-1)
      } else {
        p2 <- propo[[2]][-1]
        D2 <- NULL
        D2q <- penalization((ncat2-1),NULL,pnl.type)
        index <- rep(TRUE,ncat2-1)            
        npar <- new.par[(1+sum(xx1)):sum(xx1,xx2)]
        if ((length(xx2)>1)&(sum(!p2))) {
          for (i in 2:length(xx2)) index<-c(index,rep(propo[[2]][i],xx2[i]))
          D2 <- kronecker(diag(sqrt(lam2[-1]),sum(!p2)),D2q)
          D2 <- D2*ginv(diag(abs(npar[!index]),sum(!index) ))
        }  
        D2q <- D2q*ginv(diag(abs(npar[1:(ncat2-1)]),ncat2-1))
        
        INT <- if (!propo[[2]][1]) lam2[1]*D2q else as.matrix(0)
        PEN_PROP <- if (any(p2)) matrix(0,sum(p2),sum(p2)) else NULL  
        PEN_NOT_PROP <- if (any(!p2)) D2 else NULL
        P2 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))      
      }
      
      if (is.null(lam3)) P12 <- matrix(0,sum(xx12),sum(xx12)) else {
        p3 <- propo[[3]][-1]
        D12 <- NULL
        ncat12_1 <- (ncat1-1)*(ncat2-1)
        D12q <- penalization(ncat12_1,NULL,pnl.type)     
        npar <- new.par[(1+sum(xx1,xx2)):length(new.par)]   
        index <- rep(1,((ncat1-1)*(ncat2-1)))
        D12q <- D12q*ginv(diag(abs(npar[1:ncat12_1]),ncat12_1))                    
        if ((length(xx12)>1)&(sum(!p3)))  {
          for (i in 2:length(xx12)) index<-c(index,rep(propo[[3]][i],xx12[i]))    
         # browser()
          D12 <- kronecker(diag(lam3[-1],sum(!p3)),
                           D12q)*ginv(diag(abs(npar[!index]),sum(!index)))
        }                
        #D12qr<-kronecker(diag(1,(ncat1-1)),D1q)
        #D12qc<-kronecker(D2q,diag(1,(ncat2-1)))
        INT <- if (!propo[[3]][1])  lam3[1]*D12q else as.matrix(0)
        PEN_PROP <- if (any(p3)) matrix(0,sum(p3),sum(p3)) else NULL 
        PEN_NOT_PROP <- if (any(!p3)) D12  else NULL
        P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))        
      }              
      Mdiag(list(P1,P2,P12))
    }  
}

penalty.RIDGE <- function(lam1,lam2,lam3,propo,xx1,xx2,xx12,totcol,
                          ncat1,ncat2,pnl.type){
  if (is.null(lam1)) P1 <- matrix(0,sum(xx1),sum(xx1)) else {
    if (length(lam1)!=length(propo[[1]]))
      stop("length of lam1 does not match with length(propo[[1]])")
    
    D1q <- penalization((ncat1-1),NULL,pnl.type)      
    lamx1 <-c()
    for (i in 1:length(propo[[1]])) lamx1<-c(lamx1,if (propo[[1]][i]==0) rep(lam1[i],ncat1-1) else lam1[i])      
    P1 <- lamx1*Mdiag(lapply(propo[[1]],function(x) if (x==0) D1q else as.matrix(1)))
  }
  if (is.null(lam2)) P2 <- matrix(0,sum(xx2),sum(xx2)) else {
    if (length(lam2)!=length(propo[[2]]))
      stop("length of lam2 does not match with length(propo[[2]])")
    
    D2q <- penalization((ncat2-1),NULL,pnl.type)      
    lamx2 <-c()
    for (i in 1:length(propo[[2]])) lamx2<-c(lamx2,if (propo[[2]][i]==0) rep(lam2[i],ncat2-1) else lam2[i])      
    P2 <- lamx2*Mdiag(lapply(propo[[2]],function(x) if (x==0) D2q else as.matrix(1)))
  }
  sumxx12 <- sum(xx12)
  if (is.null(lam3)) P12 <- matrix(0,sumxx12,sumxx12) else {
    if (length(lam3)!=length(propo[[3]]))
      stop("length of lam3 does not match with length(propo[[3]])")
    
    D12q <- penalization((ncat1-1)*(ncat2-1),NULL,pnl.type)
    lamx3 <-c()
    for (i in 1:length(propo[[3]])) lamx3<-c(lamx3,if (propo[[3]][i]==0) rep(lam3[i],(ncat1-1)*(ncat2-1)) else lam3[i])      
    P12 <- lamx3*Mdiag(lapply(propo[[3]],function(x) if (x==0) D12q else as.matrix(1)))
  }
  Mdiag2(list(P1,P2,P12)) 
}


penalty.VLASSO<- function(lam1,lam2,lam3,propo, xx1,xx2,xx12,
                         totcol, ncat1,ncat2,pnl.type,new.par,iter){
 
  if (all(is.null(lam1),is.null(lam2),is.null(lam3))) 
    matrix(0,totcol,totcol) else {    
      if (is.null(lam1))  {
        P1 <- matrix(0,sum(xx1),sum(xx1)) 
        D1q<-matrix(0,ncat1-2,ncat1-1)
      } else {
        p1 <- propo[[1]][-1]
        D1 <- NULL
        D1q <- penalization((ncat1-1),1,pnl.type)
        npar <- new.par[1:sum(xx1)]  
        INV10 <- diag(ginv(diag(abs(npar[1:(ncat1-1)]),ncat1-1)))
        D1q <- crossprod(D1q)*INV10
        index<-rep(1,ncat1-1)   
        if ((length(xx1)>1)&(sum(!p1))) {
          for (i in 2:length(xx1)) index<-c(index,rep(propo[[1]][i],xx1[i]))
          INV1 <- diag(ginv(diag(abs(npar[!index]),sum(!index) )))
          D1 <- INV1*crossprod(kronecker(diag(lam1[-1],(sum(!p1))),penalization(ncat1-1,1,pnl.type)))
        }  
        INT <- if (!propo[[1]][1]) lam1[1]*D1q else as.matrix(0)  
        PEN_PROP <- if (any(p1)) matrix(0,sum(p1),sum(p1)) else NULL  
        PEN_NOT_PROP <- if (any(!p1)) D1 else NULL
        P1 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))  
        
      }
      #if (iter>1) browser()
      
      if (is.null(lam2))  {
        P2 <- matrix(0,sum(xx2),sum(xx2)) 
        D2q <- matrix(0,ncat2-2,ncat2-1)
      } else {
        p2 <- propo[[2]][-1]
        D2 <- NULL
        npar <- new.par[(1+sum(xx1)):sum(xx1,xx2)]
        D2q <- penalization((ncat2-1),1,pnl.type)
        INV20 <- diag(ginv(diag(abs(npar[1:(ncat2-1)]),ncat2-1)))
        D2q <- crossprod(D2q)*INV20
        index <- rep(TRUE,ncat2-1)            
        if ((length(xx2)>1)&(sum(!p2))) {
          for (i in 2:length(xx2)) index<-c(index,rep(propo[[2]][i],xx2[i]))
          #browser()
          INV2 <- diag(ginv(diag(abs(npar[!index]),sum(!index) )))
          #D2 <- lam2[-1]*INV2*crossprod(D2q)
          D2 <- INV2*crossprod(kronecker(diag(lam2[-1],sum(!p2)),penalization(ncat2-1,1,pnl.type)))
          
        }  
        INT <- if (!propo[[2]][1]) lam2[1]*D2q else as.matrix(0)
        PEN_PROP <- if (any(p2)) matrix(0,sum(p2),sum(p2)) else NULL  
        PEN_NOT_PROP <- if (any(!p2)) D2 else NULL
        P2 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))      
      }
      
      if (is.null(lam3)) P12 <- matrix(0,sum(xx12),sum(xx12)) else {
        p3 <- propo[[3]][-1]
        D12 <- NULL
        ncat12_1 <- (ncat1-1)*(ncat2-1)
        D12q <- penalization(ncat12_1,1,pnl.type)     
        npar <- new.par[(1+sum(xx1,xx2)):length(new.par)]   
        index <- rep(1,ncat12_1)
        INV120 <- diag(ginv(diag(abs(npar[1:ncat12_1]),ncat12_1)))
        D12q <- INV120*crossprod(D12q)               
        if ((length(xx12)>1)&(sum(!p3)))  {
          for (i in 2:length(xx12)) index<-c(index,rep(propo[[3]][i],xx12[i]))    
          INV12 <- diag(ginv(diag(abs(npar[!index]),sum(!index))))
          #browser()
          D12 <- INV12*crossprod(kronecker(diag(lam3[-1],sum(!p3)),penalization(ncat12_1,1,pnl.type)))
        }                
        INT <- if (!propo[[3]][1])  lam3[1]*D12q else as.matrix(0)
        PEN_PROP <- if (any(p3)) matrix(0,sum(p3),sum(p3)) else NULL 
        PEN_NOT_PROP <- if (any(!p3)) D12  else NULL
        P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))        
      }      
      #browser()
      Mdiag(list(P1,P2,P12))
    }  
}







penalty.marginal.ARC <- function(lam,xx,nc,d,prop,pnl.type){
  if (is.null(lam)|(all(lam==0))) P <- matrix(0,sum(xx),sum(xx)) else {
    if (nc==2) P <- matrix(0,sum(xx),sum(xx))  else {
      Dq <- penalization((nc-1),d[1],pnl.type)
      if (sum(prop[-1]==0)>0){
        D <- list()
        for (i in 1:(sum(prop[-1]==0))) D[[i]] <- penalization((nc-1),d[i+1],pnl.type)
        D <- Mdiag(D)
      }          
      INT <- if (prop[1]==0) lam[1]*crossprod(Dq) else as.matrix(0)
      PEN_PROP <- if (any(prop[-1]==1)) matrix(0,sum(prop[-1]==1),sum(prop[-1]==1)) else NULL
      PEN_NOT_PROP <- if (any(prop[-1]==0)) lam[-1]*crossprod(D) else NULL
      P <- Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))
    }
  }
  P
}

penalty.ARC1 <- function(lam1,lam2,lam3,propo,xx1,xx2,xx12,totcol,
                         ncat1,ncat2,d1,d2,d3,pnl.type){
  P1 <- penalty.marginal.ARC(lam1,xx1,ncat1,d1,propo[[1]],pnl.type)
  P2 <- penalty.marginal.ARC(lam2,xx2,ncat2,d2,propo[[2]],pnl.type)  
  sumxx12 <- sum(xx12)
  if (is.null(lam3)|all(lam3==0)) P12 <- matrix(0,sumxx12,sumxx12) else {
    if ((ncat1==2)&(ncat2==2))  P12 <- matrix(0,sumxx12,sumxx12)  else {
      if (pnl.type=="ARC1"){
        D12q <- penalization((ncat1-1)*(ncat2-1),d3[1],pnl.type)
        D12 <- list()
        for (i in 1:length(d3)) D12[[i]] <- penalization((ncat1-1)*(ncat2-1),d3[i],pnl.type)
        D12 <- Mdiag(D12)        
        if (any(propo[[3]][-1]==0)) {
          non_prop <- !propo[[3]][-1][which(propo[[3]][-1]==0)] 
          D12 <- if (propo[[3]][1]==0) kronecker(diag(sqrt(lam3[-1][propo[[3]][-1]==0]*non_prop),sum(propo[[3]][-1]==0)),D12q) else 
            kronecker(diag(sqrt(lam3[-1]*non_prop),sum(propo[[3]][-1]==0)),D12q)
        } else D12 <- 0*D12
      } 
      INT <- if (propo[[3]][1]==0) lam3[1]*crossprod(D12q) else as.matrix(0)
      PEN_PROP <- if (any(propo[[3]][-1]==1)) matrix(0,sum(propo[[3]][-1]==1),sum(propo[[3]][-1]==1)) else NULL
#browser()
      PEN_NOT_PROP <- if (any(propo[[3]][-1]==0)) crossprod(D12) else NULL
      P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))
    } 
  }    
  Mdiag(list(P1,P2,P12))
} 



penalty.ARC2 <- function(lam1,lam2,lam3,lam4,propo,xx1,xx2,xx12,
                         totcol,ncat1,ncat2,d1,d2,d3,d4,pnl.type){
  P1 <- penalty.marginal.ARC(lam1,xx1,ncat1,d1,propo[[1]],pnl.type)
  P2 <- penalty.marginal.ARC(lam2,xx2,ncat2,d2,propo[[2]],pnl.type)  
  sumxx12 <- sum(xx12)
  if ((is.null(lam3)|all(lam3==0))&(is.null(lam4)|all(lam4==0)))
    P12 <- matrix(0,sumxx12 ,sumxx12) else {
      if ((ncat1==2)&(ncat2==2))  P12 <- matrix(0,sumxx12,sumxx12)  else {
        D12qr <- kronecker(diag(1,(ncat2-1)),penalization((ncat1-1),d3[1],pnl.type))
        D12qc <- kronecker( penalization((ncat2-1),d4[1],pnl.type),diag(1,(ncat1-1)))            
        if (sum(propo[[3]][-1]==0)>0) {
          D3 <- list()
          D4 <- list()
          for (i in 1:(sum(propo[[3]][-1]==0))) {            
            D3[[i]] <- kronecker(diag(1,(ncat2-1)),penalization((ncat1-1),d3[i+1],pnl.type))            
            D4[[i]] <- kronecker(penalization((ncat2-1),d4[i+1],pnl.type),diag(1,(ncat1-1)))
          }                         
          D12r <- Mdiag(D3)
          D12c <- Mdiag(D4)           
        } 
        INT <- if (propo[[3]][1]==0) lam3[1]*crossprod(D12qr)+lam4[1]*crossprod(D12qc) else as.matrix(0)
        PEN_PROP <- if (any(propo[[3]][-1]==1)) matrix(0,sum(propo[[3]][-1]==1),sum(propo[[3]][-1]==1)) else NULL
        PEN_NOT_PROP <- if (any(propo[[3]][-1]==0)) lam3[-1]*crossprod(D12r)+lam4[-1]*crossprod(D12c) else NULL
        P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))
      } 
    }    
  Mdiag(list(P1,P2,P12))
} 



# penalty.lassoV <- function(lam1,lam2,lam3,propo,xx1,xx2,xx12,totcol,
#                          ncat1,ncat2,d1,d2,d3,pnl.type){
#   P1 <- penalty.marginal.ARC(lam1,xx1,ncat1,d1,propo[[1]],pnl.type)
#   P2 <- penalty.marginal.ARC(lam2,xx2,ncat2,d2,propo[[2]],pnl.type)  
#   sumxx12 <- sum(xx12)
#   if (is.null(lam3)|all(lam3==0)) P12 <- matrix(0,sumxx12,sumxx12) else {
#     if ((ncat1==2)&(ncat2==2))  P12 <- matrix(0,sumxx12,sumxx12)  else {
#       if (pnl.type=="lassoV"){
#         D12q <- penalization((ncat1-1)*(ncat2-1),d3[1],pnl.type)
#         D12 <- list()
#         for (i in 1:length(d3)) D12[[i]] <- penalization((ncat1-1)*(ncat2-1),d3[i],pnl.type)
#         D12 <- Mdiag(D12)        
#         if (any(propo[[3]][-1]==0)) {
#           non_prop <- !propo[[3]][-1][which(propo[[3]][-1]==0)] 
#           D12 <- if (propo[[3]][1]==0) kronecker(diag(sqrt(lam3[propo[[3]][-1]==0]*non_prop),sum(propo[[3]][-1]==0)),D12q) else 
#             kronecker(diag(sqrt(lam3*non_prop),sum(propo[[3]][-1]==0)),D12q)
#         } else D12 <- 0*D12
#       } 
#       INT <- if (propo[[3]][1]==0) lam3[1]*crossprod(D12q) else as.matrix(0)
#       PEN_PROP <- if (any(propo[[3]][-1]==1)) matrix(0,sum(propo[[3]][-1]==1),sum(propo[[3]][-1]==1)) else NULL
#       PEN_NOT_PROP <- if (any(propo[[3]][-1]==0)) crossprod(D12) else NULL
#       P12 <-  Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))
#     } 
#   }    
#   Mdiag(list(P1,P2,P12))
# } 




penalty.equality <- function(lam,xx,nc,d,prop,pnl.type,xx12){
  lam <- sqrt(lam)
  if (is.null(lam)|(all(lam==0))) P1 <- matrix(0,sum(xx),sum(xx)) else {
    Dq <- penalization((nc-1),d[1],pnl.type)
    if (sum(prop[-1]==0)>0){
      D <- list()
      for (i in 1:(sum(prop[-1]==0))) D[[i]] <- penalization((nc-1),d[i+1],pnl.type)
      D <- Mdiag(D)
    }          
    INT <- if (prop[1]==0) lam[1]*crossprod(Dq) else as.matrix(0)
    PEN_PROP <- if (any(prop[-1]==1)) lam*diag(sum(prop[-1]==1)) else NULL
    PEN_NOT_PROP <- if (any(prop[-1]==0)) lam*crossprod(D) else NULL
    P1 <- Mdiag(list(INT,PEN_PROP,PEN_NOT_PROP))
  }
  P2 <- -P1
  sumxx12 <- sum(xx12)
  P12 <- matrix(0,sumxx12,sumxx12)
  Mdiag(list(crossprod(cbind(P1,P2)),P12))
} 



pbs <- function (x, df = 3, lambda = NULL, ps.intervals = 20, degree = 3, order = 3) 
{
  scall <- deparse(sys.call())
  if (is.null(lambda) & is.null(df)) 
    stop("df or lambda should be set \n")
  number.knots <- ps.intervals + 2 * degree + 1
  x.domain <- as.vector(x)
  xl <- min(x.domain)
  xr <- max(x.domain)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  dx <- (xmax - xmin)/ps.intervals
  nx <- names(x.domain)
  nax <- is.na(x.domain)
  if (nas <- any(nax)) 
    x.domain <- x[!nax]
  sorder <- degree + 1
  Aknots <- range(x.domain)
  nAknots <- ps.intervals - 1
  if (nAknots < 1) {
    nAknots <- 1
    warning(paste("ps.intervals was too small; have used 2"))
  }
  if (nAknots > 0) {
    Aknots <- seq(from = xmin - degree * dx, to = xmax + 
                    degree * dx, by = dx)
  } else knots <- NULL
  basis <- splineDesign(Aknots, x.domain, sorder, 0 * x.domain)
  n.col <- ncol(basis)
  dimnames(basis) <- list(1:nrow(basis), 1:n.col)
  if ((order - n.col + 1) > 0) {
    order <- n.col - 1
    warning(paste("order was too large; have used ", n.col - 1))
  }
  if (df < 1) warning("the df are set to 1")
  df <- if (df < 1) 1 else df + 2
  if (!is.null(lambda)) {
    if (lambda < 0) {
      lambda <- 0
      warning(paste("lambda was negative; have used ", lambda))
    }
  }
  aug <- diag(n.col)
  if (order != 0) {
    for (tt in 1:order) {
      aug <- diff(aug)
    }
  }
  pen.aug <- aug
  xvar <- x
  attr(xvar, "smooth.name") <- "pbs"
  attr(xvar, "knots") <- Aknots
  attr(xvar, "D") <- pen.aug
  attr(xvar, "X") <- basis
  attr(xvar, "call") <- substitute(gamlss.ps(data[[scall]], z, w))
  attr(xvar, "lambda") <- lambda
  attr(xvar, "df") <- df
  attr(xvar, "order") <- order
  xvar
}

pb <- function (x, df = NULL, lambda = NULL, control = pb.control(...), 
                ...){
  bbase <- function(x, xl, xr, ndx, deg, quantiles = FALSE) {
    tpower <- function(x, t, p) (x - t)^p * (x > t)
    dx <- (xr - xl)/ndx
    if (quantiles) {
      knots <- sort(c(seq(xl - deg * dx, xl, dx), 
                      quantile(x, prob = seq(0, 1, length = ndx)), 
                      seq(xr, xr + deg * dx, dx)))
      B <- splineDesign(knots, x = x, outer.ok = TRUE, ord = deg + 1)
      return(B)
    } else {
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      P <- outer(x, knots, tpower, deg)
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
      B <- (-1)^(deg + 1) * P %*% t(D)
      B
    }
  }
  scall <- deparse(sys.call())
  lx <- length(x)
  control$inter <- if (lx < 99) 10 else control$inter
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  X <- bbase(x, xmin, xmax, control$inter, control$degree, 
             control$quantiles)
  r <- ncol(X)
  D <- if (control$order == 0) diag(r) else 
          diff(diag(r), diff = control$order)
  if (!is.null(df)) {
    if (df > (dim(X)[2] - 2)) {
      df <- 3
      warning("The df's exceed the number of columns of the design matrix", 
              "\n", "   they are set to 3")
    }
    if (df < 0) 
      warning("the extra df's are set to 0")
    df <- if (df < 0) 2 else df + 2
  }
  xvar <- x
  attr(xvar, "smooth.name") <- "pb"
  attr(xvar, "control") <- control
  attr(xvar, "D") <- D
  attr(xvar, "X") <- X
  attr(xvar, "df") <- df
  attr(xvar, "lambda") <- lambda
  attr(xvar, "class") <- "smooth"
  xvar
}


pb.control <- function (inter = 20, degree = 3, order = 2, quantiles = FALSE, ...){
  if (inter <= 0) {
    warning("the value of inter supplied is less than 0, the value of 10 was used instead")
    inter <- 10
  }
  if (degree <= 0) {
    warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
    degree <- 3
  }
  if (order < 0) {
    warning("the value of order supplied is zero or negative the default value of 2 was used instead")
    order <- 2
  }
 
  list(inter = inter, degree = degree, order = order, 
       quantiles = as.logical(quantiles)[1])
}

expand.pblm <- function(x){
  x <- x[,ncol(x):1]
  lev1 <- levels(factor(x[,1]))
  lev2 <- levels(factor(x[,2]))
  grid <- expand.grid(lev1,lev2)
  a <- t(apply(grid,1,function(x)paste(x,collapse="")))
  b <- apply(x,1,function(xx)paste(xx,collapse=""))
  t(sapply(b,function(x)(x==a)*1))
}


pvc <- function (x, df = NULL, lambda = NULL, by = NULL, control = pvc.control(...), 
                 ...){
  bbase <- function(x, xl, xr, ndx, deg, quantiles = FALSE) {
    tpower <- function(x, t, p) (x - t)^p * (x > t)
    dx <- (xr - xl)/ndx
    if (quantiles) {
      knots <- sort(c(seq(xl - deg * dx, xl, dx), 
                      quantile(x, prob = seq(0, 1, length = ndx)), seq(xr, xr + deg * dx, dx)))
      B <- splineDesign(knots, x = x, outer.ok = TRUE, 
                        ord = deg + 1)
      return(B)
    }
    else {
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      P <- outer(x, knots, tpower, deg)
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * 
                                            dx^deg)
      B <- (-1)^(deg + 1) * P %*% t(D)
      B
    }
  }
  no.dist.val <- length(table(x))
  lx <- length(x)
  control$inter <- if (lx < 99) 
    10
  else control$inter
  control$inter <- if (no.dist.val <= control$inter) 
    no.dist.val else control$inter
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  X <- bbase(x, xmin, xmax, control$inter, control$degree, 
             control$quantiles)
  r <- ncol(X)
  D <- if (control$order == 0) 
    diag(r)
  else diff(diag(r), diff = control$order)
  if (is.null(by)) {
    by.var <- NULL
    xvar <- x
    if (!is.null(df)) {
      if (df[1] > (dim(X)[2] - 2)) {
        df <- 3
        warning("The df's exceed the number of columns of the design matrix", 
                "\n", "   they are set to 3")
      }
      df <- if (any(df < 1)) 
        1
      else df
      if (any(df < 1)) 
        warning("the df are set to 1")
    }
  }
  else {
    if (is(by, "factor")) {
      by.var <- by
      if (!is.null(df[1])) {
        if (length(df) != nlevels(by)) 
          df <- rep(df[1], nlevels(by))
        if (any(df > (dim(X)[2] - 2))) {
          df <- ifelse((df > (dim(X)[2] - 2)), 3, df)
          warning("The df's exceed the number of columns of the design matrix", 
                  "\n", "   they are set to 3")
        }
        df <- ifelse(df < 1, 1, df)
        if (any(df < 1)) 
          warning("the df are set to 1")
      }
    }
    else {
      by.var <- by - mean(by)
      if (!is.null(df)) {
        if (df > (dim(X)[2] - 2)) {
          df <- 3
          warning("The df's exceed the number of columns of the design matrix", 
                  "\n", "   they are set to 3")
        }
        df <- if (df[1] < 1) 
          1
        else df + 2
        if (df[1] < 1) 
          warning("the df are set to 1")
      }
    }
    xvar <- rep(0, length(x))
  }
  attr(xvar, "control") <- control
  attr(xvar, "D") <- D
  attr(xvar, "X") <- X
  attr(xvar, "df") <- df
  attr(xvar, "by") <- by.var
  attr(xvar, "xorig") <- x
  attr(xvar, "lambda") <- lambda
  attr(xvar, "class") <- "smooth"
  attr(xvar,"smooth.name") <-"pvc"
  xvar
}

pvc.control <- function(inter = 20, degree = 3, order = 2, start = 10, quantiles = FALSE, 
                        k = 2, ...) 
{
  if (inter <= 0) {
    warning("the value of inter supplied is less than 0, the value of 10 was used instead")
    inter <- 10
  }
  if (degree <= 0) {
    warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
    degree <- 3
  }
  if (order < 0) {
    warning("the value of order supplied is zero or negative the default value of 2 was used instead")
    order <- 2
  }
  if (k <= 0) {
    warning("the value of GAIC/GCV penalty supplied is less than zero the default value of 2 was used instead")
    k <- 2
  }
  list(inter = inter, degree = degree, order = order, start = start, 
       quantiles = as.logical(quantiles)[1], k = k)
}


multicolumn <- function(formula,data){
  
  MF <- model.frame(formula=formula,data=data)
  a <- xtabs(formula=formula,data=MF)
  lev <- attributes(terms(formula))$"term.labels"
  lenv <- length(lev)
  ex1 <- paste("terms(formula)[[3]]",paste(rep("[[2]]",lenv-2),collapse=""),"[[3]]",sep="")
  ex2 <- paste("terms(formula)[[3]]",paste(rep("[[2]]",lenv-2),collapse=""),"[[2]]",sep="")
  y1 <- eval(parse(text=ex1))
  y2 <- eval(parse(text=ex2))
  lev1 <- nlevels(as.factor(MF[,colnames(MF)==y1]))
  lev2 <- nlevels(as.factor(MF[,colnames(MF)==y2]))
  if (lenv>2) {
    da <- apply(a,3:lenv,as.vector) 
    da <- matrix(da,ncol=lev1*lev2,byrow=TRUE)
  } else {
    da <- as.data.frame(t(as.vector(a)))
  }
  grid <- expand.grid(1:lev1,1:lev2)[,2:1]
  index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
  colnames(da) <- paste("n",index.lab,sep="")  
  if (lenv>2){ 
    dat <- MF[,which(colnames(MF)%in%lev[-(1:2)])]
    if (lenv==3) { 
      leve <- levels(factor(dat)) 
      xgrid <- data.frame(factor(leve))
      colnames(xgrid) <- lev[3]
    } else  {
      leve <- lapply(dat,function(x)levels(factor(x))) 
      xgrid <- eval(parse(text=paste("expand.grid(",paste(lev[-(1:2)],"=",
                  "leve[[",1:(lenv-2),"]]",collapse=",",sep=""),")",sep=""))) 
    }
    da <- cbind(da,xgrid)
  }
  attributes(da)$"resp" <- if (lenv>2) paste("cbind(",paste(colnames(da)[1:length(index.lab)],
                                                            collapse=","),")",sep="") else
                                                              paste("cbind(",paste(names(da)[1:length(index.lab)],
                                                                                   collapse=","),")",sep="")           
  attributes(da)$"Y1" <- y1
  attributes(da)$"ncat1" <- nlevels(as.factor(MF[,colnames(MF)==y1]))
  attributes(da)$"Y2" <- y2
  attributes(da)$"ncat2" <- nlevels(as.factor(MF[,colnames(MF)==y2]))
  return(da)
}

which.eq.smooth <- function(f,ik,ret=NULL){
  if (!is.null(f)) { 
    if (is.list(f)) {
      if (length(f)>=ik) ncol(attributes(f[[ik]])$"X") else ret  
    } else {
      if (ik == 1) ncol(attributes(f)$"X") else NULL 
    }
  }  else ret
}


make.Psmooth <- function(fsmooth,lc){
  if (!is.null(fsmooth)){  
    pslist <- mapply(function(x1,x2)  if (!((is.null(x1))&&(is.null(x2))))
                     lapply(1:length(x2), function(y) 
                       if (x2[y]>0) kronecker(diag(x2[y]),attributes(x1[[y]])$"D") else {
                                    if (is.list(x1)) attributes(x1[[y]])$"D" else attributes(x1)$"D" }),
                     fsmooth,lc,SIMPLIFY=FALSE)
  } else pslist <- as.matrix(0)
  return(pslist)
}

make.Psmooth2 <-  function(fsmooth,las,prop.smooth,nc,iB){
  if (!is.null(fsmooth)){       
    pslist <- lapply(iB,function(y) 
      mapply(function(x1,x2,x3,x4) {
        if ((!is.null(x2))&&(!is.na(x2[y]))){ 
          if (is.list(x1)) x2[y]*(attributes(x1[[y]])$"D") else       
          {if (is.null(x2)) as.matrix(0) else x2*(if (x3) attributes(x1)$"D" else 
            kronecker(diag(x4),attributes(x1)$"D"))}
        } else NULL},
        fsmooth,las,prop.smooth,nc,SIMPLIFY=FALSE))        
  } else pslist <- as.matrix(0)
  return(pslist)
}

invert.list <- function(x,maxB){
  which.convert <- sapply(x,function(y) (length(y)>0))
  new.x <- vector(mode="list",maxB) 
  for (j in 1:maxB){     
    for (i in 1:length(x)) {
      if (which.convert[i]) { 
        if ((is.list(x[[i]]))&&(length(x[[i]])>=j)) 
          new.x[[j]][[i]] <- x[[i]][[j]] else 
            if (j==1) new.x[[j]][[i]] <- x[[i]]
      }          
    } 
  }
  new.x
}

invert.list2 <- function(x,cond=TRUE){
  foo <- function(cond,ret) if (cond) (ret) else TRUE
  max.length <- max(sapply(x,length))
  new.x <-  vector(mode="list",max.length) 
  new.x <- lapply(new.x,c)
  for (j in 1:max.length){
    for (i in length(x):1){
      if ((!is.null(x[[i]]))&&foo(TRUE,!is.na(x[[i]][j]))) 
        new.x[[j]][i] <- x[[i]][j]
    } 
  }
  new.x
}


edf_df.pblm <- function(lams.int,bwb,pmat,cw,dl,wh,...){
  lams.int2 <- unlist(mapply(rep,lams.int,wh,SIMPLIFY=F))
  pbwb <- bwb + lams.int2*pmat
  diag.int <- colSums(solve(pbwb,tol=.Machine$double.eps^2)*bwb)
  edf_dl <- as.vector(diag.int %*% cw) - dl
  sum(edf_dl*edf_dl)
}


sym.tri <- function(x,q){
  X <- matrix(0,q,q)
  X[lower.tri(X,T)] <- x 
  tX <- t(X)
  utx <- upper.tri(X)
  X[utx] <- tX[utx]
  X
}



init.lams <- function(fsmooth,lams.name,prop,nc,wsnp){  
  if (!is.null(fsmooth)) {
    lams <- c()
    df <- c()
    maxS <- if (is.list(fsmooth)) length(fsmooth) else 1      
    k <- 0
    for (i in 1:maxS){  
      k <- k + 1
      a <- if (is.list(fsmooth)) attributes(fsmooth[[i]]) else attributes(fsmooth)
      if (!is.null(a$"lambda")){  
        .lams <- a$"lambda" 
        .df <- NA
      } else {
        .df <- if (is.null(a$"df")) 5 else a$"df"
        .lams <- if (!is.null(a$smooth.name)) 0 else
          log( abs((ncol(a$"X")-.df+1e-10)/(.df-(2-1e-10))))
      }
      if ((any(!wsnp))&&(!prop[!wsnp])[i]) {
        lams[k:(k+nc-1)] <- rep(.lams,nc)
        df[k:(k+nc-1)] <- rep(.df,nc)
        k <- (k+nc-1)
      } else {                               
        if (length(.lams)>1|length(.df)>1){ 
          lams <- .lams                     
          df <- .df                          
        } else {
          lams[k] <- .lams
          df[k] <- .df
       }
    }  }  
    names(lams) <- paste(lams.name," ",1:length(lams),sep="")
    attr(lams,"df") <- df
  }  else lams <- NULL
  return(lams)
}

propX<-function(.X,.xx,.pro,.nc,.ccv,.mat,.mat.RC,totrow){
  
  cxx <- cumsum(.xx)
  XN <- matrix(0,totrow,sum(.xx))
  XN[,1:cxx[1]] <- if (is.null(.mat.RC)) {
    if (.pro[1]) kronecker(.X[,1],.ccv) else kronecker(.X[,1],.mat)
  } else kronecker(.X[,1],.mat.RC) 
  if (any(.pro[-1]))   XN[,cxx[.pro]] <- kronecker(.X[,.pro],.ccv)
  np <- cxx[.pro]
  #if (length(.xx)==2) browser()
  if (any(!.pro[-1])) 
    XN[,c(rep(FALSE,cxx[1]),rep(!.pro[-1],diff(cxx)))] <- kronecker(as.matrix(.X[,-1])[,!.pro[-1]],.mat)
  return(XN)
}   

ridge <- function(x, df=1,lambda=NULL){
  if (!is.null(lambda)){
    if ((lambda<0)) {
      warnings("lambda cannot be negative, lambda=0 is chosen instead")
      lambda <- 0
    }
  }
  x.names <- colnames(x)
  x <- as.matrix(x)
  if (!is.null(df)){
    if ((df>ncol(x))||(df<=0)) {
      warnings("df out of the range (0,ncol(x)], df=1 is chosen instead")
      df <- 1
    }
  }
  D <- diag(ncol(x))
 xvar <- rep(0,nrow(x))    
 attr(xvar,"smooth.name") <- "ridge"  
 attr(xvar,"lambda") <- lambda
 attr(xvar,"df") <- df 
 attr(xvar, "D") <- D
 attr(xvar,"X") <- x
 attr(xvar,"names") <- x.names
 xvar
}

lasso <- function(x, df=1,lambda=NULL){
  if (!is.null(lambda)){
    if ((lambda<0)) {
      warnings("lambda cannot be negative, lambda=0 is chosen instead")
      lambda <- 0
    }
  }
  x.names <- colnames(x)
  x <- as.matrix(x)
  if (!is.null(df)){
    if ((df>ncol(x))||(df<=0)) {
      warnings("df out of the range (0,ncol(x)], df=1 is chosen instead")
      df <- 1
    }
  }
  D <- diag(ncol(x))
  xvar <- rep(0,nrow(x))    
  attr(xvar,"smooth.name") <- "lasso"  
  attr(xvar,"lambda") <- lambda
  attr(xvar,"df") <- df
  attr(xvar, "D") <- D
  attr(xvar,"X") <- x
  attr(xvar,"names") <- x.names
  xvar
}

model.frame.pblm <- function(formula,data,response=TRUE,center,scale,type,
                             weights,pro,ncat1,ncat2,fit=TRUE,ato,...){
  spec.smoothers <- c("pbs","pb","pvc","ridge","lasso")
  spec.functions <- c("arc2")
  fo.term <- terms(formula,specials=c("pbs","pb","pvc","ridge","lasso")) 
  X <- model.frame(as.formula(fo.term),data)
  res <- list()
  if (response){
    res$resp <- model.response(X)
    nresp <- ncol(res$resp)
    res$ynames <- colnames(res$resp)
    res$ncat1 <- if (nresp==2) 
                   nlevels(factor(res$resp[,1])) else {
                     if (!is.null(ncat1)) 
                       ncat1 else 
                       stop("argument ncat1 is required") 
                   }
    res$ncat2 <- if (nresp==2) 
                   nlevels(factor(res$resp[,2])) else 
                   as.integer(nresp/ncat1)
    if (nresp>2) 
      res$ta <- colSums(res$resp) else {
        res$ta <- if (is.null(weights)) 
                    table(res$resp[,1],res$resp[,2]) else
                    xtabs(as.formula(paste("weights~",
                                           paste(res$ynames,collapse="+"))),
                          data)
      }
    res$ncat1_1 <- res$ncat1-1
    res$ncat2_1 <- res$ncat2-1
    res$ncat <- res$ncat1*res$ncat2
    res$ncat_1 <- res$ncat-1
    res$ncat12 <- res$ncat1 + res$ncat2
    res$ncat12_1 <- res$ncat12-1
    res$ncat_1_2 <- (res$ncat1-1)*(res$ncat2-1)
    Y <- if (nresp==2) 
           expand.pblm(res$resp) else 
           res$resp
  } else {
    Y <- NULL
    if (fit==F){
      res$resp <- model.response(X)
      nresp <- ncol(res$resp)
      if(!is.null(res$resp)){
        res$ynames <- colnames(res$resp)
        res$ncat1 <- if (nresp==2) 
                       nlevels(factor(res$resp[,1])) else 
                       ncat1
        res$ncat2 <- if (nresp==2) 
                       nlevels(factor(res$resp[,2])) else 
                       as.integer(nresp/ncat1)
      } else{     
        res$ncat1 <- ncat1
        res$ncat2 <- ncat2
      }
      res$ncat1_1 <- res$ncat1-1
      res$ncat2_1 <- res$ncat2-1
      res$ncat <- res$ncat1*res$ncat2
      res$ncat_1 <- res$ncat-1
      res$ncat12 <- res$ncat1 + res$ncat2
      res$ncat12_1 <- res$ncat12-1
      res$ncat_1_2 <- (res$ncat1-1)*(res$ncat2-1)
    }
  }  
  u <- unlist(attr(fo.term,"specials"))
  if (!is.null(u))  u <- sort(u)
  us <- unlist(attr(fo.term,"specials")[spec.smoothers])
  uf <- unlist(attr(fo.term,"specials")[spec.functions])
  at <- rownames(attributes(attributes(X)$terms)$"factors")[u]
  atf <- rownames(attributes(attributes(X)$terms)$"factors")[uf]
  ats <- rownames(attributes(attributes(X)$terms)$"factors")[us]
  sm <- length(at)
  smf <- length(atf)
  sms <- length(ats)
  spec.name <- if (sm!=0) at else NULL
  spec.smooth <- if (sms!=0) ats else NULL
  spec.fun <- if (smf!=0) atf else NULL
  fsmooth <- if (!is.null(spec.smooth)) 
               X[,(colnames(X) %in% spec.smooth)] else
               NULL
  nby <- NULL
  X.pvc <- function(x){
    if (is.factor(attributes(x)$by)){
      fact <- attributes(x)$by
      lev <- levels(fact)
      nby <<- nlevels(fact)
      by <- model.matrix(~-1+fact)
      by <- split(by, rep(1:ncol(by), each = nrow(by)))
      Xlist <- lapply(1:nby,function(i)x)
      Xlist <- mapply(function(a,b){
                        attributes(a)$"X" <- attributes(a)$"X"*b;
                        a}, Xlist,by,SIMPLIFY=F)  
      Xlist <- lapply(1:nby,
                      function(i){
                        name <- attributes(Xlist[[i]])$"name";
                        name <- paste(name,lev[[i]],sep="");
                        attributes(Xlist[[i]])$"name" <- name;
                        Xlist[[i]]
                      })
      return(Xlist)
    } else {
      nby <<- 1
      attributes(x)$"Xor" <- attributes(x)$"X"
      attributes(x)$"X" <- attributes(x)$"X"*attributes(x)$by
      Xlist <- vector("list",1)
      Xlist[[1]] <- x
      return(Xlist)
    }
  }
  by.var <- list()
  BY.VAR <- if (response) as.list(rep(NA,length(colnames(X)[-1]))) else 
                          as.list(rep(NA,length(colnames(X))))
  names(BY.VAR) <- if (response) colnames(X)[-1] else colnames(X)
  by.var <- BY.VAR
   if (is.list(fsmooth)) {
     by.var <- lapply(fsmooth,function(x) 
                                if (is.null(attributes(x)$"by")) 
                                  NA else 
                                  attributes(x)$"by")
     for (ind in (names(BY.VAR))){
       if (ind %in% names(by.var)) 
         BY.VAR[[ind]] <- by.var[[ind]]
       
     }
     by.var <- BY.VAR
     names.fsmooth <- names(fsmooth)
     
     fsmooth <- lapply(1:length(fsmooth),
                       function(i){
                         attributes(fsmooth[[i]])$"name"<-names.fsmooth[i];
                         fsmooth[[i]]})
     fsm <- list()   
     nlevby <- c()
     k <- 0
     fsm1 <- name.pvc <- list()
     for (i in 1:length(fsmooth)) {
       if (is.null(attributes(fsmooth[[i]])$by)) {
         fsm0 <- list()
         fsm0[[1]] <- fsmooth[[i]] 
         names(fsm0) <- sapply(fsm0,
                               function(x) attributes(x)$"name")
         fsm <- c(fsm,fsm0)
       } else { 
         k <- k + 1
         fsm1[[k]] <- X.pvc(fsmooth[[i]])
         names(fsm1[[k]]) <- sapply(fsm1[[k]],function(x) attributes(x)$"name")
         
         name.pvc[[k]] <- names(fsm1[[k]])
         fsm <- c(fsm,fsm1[[k]])
         nlevby[length(nlevby)+1] <- nby 
       }
     }
     fsmooth <- fsm
   } else {

     by.var <- if (is.null(attributes(fsmooth)$"by")) 
                 NA else 
                 attributes(fsmooth)$"by"
     for (ind in (names(BY.VAR))){
       if (substr(ind,1,3)=="pvc")
         BY.VAR[[ind]] <- by.var

     }
     by.var <- BY.VAR
     nlevby <- c()
     attributes(fsmooth)$"name"<- spec.smooth
     if (is.null(attributes(fsmooth)$by)) {
       fsm <- list()
       fsm[[1]] <- fsmooth
       names(fsm) <- sapply(fsm,function(x)attributes(x)$"name")
       name.pvc <- names(fsm) 
     } else { 
       fsm <- list()
       fsm <- X.pvc(fsmooth)
       names(fsm) <- if (nby>1) 
                       sapply(fsm,function(x) attributes(x)$"name") else 
                       attributes(fsm)$"name"
       name.pvc <- names(fsm) 
       nlevby[length(nlevby)+1] <- nby 
     }    
     fsmooth <- fsm
   }
  ffun <- if (!is.null(spec.fun)) 
            X[,(colnames(X) %in% spec.fun)] else 
            NULL
  PRED <- (length(attr(X,"names"))>(response)) 
  mf <- model.frame(fo.term, data) 
  con <- colnames(mf)
  if (response) {
    con <- con[-1]
    mf <- subset(mf,select=con)
  }
 
  if (length(con)>0){
    is.fac <- sapply(1:length(con), 
                     function(i) (is.factor(mf[, i])||length(unique(mf[,i]))==2))
     fac.lev <- lapply(1:length(con), function(i) 
                                        if (is.factor(mf[, i])||(length(unique(mf[,i]))==2)) 
                                          levels(as.factor(mf[, i])) else 
                                          NA)  
    names(fac.lev) <- con
  } else {
    is.fac <- NULL
    fac.lev <- NULL  
  }
  X <- model.matrix(fo.term, data) 
  nlev <- attr(X, "assign")
  term <- terms(fo.term)
  termlabs <- attr(term,"term.labels")
  nlev0 <- lapply(1:length(unique(nlev)),
                  function(x) nlev[nlev %in% unique(nlev)[x]])
  i1 <- 0
  for (i in unique(nlev[-1])){
    if (substr(termlabs[i],1,3)=="pvc"){
      i1 <- i1 + 1
      nlev0[[i+1]] <- rep(unique(nlev[-1])[i],nlevby[i1])
    } 
  }
  nlev0 <- unlist(nlev0)

  if (!is.null(pro)){
    if (any(colnames(X)[!pro]%in% spec.name)) 
      stop("Sorry, category-dependent additive terms are not allowed")
    if (length(pro)!=ncol(X)) 
      stop("One component specified in argument <proportional> has wrong length")
    prop.smooth <- any(pro[colnames(X) %in% spec.smooth])
    prop.fun <- any(pro[colnames(X) %in% spec.fun])
  } else {
    prop.smooth <- TRUE
    prop.fun <- TRUE
  }

  if (PRED){
    nomi <- as.list(colnames(X))
    nomi2 <- nomi
    ss <- spec.smooth
    spec.smooth <- as.list(spec.smooth)
    spec.smooth2 <- spec.smooth
    k <- k0 <- k1 <- 0
    for (i in 1:length(nomi)){
      if (any(nomi[[i]] %in% ss)){
        if ((substr(nomi2[[i]],1,4) %in% "pvc(")){
          k1 <- k1 + 1
          k0 <- k0 + 1
          if (nlevby[k1]>1){
            k <- k + 1
            nomi[[i]] <-  if (is.list(name.pvc)) 
                            name.pvc[[k]] else 
                            name.pvc
            nomi2[[i]] <- rep(nomi2[[i]],nlevby[k])  
          }
        } else k0 <- k0 + 1
          spec.smooth[[k0]] <-  nomi[[i]] 
          spec.smooth2[[k0]] <- nomi2[[i]]
        } 
      }
    nomi2 <- unlist(nomi2)
    nomi <- unlist(nomi)
    spec.smooth <- unlist(spec.smooth)
    spec.smooth2 <- unlist(spec.smooth2)
    names(nomi2) <- NULL     
    names(nomi) <- NULL     
    names(spec.smooth2) <- NULL 
    names(spec.smooth) <- NULL  
    if (fit){
      X <- X[,nomi2,drop=F]
      eg <- qr(X)
      ok <- eg$pivot[1:eg$rank]      
    }else {
      X <- X[,nomi2,drop=F]
      ok <- 1:ncol(X[,nomi2,drop=F])
    } 
    X <- as.matrix(as.matrix(X)[,ok])
    colnames(X) <- colnames(X)[ok]
    which.ns <- !(nomi2 %in% ss)
    which.nf <- !(nomi2 %in% spec.fun)
    if (!is.null(dim(X)[2])) {
      X[,-1] <- scale(X[,-1],center=center,scale=scale)
    } else {
      X <- as.matrix(X)
      colnames(X) <- nomi
    }
  } else {
    ss <- spec.smooth
    nomi <- colnames(X)
    nomi2 <- unlist(sapply(nomi,
                           function(x) 
                             if (substr(x,1,4) %in% "pvc(") 
                               rep(x,nby) else 
                               x))
    names(nomi2) <- NULL     
    nomi <- unlist(sapply(nomi,
                          function(x) 
                            if (substr(x,1,4) %in% "pvc(") 
                              name.pvc else 
                              x))
    names(nomi) <- NULL     

    spec.smooth2 <- unlist(sapply(spec.smooth,
                                  function(x) 
                                    if (substr(x,1,4) %in% "pvc(") 
                                      name.pvc else 
                                      x,
                                  simplify=F))
    names(spec.smooth2) <- NULL 
    
    spec.smooth <- unlist(sapply(spec.smooth,
                                 function(x) 
                                   if (substr(x,1,4) %in% "pvc(") 
                                     rep(x,nby) else 
                                     x,
                                 simplify=F))
    names(spec.smooth) <- NULL  
    which.ns <- !(nomi %in% ss)
    which.nf <- !(nomi %in% spec.fun)   
    eg <- qr(X)
    ok <- eg$pivot[1:eg$rank]
    X <- as.matrix(as.matrix(X)[,ok])
    colnames(X) <- nomi[ok]
  }
  res1 <- c(list("X"=X,"Y"=Y,"ok"=ok,"nomi"=nomi,"fsmooth"=fsmooth,"by.var"=by.var,
                 "fo.term"=fo.term,"pivot"=if (fit) eg$pivot else NULL,
                 "which.ns"=which.ns,"spec.name"=spec.name, "nomi2"=nomi2,
                 "prop.smooth"=prop.smooth,"rank"=if (fit) eg$rank else NULL,
                 "pro"=pro,"spec.smooth"=spec.smooth2,"spec.fun"=spec.fun,
                 "prop.fun"=prop.fun,"ffun"=ffun,"which.nf"= which.nf,"nlev"=nlev,
                 "nlev0"=nlev0,"is.fac"=is.fac,"fac.lev"=fac.lev),res,"terms"=term)
  return(res1)
}

howmany <- function(x){
  hw <- c()
  ux <- unique(x)
  for (i in 1:length(ux)) hw[i] <- sum(x %in% ux[i])
  hw
}


make.mat <- function(nc,smooth=F){
  if (smooth) 
    list(c(rep(1,nc[[1]]),rep(0,Reduce('+',nc[2:3]))),
         c(rep(0,nc[[1]]),rep(1,nc[[2]]),rep(0,nc[[3]])),
         c(rep(0,Reduce('+',nc[1:2])),rep(1,nc[[3]]))) else
    list(rbind(diag(nc[[1]]),matrix(0,Reduce('+',nc[2:3]),nc[[1]])),
         rbind(matrix(0,nc[[1]],nc[[2]]),diag(nc[[2]]),matrix(0,nc[[3]],nc[[2]])),
         rbind(matrix(0,Reduce('+',nc[1:2]),nc[[3]]),diag(nc[[3]]))) 
  
}


IM2x2 <- function(Mp){
  Va <- Mp[2]*Mp[3]
  Vb <- Mp[4]*Mp[5]
  Vab <- 1/(sum(1/Mp[6:9]))
  Delta_p <- Mp[7]*Mp[8]-Mp[6]*Mp[9]
  VaVb <- Va*Vb
  Delta <- (VaVb-Delta_p*Delta_p)/VaVb
  D1 <- Va/Delta
  D2 <- Vb/Delta
  D12 <- Delta_p/Delta 
  matrix(c(D1,D12,0,D12,D2,0,0,0,Vab),3,3)
}

MIM2x2 <- function(Mp,w){
  Va <- Mp[,2]*Mp[,3]
  Vb <- Mp[,4]*Mp[,5]
  Vab <- 1/(rowSums(1/Mp[,6:9,drop=F]))
  Delta_p <- Mp[,7]*Mp[,8]-Mp[,6]*Mp[,9]
  VaVb <- Va*Vb
  Delta <- (VaVb-Delta_p*Delta_p)/VaVb
  D1 <- Va/Delta
  D2 <- Vb/Delta
  D12 <- Delta_p/Delta 
  mapply(function(x1,x12,x2,x3,x4) x4*matrix(c(x1,x12,0,x12,x2,0,0,0,x3),3,3),D1,D12,D2,Vab,w,SIMPLIFY=F)
}


dlogMult <- function (x,prob,N)
{
  r <- sum(lgamma(N + 1)) + sum(x * log(prob) - lgamma(x + 1))   
  r
}

C_Cdqrls <- getNativeSymbolInfo('Cdqrls', PACKAGE=getLoadedDLLs()$stats)

lamfun <- function(LA01,LA,LA1){
  LAM <- c(LA01,LA[names(LA) %in% names(LA1)])
  LAM <-LAM[order(names(LAM))]
  LAM
}    

make.bases <- function(fsmooth,maxB,m,ncat,mat.smooth,mat,prop.smooth,...){
  at.fsmooth <- vector(mode="list",maxB) 
  for (i in 1:maxB){
    at.fsmooth[[i]] <- lapply(fsmooth, function(x)
      if (is.list(x)){
        if (length(x)>=i) attributes(x[[i]])$"X" else NULL
      } else if (i==1) attributes(x)$"X")
  }
  
  af <- lapply(at.fsmooth, function(x){ 
    mapply(function(x1,x2,x3,x4){ 
      if (!is.null(x1)) {
        if (x4) kronecker(x1,x2) else 
          do.call("rbind",lapply(1:nrow(x1),function(i) kronecker(x3,t(x1[i,]))))
      } else NULL
    }, x, mat.smooth,mat,prop.smooth, SIMPLIFY=FALSE)
  }) 

  Bmat <- lapply(af,function(x) do.call("cbind",x))
  Bs <- lapply(Bmat,
               function(y){
                 lapply(split(y,rep(1:m, each=(ncat-1))),
                        function(x) matrix(x,ncat-1,ncol(y)))
               })
  ok <- lapply(Bmat,function(x)(colSums(abs(x))>0))
  Bmat <- mapply(function(x,y)x[,y],Bmat,ok,SIMPLIFY=FALSE)
  Bmat <- lapply(Bmat,function(x) Matrix(x,sparse=T))
  Bs <- mapply(function(y1,y2) lapply(y1,function(x)x[,y2]),Bs,ok,SIMPLIFY=FALSE)
  list("Bmat"=Bmat,"Bs"=Bs,"ok"=ok)
}

wh.eq.fun <- function(iB,fsmooth,prop.smooth,nc,ret=NULL,...){
  lapply(iB,function(w) 
    unlist(mapply(function(x,y,z,n) 
                    rep(which.eq.smooth(x,y,ret), if (z) 1 else n), 
                  fsmooth,w,prop.smooth,nc,SIMPLIFY=FALSE))) 
}   


inner.fit <- function(BM,NS,PW,OS,BW,eb,w.res,max.backfitting,rss.tol,...){
  w.res <- w.res - eb
  rss0 <- sum(w.res*w.res)
  for (kk in 1:max.backfitting){
    NS <- solve(PW+OS,BW%*%w.res)
    EL <- BM%*%NS
    etasm <- drop(eb + EL)
    rss <- sum((w.res-etasm)*(w.res-etasm))
    if (abs(rss-rss0)<(rss.tol*rss)) break
    if (kk==max.backfitting)
      warnings("Max number of iterations achieved in an inner cycle","\n")
    rss0 <- rss
  }
  list("EL"=EL,"NS"=NS,"etasm"=etasm,"w.res"=w.res)
}


additive.fit <- function(smoother,Bmat,new.smooth,PBWB,one.smooth,prodBW,etalist,w.res,
                         max.backfitting,rss.tol,maxB,m,ncat){
  if (smoother){
    for (ik in 1:maxB){
      eta.backfit <- if (length(etalist)>1) drop(Reduce("+",etalist[-ik])) else 0
      afit <- inner.fit(Bmat[[ik]],new.smooth[[ik]],PBWB[[ik]],one.smooth[[ik]],
                        prodBW[[ik]],eta.backfit,w.res,max.backfitting,rss.tol)
      new.smooth[[ik]] <- afit$NS
      etalist[[ik]] <- afit$EL
      w.res <- afit$w.res
    }
    etasmooth <- matrix(afit$etasm,m,ncat-1,byrow=TRUE)
    list("etalist"=etalist,"new.smooth"=new.smooth,"w.res"=w.res,"etasmooth"=etasmooth)
  } else list("etalist"=etalist,"new.smooth"=new.smooth,"w.res"=w.res,"etasmooth"=0)
}

make.ppmat <- function(BWB,ppmat,cwe,dflist,wh.eq,which.fixed.lams,lams.int,lams.start,
                      pgtol.df,lmm.df,factr.df,parscale.df,iter,...){
  
  if (iter>1) lams.start <- mapply(function(x1,x2){ x2[!is.na(x2)] <- x1; x2},
                                   lams.int,lams.start,SIMPLIFY=FALSE)        
  lams.int <- mapply(function(x1,x2,x3,x4,x5,x6,x7){
    if (length(x1[x7])==1){ 
      x1[x7] <- optimize(edf_df.pblm,lower=0,upper=1e5,
                         bwb=x2,pmat=x3,cw=x4,dl=x5[x7],wh=x6,
                         control=list(pgtol=pgtol.df,lmm=lmm.df,
                                      factr=factr.df,parscale=parscale.df))$minimum            
    } else {     
      if (length(x1[x7])>0){  
        x1[x7] <- optim(x1[x7],edf_df.pblm,method="L-BFGS-B",lower=0,upper=1e5,
                        bwb=x2,pmat=x3,cw=x4,dl=x5[x7],wh=x6,
                        control=list(pgtol=pgtol.df,lmm=lmm.df,factr=factr.df,
                                     parscale=rep(parscale.df,length(x1[x7]))))$par              
      } 
    };
    return(na.omit(x1))
  },lams.start, BWB, ppmat, cwe, dflist, wh.eq, which.fixed.lams,SIMPLIFY=FALSE)
  
  lams.vet <- mapply(unlist,mapply(rep, lams.int,na.omit(wh.eq),SIMPLIFY=FALSE),
                     SIMPLIFY=FALSE)
  PPmat <- mapply('*',lams.vet,ppmat,SIMPLIFY=F)
  list("lams.int"=lams.int,"PPmat"=PPmat)
}

 

pblm <- function(fo1=NULL, fo2=NULL, fo12=NULL, RC.fo=NULL, data, weights=NULL, 
                 type="gg", verbose=FALSE, contrasts=NULL, start=NULL, x=FALSE, 
                 center=FALSE, scale=FALSE, plackett=NULL, ncat1=NULL, ncat2=NULL, 
                 fit=TRUE, proportional=pblm.prop(...), penalty=pblm.penalty(...), 
                 control=pblm.control(...), ...){
  call <- match.call()
  if (missing(data)) data <- environment(fo1)
  mf <- match("weights",names(call))
  weights <- eval(call[[mf]],data)
  set0 <- pblm.prop()
  nmsC <- names(set0)
  set0[(namc <- names(proportional))] <- proportional
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in proportional: ", paste(noNms,collapse = ", "))     
  set1 <- pblm.control()
  nmsC <- names(set1)
  set1[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ", paste(noNms,collapse = ", "))  
  set2 <- pblm.penalty()
  nmsC <- names(set2)
  set2[(namc <- names(penalty))] <- penalty
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in penalty: ", paste(noNms,collapse = ", ")) 
  if (is.null(fo1)) stop("A model formula must be specified for fo1")
  conv.crit <- set1$conv.crit[1]
  vassign <- function(...,values,envir=parent.frame()){
    vars <- as.character(substitute(...()))
    for (i in seq_along(vars)){
      assign(vars[[i]],values,envir)
    }
  }
  vassign(X,Y,ok,nomi,fsmooth,ffun,spec.smooth,spec.fun,
          which.ns,which.nf,spec.name,prop.smooth,nomi2,
          nlev,nlev0, prop.fun,propos,pivot,rank,is.fac,
          fac.lev,term,values=vector(mode="list",3))
  by.var <- list(NA,NA,NA)
  fo.list <- list(fo1, fo2, fo12)
  vec.resp <- c(ifelse(fit,T,F),F,F)  
  if (!is.null(attr(data,"ncat1"))) {
    ncat1 <- attributes(data)$ncat1
    ncat2 <- attributes(data)$ncat2
  }  
  label1 <- attributes(data)$Y1
  label2 <- attributes(data)$Y2

  for (i in 1:3){
    if (!is.null(fo.list[[i]])){
      mfp <- model.frame.pblm(formula=fo.list[[i]],data=data,response=vec.resp[i],
                              center=center, scale=scale, type=type, weights=weights,
                              pro=set0[[i]],ncat1=ncat1,ncat2=ncat2,fit=fit,ato=i)
      X[[i]] <- mfp$X
      ok[[i]] <- mfp$ok
      nomi[[i]] <- mfp$nomi
      which.ns[[i]] <- mfp$which.ns
      which.nf[[i]] <- mfp$which.nf
      spec.name[[i]] <- mfp$spec.name
      spec.smooth[[i]] <- mfp$spec.smooth
      spec.fun[[i]] <- mfp$spec.fun
      pivot[[i]] <- mfp$pivot
      rank[[i]] <- mfp$rank
      nomi2[[i]] <- mfp$nomi2
      prop.smooth[[i]] <- mfp$prop.smooth      
      prop.fun[[i]] <- mfp$prop.fun
      is.fac[[i]] <- mfp$is.fac
      fac.lev[[i]] <- mfp$fac.lev
      term[[i]] <- mfp$terms
      nlev[[i]] <- mfp$nlev
      nlev0[[i]] <- mfp$nlev0
      if (length(mfp$by.var)>0) by.var[[i]] <- mfp$by.var
      if ((!is.null(mfp$fsmooth))&&length(mfp$fsmooth)>0) 
         fsmooth[[i]] <- mfp$fsmooth
      if (!is.null(mfp$ffun)) ffun[[i]] <- mfp$ffun
    } else {
      if (i==1) stop("fo1 must be specified")
      X[[i]] <- X[[1]]
      ok[[i]] <- ok[[1]]
      nomi[[i]] <- nomi[[1]]
      which.ns[[i]] <- which.ns[[1]]
      which.nf[[i]] <- which.nf[[1]]
      spec.name[[i]] <- spec.name[[1]]
      spec.smooth[[i]] <- spec.smooth[[1]]
      spec.fun[[i]] <- spec.fun[[1]]
      pivot[[i]] <- pivot[[1]]
      rank[[i]] <- rank[[1]]      
      nomi2[[i]] <- nomi2[[1]]      
      prop.smooth[[i]] <- prop.smooth[[1]]
      prop.fun[[i]] <- prop.fun[[1]]
      is.fac[[i]] <- is.fac[[1]]
      fac.lev[[i]] <- fac.lev[[1]]
      term[[i]] <- term[[1]]
      nlev[[i]] <- nlev[[1]]
      nlev0[[i]] <- nlev0[[1]]
      by.var[[i]] <- by.var[[1]]
      if(length(by.var[[1]])>0) by.var[[i]] <- by.var[[1]]
      if (!is.null(fsmooth[[1]])) fsmooth[[i]] <- fsmooth[[1]] 
      if (!is.null(ffun[[1]])) ffun[[i]] <- ffun[[1]]
    }
   propos[[i]] <- if (!is.null(set0[[i]])) rep(set0[[i]],howmany(nomi2[[i]])) else 
                     c(FALSE, if (length(nomi[[i]])>1) rep(TRUE,length(nomi[[i]])-1) else NULL)
   set0[[i]] <-  if (!is.null(set0[[i]])) set0[[i]][ok[[i]]] else 
                    c(FALSE, if (dim(X[[i]])[2]>1) rep(TRUE,ncol(X[[i]])-1) else NULL) 
   
    if (i==1) {
      Y <- mfp$Y
      resp <- mfp$resp
      ynames <- mfp$ynames
      ta <- mfp$ta 
      ncat1 <- mfp$ncat1 
      ncat2 <- mfp$ncat2 
      ncat1_1 <- mfp$ncat1_1 
      ncat2_1 <- mfp$ncat2_1 
      ncat <- mfp$ncat 
      ncat_1 <- mfp$ncat_1 
      ncat12 <- mfp$ncat12 
      ncat12_1 <- mfp$ncat12_1
      ncat_1_2 <- mfp$ncat_1_2  
    }  
  }   
  fun <- any(!unlist(which.nf))
  pnl.type <- if (set2$pnl.type!="none") set2$pnl.type else NULL
  if (penalty$pnl.type == "equal"){
    if (ncat1!=ncat2) 
      stop("penalty type not allowed because the responses have different numbers of categories ")
    if ((!is.null(fo2))&(!identical(fo1,fo2))) 
      stop("penalty type not allowed because fo1 and fo2 are not equal")
    if (!identical(set0[[1]],set0[[2]])) 
      stop("penalty type not allowed because propo[[1]] and propo[[2]] are not equal")
  } 
  if (is.null(plackett)) {
    plackett <- ifelse(((type=="gg")||((ncat1==2L)&(ncat2==2L))), TRUE, FALSE)
  }
  if ((type!="gg")&(plackett==TRUE)&((ncat1>2L)|(ncat2>2L))){
    warning("plackett=TRUE is allowed only when type=gg or when both the responses 
             are binary. A Newton-Raphson scheme is used instead.")
    plackett <- FALSE
  }
  M <- make.marginals(ncat1,ncat2,type)
  C <- make.contrasts(ncat1,ncat2)
  propo <- set0
  if (!is.null(RC.fo)) propo[[3]][1] <- FALSE
  if (is.null(set2$s1)) set2$s1 <- rep(1,sum(propo[[1]]==FALSE)) else
  {if (length(set2$s1)!=sum(propo[[1]]==FALSE))
    stop("length of s1 does not match with length(propo[[1]]==FALSE)")}
  if (is.null(set2$s2)) set2$s2<-rep(1,sum(propo[[2]]==FALSE)) else
  {if (length(set2$s2)!=sum(propo[[2]]==FALSE))
    stop("length of s2 does not match with length(propo[[2]]==FALSE)")}
  if (is.null(set2$s3)) set2$s3<-rep(1,sum(propo[[3]]==FALSE)) else
  {if (length(set2$s3)!=sum(propo[[3]]==FALSE))
    stop("length of s3 does not match with length(propo[[3]]==FALSE)")}
  if (is.null(set2$s4)) set2$s4<-rep(1,sum(propo[[3]]==FALSE)) else
  {if (length(set2$s4)!=sum(propo[[3]]==FALSE))
    stop("length of s4 does not match with length(propo[[3]]==FALSE)")}
  m <- nrow(X[[1]])
  nc <- list(ncat1_1,ncat2_1,ncat_1_2)
  smoother <- any(!unlist(which.ns))
  n.smoothers <- sum(sapply(which.ns,function(x)sum(!x)))
  lamsa <- mapply(init.lams,fsmooth,c("lamsa[[1]]","lamsa[[2]]","lamsa[[3]]"),
                  propos,nc,which.ns,SIMPLIFY=FALSE) 
  lamfa <- mapply(init.lams,ffun,c("lam1","lam2","lam3"),
                  propo,nc,which.nf,SIMPLIFY=FALSE)
  if (fun) {
    set2$lam1 <- rep(0,length(propo[[1]]))
    set2$lam1[!which.nf[[1]]] <- lamfa[[1]]
    set2$lam2 <- rep(0,length(propo[[2]]))
    set2$lam2[!which.nf[[2]]] <- lamfa[[2]]
    set2$lam3 <- rep(0,length(propo[[3]]))
    set2$lam3[!which.nf[[3]]] <- lamfa[[3]]
  }
  df <- lapply(lamsa,function(x) attributes(x)$df)
  un.df <- unlist(df)
  dffun <- lapply(lamfa,function(x) attributes(x)$df)
  which.fixed.df <- sapply(un.df,function(x)!is.na(x))
  which.fixed.lams <- rep(TRUE,length(which.fixed.df))
  las <- lamsa
  lc <- mapply(function(x0,x1,x2) 
               if ((length(x1)>1)&&(!all(x1))) (!x0[(!x1)])*x2 else NULL,
               propos,which.ns,nc,SIMPLIFY=F)
  mat <- make.mat(nc,smooth=FALSE)
  RC <- if (is.null(RC.fo)) NULL else make.RC(RC.fo,ncat1_1,ncat2_1,contrasts)
  mat.RC <- if (is.null(RC.fo)) vector("list",3) else 
              list(NULL,NULL,rbind(matrix(0,Reduce('+',nc[1:2]),ncol(RC)),RC)) 
  if (smoother){ 
    mat.smooth <- make.mat(nc,smooth=TRUE) 
    maxB <- max(sapply(which.ns,function(x)sum(!x)))
    iB <- as.list(1:maxB)
    ncl <- list(ncat1_1,ncat2_1,ncat_1_2)
    wh.eq2 <- wh.eq.fun(iB,fsmooth,prop.smooth,nc,ret=NA) 
    wh.eq <- wh.eq.fun(iB,fsmooth,prop.smooth,nc,ret=NULL)
    at.fsmooth <- vector(mode="list",maxB) 
    for (i in 1:maxB){
      at.fsmooth[[i]] <- lapply(fsmooth, function(x)
        if (is.list(x)){
          if (length(x)>=i) attributes(x[[i]])$"X" else NULL
        } else if (i==1) attributes(x)$"X")
    }

    Base <- make.bases(fsmooth,maxB,m,ncat,mat.smooth,mat,prop.smooth)
    Bmat <- Base$Bmat
    Bs <- Base$Bs
    ok <- Base$ok
  } else  Bmat <- maxB <- wh.eq <- wh.eq2 <- lams <- NULL
  if (is.null(weights)) weights <- rep(1,m)
  xx1<-propo[[1]]
  xx1[!xx1] <- ncat1_1
  xx2<-propo[[2]]
  xx2[!xx2] <- ncat2_1
  if (is.null(RC.fo)){
    xx12 <- propo[[3]]
    xx12[!xx12] <- ncat_1_2
  } else {
    xx12 <-propo[[3]][-1]
    xx12[!xx12] <- ncat_1_2
    xx12 <- c(ncol(RC),xx12)
  }
  xx <- list(xx1,xx2,xx12)
  xx1s<-propos[[1]]
  xx1s[!xx1s] <- ncat1_1
  xx2s<-propos[[2]]
  xx2s[!xx2s] <- ncat2_1
  if (is.null(RC.fo)){
    xx12s <- propos[[3]]
    xx12s[!xx12s] <- ncat_1_2
  } else {
    xx12s <-propos[[3]][-1]
    xx12s[!xx12s] <- ncat_1_2
    xx12s <- c(ncol(RC),xx12s)
  }
  xxs <- list(xx1s,xx2s,xx12s)
  d <- split(diag(3),1:3)
  ccv <- lapply(d, function(x) unlist(mapply(rep,as.list(x),nc))) 
  sumcat <- Reduce('+',nc)
  totrow <- m*sumcat
  Zlist <- mapply(propX,X,xx,set0,nc,ccv,mat,mat.RC,
                  MoreArgs=list(totrow),SIMPLIFY=F)
  respStr <- mapply(function(x,y) if (ncol(as.matrix(x[,1:y[1]]))>1) 
                                     rowSums(as.matrix(x[,1:y[1]])) else 
                                     as.matrix(x[,1:y[1]]), Zlist,xxs)
  Z <- do.call("cbind",Zlist)
  totcol <- sum(c(xx1,xx2,xx12))
  margins <- rep(c(rep(TRUE,(ncat1_1+ncat2_1)),rep(FALSE,ncat_1_2)),m)
  sw <- sum(weights)
  if (!is.null(set2$pnl.type)){
    if (!is.null(set2$lam1)){
      lam1 <- set2$lam1
      if (any(lam1<0)) {
        lam1[lam1<0] <- 0
        warning("lam1 cannot be negative, lam1=0 is to be fixed instead")
      }  
      names(lam1)<-paste("lam1.",1:length(lam1),sep="")
      la01 <- lam1[lam1<=set2$min.lam.fix]
      la1 <- lam1[lam1>set2$min.lam.fix]
    }   else la1 <-la01 <- lam1 <- NULL
    if (!is.null(set2$lam2)){
      lam2 <- set2$lam2
      if (any(lam2<0)) {
        lam2[lam2<0] <- 0
        warning("lam2 cannot be negative, lam2=0 is to be fixed instead")
      }      
      names(lam2)<-paste("lam2.",1:length(lam2),sep="")
      la02 <- lam2[lam2<=set2$min.lam.fix]
      la2 <- lam2[lam2>set2$min.lam.fix]
    }  else la2 <-la02 <- lam2 <- NULL    
    if (!is.null(set2$lam3)){
      lam3 <- set2$lam3
      if (any(lam3<0)) {
        lam3[lam3<0] <- 0
        warning("lam3 cannot be negative, lam3=0 is to be fixed instead")
      }  
      names(lam3)<-paste("lam3.",1:length(lam3),sep="")
      la03 <- lam3[lam3<=set2$min.lam.fix]
      la3 <- lam3[lam3>set2$min.lam.fix]
    } else la3 <-la03 <- lam3 <- NULL   
    if ((!is.null(pnl.type))&&(pnl.type=="ARC2")&(!is.null(set2$lam4))){
      lam4 <- set2$lam4
      if (any(lam4<0)) {
        lam4[lam4<0] <- 0
        warning("lam4 can not be negative, lam4=0 is to be fixed instead")
      }  
      names(lam4)<-paste("lam4.",1:length(lam4),sep="")
      la04 <- lam4[lam4<=set2$min.lam.fix]
      la4 <- lam4[lam4>set2$min.lam.fix]
      lam <- c(la1,la2,la3,la4)
    } else {
      lam4 <- NULL
      lam <- c(la1,la2,la3)
    }
  } else lam <- NULL  
  lambda <- c(lam,lams <- unlist(las))
  if (set2$constraints) constr <- DIFF(ncat1,ncat2,m)
  wY <- weights*Y  
  if (!fit) {
    rval <- if (smoother)
    list("Y"=Y,"Z"=Z,"Xlist"=X,"m"=m,"ta"=ta,"C"=C,"M"=M,"Zlist"=Zlist,
         "nlev0"=nlev0,"Bmat"=Bmat,"Bs"=Bs,"fsmooth"=fsmooth,"ok"=ok,
         "terms"=term,"nlev"=nlev,"by.var"=by.var,"label1"=label1,
         "label2"=label2,"xxs"=xxs) else
    list("Y"=Y,"Z"=Z,"Xlist"=X,"m"=m,"ta"=ta,"C"=C,"M"=M,"ok"=ok,
         "terms"=term,"nlev"=nlev,"nlev0"=nlev0,"Zlist"=Zlist,
         "label1"=label1,"label2"=label2,"xxs"=xxs)  
   return(rval)
  } else {
    Xlist <- X  
    X <- aperm(array(t(Z),c(ncol(Z),ncat_1,m)) ,c(2,1,3))
    taz <- ta+set1$zero.adj
    prob <- as.vector(t(taz/sum(taz)))
    veroeta <- crossprod(C,log(M%*%prob))
    if (smoother) {
      etalist <- vector(mode="list",maxB) 
      new.smooth <- vector(mode="list",maxB)
      select.lams <- any(which.fixed.df)
      if (select.lams){
        splist <- make.Psmooth(fsmooth,lc)
        pslist <- invert.list(splist,maxB)      
        dlf <- if (any(!unlist(prop.smooth))) 
                 lapply(un.df,function(x)x) else 
                 df
        dflist <- invert.list2(dlf, FALSE)
        which.fixed.lams <- lapply(dflist,function(x)!is.na(x))
        Pmat <- if (length(pslist)>1) 
                  lapply(pslist, Mdiag) else 
                list(Mdiag(pslist[[1]]))
        cwe <- lapply(wh.eq,function(x) 
                              Mdiag(lapply(x,function(y) 
                                               if (!is.na(y)) matrix(1,y))))
        ppmat <- lapply(Pmat,crossprod)
      }
    }    
    if (is.null(start)) {
      par1 <- c(veroeta[2:(xx1[1]+1)], if (length(xx1)>1) rep(0,sum(xx1[-1])))
      par2 <- c(veroeta[(xx1[1]+2):(xx1[1]+xx2[1]+1)], 
                if (length(xx2)>1) rep(0,sum(xx2[-1])))
      par12 <- if (!is.null(RC.fo))  rep(0,sum(xx12)) else
      {
        if (xx12[1]==1) c(mean(veroeta[(xx1[1]+xx2[1]+2):(length(veroeta))]),
                          rep(0,sum(xx12[-1]))) else
                            c(veroeta[(xx1[1]+xx2[1]+2):(length(veroeta))], 
                              if (length(xx12)>1) rep(0,sum(xx12[-1])))
      }
      start <- c(par1,par2,par12)
    }  
    GRAD <- find.score(C,M,prob)
    GRAD.INV <- array(solve(GRAD)[,-1],c(ncat,ncat_1,m))  
    GRAD <- array(GRAD,c(ncat,ncat,m))            
    etasmooth <- matrix(0,m,ncat_1)
    NOTCONV <- TRUE
    toll <- rep(1,totcol)
    zerovect <- rep(0,m)
    neweta <- cbind(zerovect,matrix(Z%*%start,m,ncat_1,byrow=TRUE))
    sum_row_Y <- rowSums(Y+set1$zero.adj) 
    propY <- (Y+set1$zero.adj)/sum_row_Y
    wei_sum_row_Y <- weights*sum_row_Y 
    Xt <- aperm(X ,c(2,1,3))
    if (smoother){  
      one.smooth <- vector("list",length(wh.eq))
      for (i in 1:length(wh.eq)){
        if (is.null(pnl.type)||(pnl.type!="equal")){
          one.smooth[[i]] <- matrix(0,sum(wh.eq[[i]]),sum(wh.eq[[i]]))
        } else {
          cb <- cbind(diag(wh.eq[[i]][1]),-diag(wh.eq[[i]][2]),
                      if (is.na(wh.eq[[i]][3])) NULL else 
                        matrix(0,wh.eq[[i]][3],wh.eq[[i]][3]))
          one.smooth[[i]] <- Mdiag(list(cb))
        }
      }
    } else one.smooth <- NULL
    bin <- (ncat1==2)&(ncat2==2)
    Ct <- t(C)
  
    if (set1$auto.select) {
      which.lambdas <- c(if (length(lam>0)) rep(T,length(lam)),unlist(which.fixed.lams)) 
      opt <- optim(lambda[which.lambdas],pblm.fit,method="L-BFGS-B",
                   lower=rep(0,length(lambda[which.lambdas])),
                   upper=rep(1e05,length(lambda[which.lambdas])),
                   pnl.type=pnl.type, lam1=set2$lam1, la01=la01, lam2=set2$lam2,
                   la02=la02, lam3=set2$lam3, la03=la03, lam4=set2$lam4, la04=la04, 
                   la1=la1,la2=la2,la3=la3,la4=la4,totcol=totcol, propo=propo, xx1=xx1,
                   xx2=xx2, xx12=xx12, ncat1=ncat1, ncat2=ncat2, s1=set2$s1, s2=set2$s2, 
                   s3=set2$s3,s4=set2$s4,constraints=set2$constraints,lamv1=set2$lamv1, 
                   lamv2=set2$lamv2, ncat12_1=ncat12_1, Z=Z, margins=margins, constr=constr, 
                   m=m, type=type, smoother=smoother, fsmooth=fsmooth, lamsa=lamsa,
                   maxB=maxB, which.fixed.df=which.fixed.df,df=df, wh.eq=wh.eq, ok=ok, 
                   prop.smooth=prop.smooth, ncl=ncl, iB=iB, start=start, veroeta=veroeta,
                   RC.fo=RC.fo, C=C, M=M, prob=prob, ncat=ncat, ncat_1=ncat_1, wY=wY,
                   l=set1$l, weights=weights, Y=Y, X=X, Bs=Bs, Bmat=Bmat, pgtol.df=set1$pgtol.df, 
                   lmm.df=set1$lmm.df, factr.df=set1$factr.df, parscale.df=set1$parscale.df,
                   max.backfitting=set1$max.backfitting,plackett=plackett,
                   ncat12=ncat12, ncat1_1=ncat1_1, ncat2_1=ncat2_1, zero.adj=set1$zero.adj,
                   min.step.l=set1$min.step.l, restore.l=set1$restore.l, acc2=set1$acc2, 
                   maxit2=set1$maxit2,lc=lc, conv.crit=conv.crit,grad.tol=set1$grad.tol,
                   acc=set1$acc, maxit=set1$maxit, verbose=verbose, auto.select=set1$auto.select, 
                   spec.name=set1$spec.name, gaic.m=set1$gaic.m, rss.tol=set1$rss.tol, 
                   control=list(maxit=set1$max.gaic.iter,pgtol=set1$pgtol.gaic,
                                lmm=set1$lmm.gaic,factr=set1$factr.gaic,
                               parscale=rep(set1$parscale,length(lambda[which.lambdas]))),
                   Pmat=Pmat,  which.fixed.lams=which.fixed.lams,
                   select.lams=select.lams,dflist=dflist, cwe=cwe,ppmat=ppmat,GRAD=GRAD,GRAD.INV=GRAD.INV,
                   etasmooth=etasmooth,NOTCONV=NOTCONV,toll=toll,zerovect=zerovect,neweta=neweta,
                   sum_row_Y=sum_row_Y,propY=propY,wei_sum_row_Y=wei_sum_row_Y,Xt=Xt,
                   one.smooth=one.smooth,bin=bin,Ct=Ct)  

      lambda[which.lambdas] <- opt$par
    }  else opt <- NULL
    fit <- pblm.fit(la=lambda, pnl.type=pnl.type, lam1=set2$lam1, la01=la01, 
                    lam2=set2$lam2, la02=la02, lam3=set2$lam3, la03=la03, lam4=set2$lam4, 
                    la04=la04,la1=la1,la2=la2,la3=la3,la4=la4, totcol=totcol, propo=propo,
                    xx1=xx1,xx2=xx2, xx12=xx12, ncat1=ncat1, ncat2=ncat2, s1=set2$s1, s2=set2$s2, 
                    s3=set2$s3,lc=lc,s4=set2$s4, constraints=set2$constraints, lamv1=set2$lamv1, 
                    lamv2=set2$lamv2, ncat12_1=ncat12_1, Z=Z, margins=margins, constr=constr, m=m, 
                    type=type, smoother=smoother, fsmooth=fsmooth, lamsa=lamsa, maxB=maxB,
                    which.fixed.df=which.fixed.df,df=df, wh.eq=wh.eq, ok=ok, ncl=ncl, iB=iB,
                    prop.smooth=prop.smooth, start=start, veroeta=veroeta, RC.fo=RC.fo,
                    C=C, M=M, prob=prob, ncat=ncat, ncat_1=ncat_1, wY=wY, l=set1$l,weights=weights,
                    Y=Y, X=X, Bs=Bs,  Bmat=Bmat, pgtol.df=set1$pgtol.df, lmm.df=set1$lmm.df,
                    factr.df=set1$factr.df,  parscale.df=set1$parscale.df,grad.tol=set1$grad.tol,
                    max.backfitting=set1$max.backfitting,
                    ncat12=ncat12, ncat1_1=ncat1_1, ncat2_1=ncat2_1, zero.adj=set1$zero.adj,
                    min.step.l=set1$min.step.l, restore.l=set1$restore.l, acc2=set1$acc2, 
                    maxit2=set1$maxit2, acc=set1$acc, maxit=set1$maxit, verbose=verbose, 
                    auto.select=set1$auto.select, spec.name=set1$spec.name, gaic.m=set1$gaic.m, 
                    rss.tol=set1$rss.tol, plackett=plackett, conv.crit=conv.crit,
                    Pmat=Pmat,which.fixed.lams=which.fixed.lams,GRAD=GRAD,GRAD.INV=GRAD.INV,
                    select.lams=select.lams,dflist=dflist, cwe=cwe,ppmat=ppmat,
                    etasmooth=etasmooth,NOTCONV=NOTCONV,toll=toll,zerovect=zerovect,neweta=neweta,
                    sum_row_Y=sum_row_Y,propY=propY,wei_sum_row_Y=wei_sum_row_Y,Xt=Xt,
                    one.smooth=one.smooth,bin=bin,Ct=Ct)
    
    xnames <- vector(mode="list",3)
    xnames[[1]] <- rep(nomi[[1]],xx1s)
    seq1 <- c()
    for (i in 1:length(xx1s)) seq1 <- c(seq1,1:xx1s[i])
    logit1 <- paste("logit1.",seq1,": ",xnames[[1]][1:sum(xx1s)],sep="")
    seq2 <- c()
    for (i in 1:length(xx2s)) seq2<-c(seq2,1:xx2s[i])
    xnames[[2]] <- rep(nomi[[2]],xx2s)
    logit2 <- paste("logit2.",seq2,": ",xnames[[2]][1:sum(xx2s)],sep="")
    seq12 <- c()
    for (i in 1:length(xx12s)) seq12 <- c(seq12,1:xx12s[i])
    if (is.null(RC.fo)){
      xnames[[3]] <- rep(nomi[[3]],xx12s)
      GOR12 <- paste("lor12.",seq12,": ",xnames[[3]][1:sum(xx12s)],sep="")
    }  else {
      col.non.rc <- sum(xx12s[-1])  
      xnames[[3]] <- c(colnames(RC),rep(nomi[[3]][-1],xx12s[-1]))
      GOR12 <- paste("lor12.",seq12,": ",xnames[[3]],sep="")
    }
    nomiresp <- colnames(resp) 
    w.res <- matrix(attributes(fit)$"w.res",m,ncat-1,byrow=TRUE) 
    z <- matrix(attributes(fit)$"z",ncol=ncat-1)
    coef <- rep(NA,length(unlist(mapply(function(x,y)rep(x,y),pivot,xxs,SIMPLIFY=T))))
    coef[unlist(mapply(function(x,y,z)rep(1:length(x) %in% x[1:y],z),pivot,rank,xxs,SIMPLIFY=F))] <- attributes(fit)$"new.par"
    names(coef) <- c(paste(substr(type,1,1),".",logit1,sep=""),
                     paste(substr(type,1,1),".",logit2,sep=""),
                     paste(type,".",GOR12,sep=""))     
    maxNpred <- max(sapply(1:3,function(j) length(attr(terms(term[[j]]),"term.labels") )))
    which.term <- matrix(FALSE,nrow=3,ncol=maxNpred)
    wt <- lapply(term,function(x) length(attributes(x)$"term.labels")) 
    for (i in 1:3){
      if (wt[[i]]>0) which.term[i,1:wt[[i]]] <- TRUE
    } 
   # browser()
    rval <- list("coef"=coef,"n"=sw,"m"=m,"p"=attributes(fit)$"p",
                 "Y"=wY,"x"=Matrix(Z,sparse=TRUE),"xx1"=xx1,"xx2"=xx2,"xx12"=xx12,
                 "xx1s"=xx1s,"xx2s"=xx2s,"xx12s"=xx12s,"ynames"=ynames,
                 "tol"=attributes(fit)$"tol","llp"=attributes(fit)$"llp",
                 "ll"=attributes(fit)$"ll","dev"=attributes(fit)$"dev",
                 "devp"=attributes(fit)$"devp","IM"=attributes(fit)$"IM",
                 "IMp"=attributes(fit)$"IMp","convergence"=!attributes(fit)$"NOTCONV",
                 "iter"=attributes(fit)$"iter", "maxB"=maxB, "maxNpred"=maxNpred, "ncat1"=ncat1,
                 "ncat2"=ncat2,"ncat"=ncat,"nc"=unlist(nc),"weights"=weights,"P"=attributes(fit)$"P",          
                 "gaic.m"=set1$gaic.m,"lam1"=attributes(fit)$"lam1",
                 "lam2"=attributes(fit)$"lam2", "lam3"=attributes(fit)$"lam3",
                 "lam4"=attributes(fit)$"lam4", "lambda"=attributes(fit)$"lambda",
                 "opt"=opt,"etasmooth"=attributes(fit)$"etasmooth",
                 "eta"=attributes(fit)$"neweta",
                 "fsmooth"= fsmooth,"one.smooth"=attributes(fit)$one.smooth,
                 "df.fix"=attributes(fit)$"df.fix","df.fix.vect"=attributes(fit)$"df.fix.vect",
                 "df.smooth"=attributes(fit)$"df.smooth","w.res"=w.res, 
                 "W2"=attributes(fit)$"W2","z"= z,"any.smoother"=smoother,
                 "Bmat"=attributes(fit)$"Bmat","wh.eq"=attributes(fit)$"wh.eq",
                 "PBWB"=attributes(fit)$"PBWB","BWB"=attributes(fit)$"BWB",
                 "etalist"=attributes(fit)$"etalist","spec.name"=spec.name,
                 "spec.smooth"=spec.smooth,"beta.smooth"=attributes(fit)$"new.smooth", 
                 "n.smoothers"=n.smoothers,"PPmat"=attributes(fit)$"PPmat",
                 "wh.eq2"=wh.eq2,"pnl.type"=pnl.type, "spec.fun"=spec.fun,
                 "xnames"=xnames,"prop.smooth"=prop.smooth,"which.ns"=which.ns,"which.nf"=which.nf,
                 "GAIC"=attributes(fit)$"GAIC","ta"=ta,"pivot"=pivot,"rank"=rank,
                 "set0"=set0,"fo.list"=fo.list,"center"=center,"scale"=scale,"type"=type,
                 "plackett"=plackett,"RC.fo"=RC.fo,"contrasts"=contrasts,"acc2"=set1$acc2,
                 "maxit2"=set1$maxit2,"X"=Xlist,"is.fac"=is.fac,"fac.lev"=fac.lev,"terms"=term,
                 "by.var"=by.var,"label1"=label1,"label2"=label2,"which.term"=which.term,
                 "respStr"=respStr)
    class(rval) <- "pblm"
    rval$call <- call
    rval
  }
}


pblm.fit <- function(la, pnl.type, lam1, la01, lam2, la02, lam3, la03, lam4, 
                     la1, la2, la3, la4, la04, totcol, propo, xx1, xx2, xx12,                 
                     ncat1, ncat2, s1, s2, s3, s4, constraints, lamv1, lamv2, lc, 
                     ncat12_1, Z, margins, constr, m, type, smoother, fsmooth, lamsa,
                     maxB, which.fixed.df, df, wh.eq, ok,
                     Bmat, prop.smooth, ncl, iB, start, veroeta, RC.fo, C, M, prob,
                     ncat, ncat_1, wY, l, weights, Y, X, Bs, pgtol.df, lmm.df,
                     factr.df, parscale.df, max.backfitting, ncat12, ncat1_1, ncat2_1, zero.adj,
                     min.step.l, restore.l, acc2, maxit2, acc, maxit, verbose, 
                     auto.select, spec.name, gaic.m, rss.tol, plackett, conv.crit,
                     which.fixed.lams, select.lams, dflist, cwe, ppmat, GRAD, GRAD.INV,
                     etasmooth, NOTCONV, toll, zerovect, neweta, sum_row_Y, propY,
                     wei_sum_row_Y, Xt, one.smooth, bin,Ct=Ct,grad.tol,...){

  Pmat <- PBWB <- BWB <- BW <- new.smooth <- etalist <- NULL
  if (!is.null(pnl.type)){    
    if (!is.null(lam1)) lam1 <- lamfun(la01,la,la1)
    if (!is.null(lam2)) lam2 <- lamfun(la02,la,la2)
    if (!is.null(lam3)) lam3 <- lamfun(la03,la,la3)
    if (pnl.type=="ARC2"){
      lam4 <- if (!is.null(lam4)) lamfun(la04,la,la4) else lam3
    }
  } 
  new.par <- start
  PTOT1 <- if (is.null(pnl.type)||(pnl.type %in% c("lassoV","lasso")))  
              matrix(0,totcol,totcol) else
              penalty.TERM(lam1,lam2,lam3,lam4,propo, xx1,xx2,xx12,totcol,
                           ncat1,ncat2,s1,s2,s3,s4,pnl.type,new.par,1)
  
  VV <- if (constraints) penalty.INEQUALITY(lamv1,lamv2,ncat1,ncat2,
                           neweta[,2:ncat12_1],Z[margins,],constr,m,type) else 0
  P <- PTOT1 + VV
  
  lams.int <- NULL
   if (!smoother) {
     Pmat <- PENmat <- PPmat <- pm <- which.fixed.lams <- NULL
     ok <- rep(TRUE,3)   
   } else {
     for (i in 1:3) if (!is.null(fsmooth[[i]])) lamsa[[i]] <- la[names(la) %in% names(lamsa[[i]])]
     if (select.lams){
        if (sum(na.omit(unlist(lc)))>0){
          lams.start <- invert.list2(lapply(unlist(lamsa),function(x)x),TRUE)
        } else lams.start <- invert.list2(lamsa,TRUE)
     } else { 
       pslist <- make.Psmooth2(fsmooth,lamsa,prop.smooth,ncl,iB)     
       Pmat <- lapply(pslist,function(x)t(Mdiag(x)))
       PPmat <- mapply(function(x,y)tcrossprod(x[y,]),Pmat,ok,SIMPLIFY=FALSE)
       PENmat <- Mdiag(list(P,Mdiag(PPmat)))
     }
  }    

  p <- matrix(prob,ncol=ncat,nrow=m,byrow=TRUE)
  lwY <- log(ifelse(wY==0,1,wY))
  rwY <- rowSums(wY)
  ll <- dlogMult(wY, p,rwY)
  llp <- ll 
  iter <- 0 
  etasmooth <- 0
  if (is.null(l)) l <- 1
  if (smoother){  
    if (!((is.null(pnl.type)||(pnl.type!="equal"))))
      one.smooth <- lapply(1:length(wh.eq),function(i)crossprod(lam1*one.smooth[[i]]))
  }
 
  while((NOTCONV)&(iter<maxit)){ 
    iter <- iter + 1     
    old.par <- new.par
    llp.old <- llp
    ll.old <- ll
    eta2 <- subset(neweta,select=-zerovect)
    if (!bin){ 
      swp_r <- matrix(sqrt(wei_sum_row_Y/p),ncol=ncol(neweta))
      
      W2 <- lapply(1:m,function(i)crossprod(swp_r[i,]*GRAD.INV[,,i]))
    } else {
      Mp <- if ((iter>1)&(plackett==TRUE)) attributes(GRAD)$"Mp" else tcrossprod(p,M)
      W2 <- MIM2x2(Mp,wei_sum_row_Y)
    }
    z <- lapply(1:m,function(i) eta2[i,] + GRAD[-1,,i]%*%propY[i,])
    XtW2 <- lapply(1:m,function(i) Xt[,,i]%*%W2[[i]])
    IM <- Reduce('+',lapply(1:m,function(i) XtW2[[i]]%*%X[,,i]))
    ZWz <- Reduce('+',mapply("%*%",XtW2,z,SIMPLIFY=F)) 
    VV <- if (constraints) 
      penalty.INEQUALITY(lamv1,lamv2,ncat1,ncat2,neweta[,2:ncat12_1],
                         Z[margins,],constr,m,type) else 0
    if ((!is.null(pnl.type))&& pnl.type %in% c("lassoV","lasso")) 
       PTOT1 <- penalty.TERM(lam1,lam2,lam3,lam4,propo, xx1,xx2,xx12,
                   totcol,ncat1,ncat2,s1,s2,s3,s4,pnl.type,new.par,iter)
     P <- PTOT1 + VV  
    IMp <- IM + P
    half.step <- TRUE
    z <- unlist(z)
    #browser()
    if (smoother){
      BW <- lapply(Bs,function(x) mapply(crossprod,x,W2,SIMPLIFY=FALSE))      
      BWB <- mapply(function(x,y) Reduce('+',mapply("%*%",x,y,SIMPLIFY=FALSE)), 
                    BW,Bs,SIMPLIFY=F)    
      etalist <- lapply(iB,function(x) matrix(0,length(z),1))
      prodBW <- mapply(function(x) do.call("cbind",x),BW,SIMPLIFY=FALSE)
      if (select.lams){
        pm <- make.ppmat(BWB,ppmat,cwe,dflist,wh.eq,which.fixed.lams,lams.int,
                         lams.start,pgtol.df,lmm.df,factr.df,parscale.df,iter)
        lams.int <- pm$lams.int
        PPmat <- pm$PPmat
      }
      PBWB <- mapply('+',BWB,PPmat,SIMPLIFY=F)
    }  
    out <- NULL
    while (half.step){
      new.par <- (1 - l) * old.par + l * solve(IMp,ZWz)
      etas <- Z%*%new.par
      w.res <- z - etas
      out <- additive.fit(smoother,Bmat,new.smooth,PBWB,one.smooth,prodBW,
                          etalist,w.res,max.backfitting,rss.tol,maxB,m,ncat)
      
      neweta <- cbind(zerovect,matrix(etas,m,ncat-1,byrow=TRUE) + out$etasmooth)
      
      if (plackett){
        p <- dplackett.pblm(neweta,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,type)
        GRAD <- find.score2(C,M,p)
      } else { 
        ris <- mapply(function(a1,a2) t(find.pi3(a1,a2,M,C,Ct,acc2,maxit2,rep(1,length(neweta)))),
                      split(neweta,1:m), split(log(p),1:m), SIMPLIFY=FALSE)
        p <- do.call("rbind",ris)
        GRAD <- array(sapply(ris,function(x) attributes(x)$der1),c(ncat,ncat,m))         
      }
      
      if (!bin){ 
          GRAD.INV <- try(apply(GRAD,3,solve))
            if (class(GRAD.INV)!="try-error")
               GRAD.INV <- array(GRAD.INV,c(ncat,ncat,m))[,-1,,drop=FALSE] else {
               stop(paste("execution failed at iteration ",iter))
            }
      } 
      if (any(p<0)) {
        l <- l/2 
        half.step=TRUE
      } else half.step=FALSE
      if ((restore.l)&&((iter>1)&(l < min.step.l))) l <- 1
      if (l < min.step.l) 
        stop("step length has achieved a too small value.\n")
    }
    ll <- dlogMult(wY, p, rwY)
    Yp <- rwY*p
    dev <- 2*sum(wY*(lwY-log(ifelse(Yp==0,1,Yp)))) 
    penl <- 0
    if (!smoother)  penl <- drop(.5*crossprod(new.par,P%*%new.par)) else {
      newsm <- unlist(out$new.smooth)
      param <- c(new.par,newsm)
      if (select.lams) PENmat <- Mdiag(list(P,Mdiag(PPmat)))
      penl <- drop(.5*crossprod(param,PENmat%*%param))
    }
    llp <- ll - penl
    devp <- dev + 2*penl
    if (conv.crit=="dev"){if (ll<ll.old) l <- l/2} else { if (llp<llp.old) l <- l/2}
    if ((restore.l)&&((iter>1)&(l < min.step.l))) l <- 1
    if (l < min.step.l) stop("step length has achieved a too small value")
    toll <- switch(conv.crit,
                   "dev"=abs(ll-ll.old)/abs(ll.old),
                   "pdev"=abs(llp-llp.old)/abs(llp.old)
    )                
    NOTCONV <- (toll > acc)
    if (iter==1) NOTCONV <- TRUE    
    if ((verbose)&(!auto.select)){
      cat("iteration = ",iter, if (conv.crit=="dev") "  logLik = " else 
        "  penalized logLik = ",round(ll,digits=7),           
        "  accuracy = ",format(max(toll),scientific=TRUE,digits=4),
        " step.l= ",l,"\n")
    }
  }
  W <- unlist(lapply(W2,diag)) 
  WZP <- crossprod(sqrt(W)*Z)+P
  Ac <- chol(t(solve(WZP,diag(ncol(WZP)))))
  Zac <- tcrossprod(Z,Ac)
  df.fix.vect <- colSums(Zac*Zac*W)
  df.fix <- sum(df.fix.vect)
  if (smoother){
    df.smooth <- c()
    inv.fsmooth <- invert.list(fsmooth,maxB)
    for (j in 1:maxB){
      a <- 1:length(wh.eq[[j]])
      df.smooth.vect <-  ifelse(attributes(inv.fsmooth[[j]][[1]])$"smooth.name" %in% c("ridge","lasso"),0,2)*(1:3 %in% a)
      nonsmooth.npar <- sum(df.smooth.vect)
      df.smooth[j] <- sum(diag(solve(PBWB[[j]]+one.smooth[[j]])%*%BWB[[j]])) - nonsmooth.npar       
    }  
  } else df.smooth <- 0
  df.tot <- df.fix + sum(df.smooth)
  gAIC <- - 2 * ll + gaic.m * df.tot
  if (auto.select){
    cat(paste("GAIC(",gaic.m,"):",round(gAIC,digits=4),sep=""),
        if (!is.null(lam1)) " lam1=",lam1,
        if (!is.null(lam2)) " lam2=",lam2,
        if (!is.null(lam3)) " lam3=",lam3,
        if ((!is.null(lam4))&((!is.null(pnl.type))||(pnl.type=="ARC2"))) " lam4=",lam4,
        if (!is.null(lamsa[[1]])) " lambda[[1]]=",lamsa[[1]],
        if (!is.null(lamsa[[2]])) " lambda[[2]]=",lamsa[[2]],
        if (!is.null(lamsa[[3]])) " lambda[[3]]=",lamsa[[3]],"\n")                                                   
  } 
  attr(gAIC,"GAIC") <- gAIC
  attr(gAIC,"df.fix") <- df.fix    
  attr(gAIC,"df.smooth") <- df.smooth
  attr(gAIC,"df.tot") <- df.tot     
  attr(gAIC,"p") <- p
  attr(gAIC,"P") <- P
  attr(gAIC,"IM") <- IM
  attr(gAIC,"IMp") <- IMp
  attr(gAIC,"new.par") <- new.par
  attr(gAIC,"neweta") <- neweta
  attr(gAIC,"new.smooth") <- out$new.smooth
  attr(gAIC,"etasmooth") <- out$etasmooth
  attr(gAIC,"etalist") <- out$etalist
  attr(gAIC,"NOTCONV") <- NOTCONV
  attr(gAIC,"ZWz") <- ZWz 
  attr(gAIC,"W2") <- sapply(W2,function(x)x[lower.tri(x,TRUE)])
  attr(gAIC,"z") <- z      
  attr(gAIC,"Bmat") <- Bmat
  attr(gAIC,"ok") <- ok
  attr(gAIC,"w.res") <- out$w.res
  attr(gAIC,"BWB") <- BWB 
  attr(gAIC,"PBWB") <- PBWB 
  attr(gAIC,"lamsa") <- lamsa
  attr(gAIC,"lam1") <- lam1
  attr(gAIC,"lam2") <- lam2
  attr(gAIC,"lam3") <- lam3
  attr(gAIC,"lam4") <- lam4
  attr(gAIC,"Y") <- Y 
  attr(gAIC,"wh.eq") <- wh.eq
  attr(gAIC,"maxB") <- maxB
  attr(gAIC,"spec.name") <- spec.name
  attr(gAIC,"PPmat") <- PPmat
  attr(gAIC,"tol") <- toll
  attr(gAIC,"iter") <- iter
  attr(gAIC,"llp") <- llp
  attr(gAIC,"ll") <- ll
  attr(gAIC,"dev") <- dev
  attr(gAIC,"devp") <- devp
  attr(gAIC,"one.smooth") <- one.smooth
  attr(gAIC,"df.fix.vect") <- df.fix.vect
  attr(gAIC,"lambda") <- lams.int
  return(gAIC)  
}


pblm.prop <- function(prop1=NULL, prop2=NULL, prop12=NULL){
  list("prop1"=prop1, "prop2"=prop2, "prop12"=prop12)
}


pblm.penalty <- function(pnl.type=c("none","ARC1","ARC2","ridge","lasso","lassoV","equal"), 
                         lam1=NULL, lam2=NULL, lam3=NULL, lam4=NULL, s1=NULL, s2=NULL, 
                         s3=NULL, s4=NULL, min.lam.fix=0.1, constraints=FALSE, 
                         lamv1=1e10, lamv2=1e10){
  
  list("pnl.type"=match.arg(pnl.type), "lam1"=lam1, "lam2"=lam2, "lam3"=lam3, "lam4"=lam4, 
       "s1"=s1, "s2"=s2, "s3"=s3, "s4"=s4, "min.lam.fix"=min.lam.fix, "constraints"=constraints, 
       "lamv1"=lamv1, "lamv2"=lamv2)
}

pblm.control <- function(maxit=30, maxit2=200, acc=1e-07, acc2=1e-06,zero.adj=1e-06, 
                         l=NULL, restore.l=FALSE, min.step.l=1e-4, auto.select=FALSE,
                         gaic.m=2, rss.tol=1e-06, max.backfitting=10, pgtol.df=1e-2, 
                         factr.df=1e7, lmm.df=5,parscale.df=1,max.gaic.iter=500, 
                         pgtol.gaic=1e-05,grad.tol=1e-07,
                         factr.gaic=1e7, lmm.gaic=5, parscale=1,conv.crit=c("dev","pdev")){
  
  list("maxit"=maxit, "maxit2"=maxit2, "acc"=acc, "acc2"=acc2, "zero.adj"=zero.adj, 
       "l"=l, "restore.l"=restore.l, "min.step.l"=min.step.l, "auto.select"=auto.select,
       "gaic.m"=gaic.m, "rss.tol"=rss.tol, "max.backfitting"= max.backfitting, 
       "pgtol.df"=pgtol.df, "factr.df"=factr.df, "lmm.df"=lmm.df,"parscale.df"=parscale.df, 
       "max.gaic.iter"=max.gaic.iter, "pgtol.gaic"=pgtol.gaic, "factr.gaic"=factr.gaic,
       "lmm.gaic"=lmm.gaic, "conv.crit"=match.arg(conv.crit),"parscale"=parscale,
       "grad.tol"=grad.tol) 
  
}

getSmo <- function(x,which.eq = 1, which.var = 1){
  if (class(x)!="pblm") stop("This is not a pblm object")
  if (x$any.smoother){
    index <- unlist(x$wh.eq2[[which.var]]) #new
    if (is.na(index[which.eq])) 
      stop("no smoothers from the selected equation/variable")
    na.num <- sum(is.na(index[1:which.eq]))
    na.num <- ifelse(which.eq > na.num,na.num,0)
    index <- na.omit(index)
    cindex <- c(0,cumsum(index))
    wh.smo <- rep(FALSE,sum(index))
    wh.smo[(cindex[which.eq-na.num]+1):cindex[which.eq-na.num+1]] <- TRUE
    smo <- unlist(x$beta.smooth[which.var])[wh.smo]
    attr(smo,"lambda") <- unlist(x$lambda[which.var])[which.eq-na.num]
    attr(smo,"df") <- edf.pblm(x,which.var,which.eq)
    smo
  } else NULL
}


plot.pblm <- function(x, which.eq = 1:3, which.var = 1:x$maxNpred, add.bands = TRUE,
                      col.bands = "lightblue", col.line = "black", type = "l", 
                      dashed.bands = FALSE, pause = FALSE, ylim, xlim, ylab, xlab, 
                      main, ...) 
{
  if (class(x) != "pblm") 
    stop("This is not a pblm object")
  YLIM <- ifelse(missing(ylim), TRUE, FALSE)
  XLIM <- ifelse(missing(xlim), TRUE, FALSE)
  main.text <- if (missing(main)) 
    c("Marginal 1","Marginal 2","Association") else main
  mo.call <- x$call 
  mo.call$fit <- FALSE
  mo.call$ncat1 <- x$ncat1
  mo.call$ncat2 <- x$ncat2
  mo <- eval(mo.call) 
  pred <- predict(x, type="terms", se.fit=add.bands)
  se.fit <- if (add.bands) pred$se.fit else NULL
  pred <- if (add.bands) pred$fitted.values else pred
  CILower <- if (add.bands) pred-2*se.fit else NULL
  CIUpper <- if (add.bands) pred+2*se.fit else NULL
  co <- a <- list()
  k<-0
  spm <- x$spec.smooth
  nc <- x$nc
  simplifyList <- function(l){

    k <- 0
    new.list <- list()
    for (i in 1:length(l)){
      if (is.list(l[[i]])){
        for (j in 1:length(l[[i]])){
          k <- k + 1
          new.list[[k]] <- l[[i]][[j]]
        } 
      } else  {
        k <- k + 1
        new.list[[k]] <- l[[i]]
      }
    }
    
    return(new.list)
  }
  x.by.var <- simplifyList(x$by.var)
  kvet <- as.vector(t(x$which.term))
  kvet[kvet] <- 1:sum(x$which.term)
  kmat <- matrix(kvet,nrow=3,ncol=x$maxNpred,byrow=T)
  for (i in which.eq) {
    nlevint <- mo$nlev[[i]]
    which.xnames <- c("(Intercept)", attr(terms(mo)[[i]],"term.labels"))
    aaa <- factor(nlevint, labels = which.xnames)
    aaa2 <- factor(mo$nlev0[[i]], labels = which.xnames)
    asgn <- split(order(mo$nlev0[[i]]), aaa2)
    asgn$"(Intercept)" <- NULL
    nterms <- length(asgn)
    if (nterms>0) {
      k2 <- 0
      for (j in which.var) {
       # if ((i==2)&(j==1)) browser()
        if (j > length(attr(terms(x)[[i]],"term.labels"))) break
         k <- kmat[i,j]
        if (!(substr(which.xnames[j+1],1,3) %in% c("pbs","pb(","pvc"))){
          if (x$is.fac[[i]][j]){
            xm <- cbind(0,as.matrix(mo$Xlist[[i]][,asgn[[j]]]))
            colnames(xm) <- x$fac.lev[[i]][[j]]
            co[[k]] <- as.factor(apply(xm,1,function(foo){return(colnames(xm)[which.max(foo)])}))
          } else co[[k]] <- mo$Xlist[[i]][,asgn[[j]],drop=F] 
        } else { 
          k2 <- which(x$spec.name[[i]] %in% which.xnames[j+1])
          if ((!is.null(x$spec.name[[i]])) && 
              (!is.na(spm[[i]][k2]))&&
              (!(substr(spm[[i]][k2],1,5) %in% c("ridge","lasso")))) {
            if (is.list(x$fsmooth[[i]])){
              co[[k]] <- if (attributes(x$fsmooth[[i]][[k2]])$"smooth.name"=="pvc") 
                as.vector(attributes(x$fsmooth[[i]][[k2]])$"xorig") else 
                  as.matrix(x$fsmooth[[i]][[k2]])
            } else {
              co[[k]] <- if (attributes(x$fsmooth[[i]])$"smooth.name"=="pvc") 
                as.matrix(attributes(x$fsmooth[[i]])$"xorig") else 
                  as.matrix(x$fsmooth[[i]][,k2])          
            }
          }  
        }
        a[[k]] <- cbind(co[[k]],  pred[,k], CILower[,k], CIUpper[,k])
        a[[k]] <- a[[k]][order(a[[k]][,1]),]
        if (YLIM){ 
          ylim <- if (add.bands) 
            range(a[[k]][, 3:4]) + c(-0.5, 0.5)
          else range(a[[k]][, 2]) * c(1.1, 0.95)
        }
        if (substr(which.xnames[j+1],1,3) %in% c("pbs","pb(","pvc")){
          yl <- ifelse(missing(ylab),x$spec.smooth[[i]][k2],ylab) 
          spm2 <- spm[[i]][k2]
          wh.spm2 <- which(strsplit(spm2,"")[[1]] %in% c(")",","))-1
          st <- switch(substr(x$spec.smooth[[i]][k2],1,3),
                       "pbs"=substr(spm2,5,wh.spm2),
                       "pb("=substr(spm2,4,wh.spm2),
                       "pvc"=substr(spm2,5,wh.spm2))
          xl <- ifelse(missing(xlab),st,xlab)
        } else {
          if (!x$is.fac[[i]][j]){
            st <- unique(x$xnames[[i]])[j+1]#new
            yl <- ifelse(missing(ylab),paste("partial for",st),ylab)
            xl <- ifelse(missing(xlab),st,xlab)
          }
        }  
        if (x$is.fac[[i]][j]){
          ll <- x$fac.lev[[i]][[j]]
          st <- names(x$fac.lev[[i]])[j]
          yl <- ifelse(missing(ylab),paste("partial for",st),ylab)
          xl <- ifelse(missing(xlab),st,xlab) 
          xlims <- range(seq_along(ll)) + c(-0.5, 0.5)
          if (pause) readline(prompt="Press any key to continue")
          
          plot(1, 0, type = "n", ylim = ylim, xlim = xlims,
               ylab=yl,xlab=xl,main=main.text[i],col=col.line, xaxt = "n", 
               ...)
            axis(1, at = seq_along(ll), labels = ll, ...)
          for (j2 in unique(a[[k]][,1])){
            jf <- j2 + c(-0.4, 0.4)
            jj <- which(!duplicated(a[[k]][,1]))[j2]
            lines(jf, a[[k]][c(jj,jj),2],...)
            if (add.bands){ 
              xy <- a[[k]][c(jj,jj),1:2]
              NotZero <- xy[,2]!=0
              z <- xy[NotZero, 1]
              yLower <- a[[k]][c(jj,jj), 3]
              yUpper <- a[[k]][c(jj,jj), 4]
              if (dashed.bands){
                lines(jf, yLower, lty="dashed",col=col.bands,...)
                lines(jf, yUpper, lty="dashed",col=col.bands,...)
              } else {
                polygon(c(jf, jf[length(z):1]), c(yLower, yUpper[length(yUpper):1]), 
                        lwd = 2, col = col.bands, border = NA)
                lines(jf, a[[k]][c(jj,jj), 2],type=type,col=col.line,...)
              }
              
            }
          }
        } else {
          if (XLIM) xlim <- range(co[[k]])
          if ((substr(which.xnames[j+1],1,3) == "pvc")&(is.factor(x.by.var[[k]]))){
              xy <- a[[k]][,1:2]
              ord.lev <- x.by.var[[k]][order(co[[k]])]
              for (i.lab in levels(x.by.var[[k]])){
                NotZero <- ord.lev==i.lab
                if (pause) readline(prompt="Press any key to continue")
                plot(xy[NotZero,], type=type, ylim = ylim, xlim = xlim,
                     ylab=paste(yl,":",i.lab,sep=""),xlab=xl,main=main.text[i],
                     col=col.line,...)
                z <- a[[k]][NotZero, 1]
                if (add.bands) {
                  yLower <- a[[k]][NotZero, 3]
                  yUpper <- a[[k]][NotZero, 4]
                  if (dashed.bands){
                    lines(z, yLower, lty="dashed",col=col.bands,...)
                    lines(z, yUpper, lty="dashed",col=col.bands,...)
                  } else {
                    polygon(c(z, z[length(z):1]), c(yLower, yUpper[length(yUpper):1]), 
                            lwd = 2, col = col.bands, border = NA)
                    lines(z, a[[k]][NotZero, 2],type=type,col=col.line,...)
                  }
                }
                segments(z, (0.96 * (ylim[1] > 0) + 1.04 * (ylim[1] < 0)) * 
                           ylim[1], z, (0.96 * (ylim[1] > 0) + 1.04 * 
                                          (ylim[1] < 0)) * ylim[1] + 0.02 * diff(ylim))
              }
          }  else {
            xy <- a[[k]][,1:2]
            if (pause) readline(prompt="Press any key to continue")
            plot(xy, type=type, ylim = ylim, xlim = xlim,
                 ylab=yl,xlab=xl,main=main.text[i],col=col.line,...)
            z <- a[[k]][, 1]
            if (add.bands) {
              yLower <- a[[k]][, 3]
              yUpper <- a[[k]][, 4]
              if (dashed.bands){
                lines(z, yLower, lty="dashed",col=col.bands,...)
                lines(z, yUpper, lty="dashed",col=col.bands,...)
              } else {
                polygon(c(z, z[length(z):1]), c(yLower, yUpper[length(yUpper):1]),
                        lwd = 2, col = col.bands, border = NA)
                lines(z, a[[k]][, 2],type=type,col=col.line,...)
              }
            }
            segments(z, (0.96 * (ylim[1] > 0) + 1.04 * (ylim[1] < 0)) * 
                     ylim[1], z, (0.96 * (ylim[1] > 0) + 1.04 * 
                                    (ylim[1] < 0)) * ylim[1] + 0.02 * diff(ylim))
          } 
        }
        names(a)[k] <- colnames(pred)[k]
        colnames(a[[k]]) <- c(st,colnames(pred)[k],
                              if (add.bands) "C.I. lower",
                              if (add.bands) "C.I. upper")
      }
    }
  }
  return(invisible(a))
}



summary.pblm <- function(object,...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  z <- object
  val <- aic.one.model(z,z$gaic.m)
  ll <- z$ll
  llp <- z$llp
  GAIC <- val$AIC
  df.tot <- val$df
  df.smooth <- df.tot - z$df.fix 
  df.res <- if (any(z$weights!=1)) z$m*(z$ncat-1)-df.tot else z$n*(z$ncat-1)-df.tot
  AIC.m <- -2*(ll-df.tot)
  BIC.m <- -2*ll+log(z$n)*df.tot
  invF <- solve(z$IMp)
  se <- z$coef
  se[!is.na(se)] <- round(sqrt(diag(invF)),digits=7)
  z1 <- round(z$coef/se,digits=7)
  results <- data.frame("beta"=round(as.vector(z$coef),digits=7),"se"=se,"z"=as.vector(z1),
                        "p.value"=round(as.vector(2*pnorm(abs(z1),lower.tail=FALSE)),digits=5)) 
  rownames(results) <- names(z$coef)
  PR <- t(apply(resid(z,"pearson"),2,fivenum))
  colnames(PR) <- c("min","1q","median","3q","max")
  rval <- list("results"=results,"convergence"=z$convergence,"iter"=z$iter,
               "tol"=as.vector(z$tol),"logLik"=ll,"logLikp"=llp,"AIC"=AIC.m,"gAIC"=GAIC,
               "BIC"=BIC.m,"gaic.m"=z$gaic.m,"np1"=sum(z$xx1),"deviance"=deviance(z),
               "np2"=sum(z$xx2),"names"=z$ynames,"df.res"=df.res,"df.tot"=df.tot,
               "df.fix"=z$df.fix,"df.smooth"=df.smooth,"res.smooth"=val$res.smooth,
               "smoother"=z$any.smoother,"pnl.type"=z$pnl.type,"PR"=PR,"label1"=z$label1,
               "label2"=z$label2)
  rval$call <- object$call
  class(rval) <- "summary.pblm"
  rval
}





chisq.smooth.pblm <- function(obj,which.var,which.eq){
  if (class(obj)!="pblm") 
    stop("This is not a pblm object")
  if (!obj$any.smoother) 
    stop("No smoother has been used in the pblm model")
  wh.eq2 <- obj$wh.eq2[[which.var]]
  a <- (1:3)[!is.na(wh.eq2)]
  eq <- rep(a==which.eq, na.omit(wh.eq2))
  if (all(!eq)) return(NULL) else { 
    PBWB <- obj$PBWB[[which.var]]
    BWB <- obj$BWB[[which.var]]
    cPBWB <- chol(PBWB[eq,eq])
    beta <- obj$beta.smooth[[which.var]][eq]
    Cq <- drop(crossprod(cPBWB%*%beta))
    smo.eq <- obj$fsmooth[[which.eq]]
    if (is.list(smo.eq)) 
      smo.eq <- smo.eq[[which.var]]
    smo.name <- attributes(smo.eq)$"smooth.name"
    opt.rid  <- ifelse(smo.name %in% c("ridge","lasso"),0,2) 
    tr.df <- sum(diag(solve(PBWB)%*%BWB)[eq])
    df <- round(tr.df,digits=4) - opt.rid
    p.value <- if (df!=0) 1-pchisq(Cq,df) else NA
    data.frame("Chisq"=Cq,"df"=df, "p.value"=p.value)
  }
}



edf.pblm <- function(object,which.var=1,which.eq=1){
  if (class(object)!="pblm") stop("This is not a pblm object")
  if (!object$any.smoother) stop("No smoother has been used in the pblm model")
  a <- !sapply(object$wh.eq2[[which.var]],is.na)
  which.eq.is.smooth <- (1:3)[a]==which.eq
  eq <- rep(which.eq.is.smooth,object$wh.eq[[which.var]])
  if (all(!eq)) return(NULL) else{
    non.smooth.npar <- sum(2*which.eq.is.smooth)
    sum(diag(solve(object$PBWB[[which.var]])%*%object$BWB[[which.var]])[eq])-non.smooth.npar
  }
}


se.smooth.pblm <- function(object){
  if (class(object)!="pblm") stop("This is not a pblm object")
  if (!object$any.smoother) stop("No smoother has been used in the pblm model")
  invPBWB <- lapply(object$PBWB,solve)
  Bmat <- object$Bmat
  sa <- mapply(function(x,y) crossprod(x,y)%*%x,invPBWB,object$BWB,SIMPLIFY=FALSE)
   for (i in 1:length(object$Bmat)){
      for (j in 1:3){
        if (length(object$fsmooth[[j]])>=i){
          if (!is.null(attributes(object$fsmooth[[j]][[i]])$"Xor")){
            index <- unlist(object$wh.eq2[[i]]) 
            na.num <- sum(is.na(index[1:j]))
            na.num <- ifelse(j > na.num,na.num,0)
            index <- na.omit(index)
            cindex <- c(0,cumsum(index))
            col.indeces <- rep(FALSE,sum(index))
            col.indeces[(cindex[j-na.num]+1):cindex[j-na.num+1]] <- TRUE
            row.indeces <- object$respStr[,j]==1
            Xor <- attributes(object$fsmooth[[j]][[i]])$"Xor"
            Bmat[[i]][row.indeces,col.indeces] <- kronecker(Xor,rep(1,object$nc[j]))
          }
        }
      }
  }
  
  vco <- mapply(function(x,y) rowSums(x*(x%*%y)),Bmat,sa,SIMPLIFY=FALSE)
  seco <- lapply(vco,function(x)matrix(sqrt(x),,object$ncat-1,byrow=TRUE))
  smooth.inf <- mapply(function(x,y) matrix(x,object$m,object$ncat-1,byrow=TRUE)-2*y,
                       object$etalist,seco,SIMPLIFY=FALSE)
  smooth.sup <- mapply(function(x,y) matrix(x,object$m,object$ncat-1,byrow=TRUE)+2*y,
                       object$etalist,seco,SIMPLIFY=FALSE)
  list("se.smooth"=seco,"smooth.inf"=smooth.inf,"smooth.sup"=smooth.sup)
}


print.summary.pblm <- function(x,digits = max(3, getOption("digits") - 3),...){     
  cat("\n")
  cat("********* Bivariate marginal logistic regression model of class pblm **********","\n")
  if (!is.null(x$call))
    cat("Call:","\n")
  cat(deparse(x$call), fill=90)
  if (!inherits(x,"summary.pblm")) stop("object is not of class summary.pblm")
  cat("\n")
  cat(" ------------------------------------------------------------------------------","\n")
  cat("                             log-likelihood :",round(x$logLik,digits=3),"\n")
  cat("                          Residual Deviance : ",round(x$deviance,digits=3),"\n")
  cat("                                        AIC : ",round(x$AIC,digits=3),"\n")
  cat(paste("                                    GAIC(",x$gaic.m,") :  ",round(x$gAIC,digits=3),sep=""),"\n")
  cat("                                        BIC : ",round(x$BIC,digits=3),"\n")
  cat("                Residual degrees of freedom : ",round(x$df.res,digits=3),"\n")
  cat("            Degrees of freedom of the model : ",round(x$df.tot,digits=3),"\n")
if (x$smoother){
  cat("  Degrees of freedom of the parametric part : ",x$df.fix,"\n")
  cat("    Degrees of freedom of the additive part : ",round(x$df.smooth,digits=3),"\n")
}
  cat(" ------------------------------------------------------------------------------","\n")
  sig <- function(z){ if (is.na(z)) "   " else if (z<0.001) "***" else if (z<0.01) "** " else if (z<0.05) "*  "  else if (z<0.1) ".  " else "   "}
  if (length(x$results[[1]])) {
    cat("\n")
    cat("Pearson Residuals:\n")
    print(x$PR)
    cat("\n")
    cat("Coefficients:\n")
    results <- na.omit(x$results)
    if (!is.null(x$names)){
      names <- as.vector(if (length(x$names)>2) 
               c(as.character(x$label1),as.character(x$label2)) else x$names)
      sig.1 <- sapply(results$p.value,sig)
      est.1 <- cbind(format(results, digits = digits),sig.1)
      colnames(est.1)[ncol(est.1)] <- ""
      cat("Marginal logit 1: ",names[1],"\n")
      #browser()
      print(format(est.1[1:x$np1,], digits = digits))
      cat("\n")
      cat("Marginal logit 2: ",names[2],"\n")
      print(format(est.1[(x$np1+1):(x$np1+x$np2),], digits = digits))
      cat("\n")
      cat("Association: ",names[1]," vs ",names[2],"\n")
      print(format(est.1[(x$np1+x$np2+1):nrow(results),], digits = digits))
      cat("\n")

    } else print(results[[1]])
  } else cat("No coefficients","\n")
  cat("-------------------------------------------------------------------------------","\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
  cat("\n")
  if (x$smoother){
    cat("Approximate significance for the smoothers: ","\n\n")
    sig.smooth <- sapply(x$res.smooth$p.value,sig)
    est.smooth <- cbind(format(x$res.smooth, digits = digits),sig.smooth)
    colnames(est.smooth)[ncol(est.smooth)] <- ""
    print(est.smooth)
    
    cat("\n")
    cat(" ------------------------------------------------------------------------------","\n")
  }  
  if (x$convergence) cat("The algorithm converged after ",x$iter," iterations","\n") else 
    cat("The algorithm did not converge after ",x$iter," iterations","\n")
  cat("\n")  
  if (x$smoother){ 
  cat("-------------------------------------------------------------------------------","\n")
  cat("NOTE: In presence of additive smoothing terms:\n")
  cat(" - Standard errors for the linear terms might be not accurate\n")
  if ((!is.null(x$pnl.type))&&((x$pnl.type=="equal")))
  cat(" - Degrees of freedom of the model might be inexact\n")  
  }
  cat("*******************************************************************************","\n")
  invisible(x)
}


print.pblm <-function(x,digits = max(3, getOption("digits") - 3),...){  
  if (!inherits(x,"pblm")) stop("x is not of class pblm")
  cat("Bivariate logistic regression model:\n")
  if (!is.null(x$call))
    cat("\nCall: ", deparse(x$call), "\n",fill=TRUE)
  if (length(x$coef)) {
    cat("Coefficients:\n")
    print.default(format(x$coef, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}


AIC.pblm<- function(object,...,k = 2){
  if (length(list(...))) {
    object <- list(object, ...)
    ispblm <- unlist(lapply(object, function(x)(inherits(x,"pblm"))))
    if (!any(ispblm)) 
      stop("some of the objects are not pblm")
    aic <- subset(as.data.frame(t(sapply(object, function(z) aic.one.model(z,gaic.m=k)))),
                  select=c(df,AIC))
    Call <- match.call()
    Call$k <- NULL
    row.names(aic) <- as.character(Call[-1])
    aic <- aic[order(unlist(aic$AIC)), ]
    aic
  } else {
    if (!inherits(object,"pblm")) stop("object is not of class pblm")
    val <- aic.one.model(object,k)
    aic <- val$AIC
    print.default(paste(round(aic,digits=4)," with ",round(val$df,digits=4)," df",sep=""))
    attr(aic,"k") <- k
    attr(aic,"ll") <- object$ll
    attr(aic,"df.tot") <- val$df
    attr(aic,"df.fix") <- object$df.fix 
    attr(aic,"df.smooth") <- val$df-object$df.fix
    return(invisible(aic))  
  }
}

aic.one.model<- function(z,gaic.m=2){
   if (!z$any.smoother) {
     res.smooth <- NULL
     df.smooth <- 0
   } else {
     k <- 0
     res.smooth <- list()
     res.smooth.names <- c()
     for (ik in 1:z$maxB){
       for (i in 1:3){
         if (ik<=length((z$fsmooth[[i]]))){
           Chisq.test <- chisq.smooth.pblm(z,ik,i)
           if (!is.null(Chisq.test)){
             k <- k + 1
             res.smooth[[k]] <- as.data.frame(Chisq.test) 
             Chisq.test <- NULL
             res.smooth.names[k] <- z$spec.smooth[[i]][ik]
           }
         }
       }
     }

     names(res.smooth) <- res.smooth.names
     res.smooth <- do.call("rbind",res.smooth)
   }
   df.smooth <- if (z$any.smoother) sum(res.smooth$df) else 0
   df.tot <- z$df.fix + df.smooth
   GAIC <- - 2 * z$ll + gaic.m * df.tot
   list("df"=df.tot,"AIC"=GAIC,"res.smooth"=res.smooth) 
 }



logLik.pblm<- function(object,penalized=FALSE,...){
  z <- object
  if (!inherits(z,"pblm")) stop("object is not of class pblm")
  if (penalized) z$llp else z$ll
}

deviance.pblm <- function(object,penalized=FALSE,...){
  z <- object
  if (!inherits(z,"pblm")) stop("object is not of class pblm")
  dev <- if (penalized) z$devp else z$dev
  dev
}

fitted.pblm <- function(object,...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  grid <- expand.grid(1:object$ncat2,1:object$ncat1)[,2:1]
  index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
  fit <- as.matrix(object$p) 
  if (length(colnames(fit)>1)) colnames(fit) <- paste("p",index.lab,sep="")
  as.data.frame(fit)
}


predict.pblm <- function(object, newdata, type=c("link","response","terms","count"), 
                         se.fit=FALSE, digits= max(6, getOption("digits") - 3),...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  type <- match.arg(type)
  nc <- object$nc
  if (missing(newdata)){
    grid <- expand.grid(1:object$ncat2,1:object$ncat1)[,2:1]
    index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
    if (type=="response"){
      fit <- fitted(object)
      if (se.fit) warning("'se.fit=TRUE' is not allowed when 'type=response'")
      return(round(as.data.frame(fit),digits=digits))
    }      
    if (type=="count"){
      fit <- fitted(object)*rowSums(object$Y)
      if (length(colnames(fit))>1) colnames(fit)<- paste("n",index.lab,sep="") 
      if (se.fit) warning("'se.fit=TRUE' is not allowed when 'type=count'")
      return(round(as.data.frame(fit),digits=digits))
    }
    if (type=="link"){
      fit <- object$eta[,-1,drop=F]
      grid <- expand.grid(1:nc[2],1:nc[1])[,2:1]
      index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
      colnames(fit) <- c(paste("eta1:",1:nc[1],sep=""),
                         paste("eta2:",1:nc[2],sep=""),
                         paste("eta12:",index.lab ,sep=""))
    
      if (se.fit){ 
        o <- rowSums((object$x%*%vcov(object))*object$x)
        stder <- matrix(sqrt(o),nrow=nrow(fit),ncol=ncol(fit),byrow=TRUE)
        colnames(stder) <- colnames(fit)
        return(list("fit"=round(as.data.frame(fit),digits=digits),
                    "se.fit"=round(stder,digits=digits)))
      } else return(round(as.data.frame(fit),digits=digits))    
    } 
    if (type=="terms"){
    
      ncat1 <- object$ncat1 
      ncat2 <- object$ncat2
      ncat <- object$ncat
      if (length(object$coef) <= sum(nc)){
        warning("no terms to predict in the model")
        return(NULL)
      } 
      ox <- object$xnames
      oo <- object$coef
      os <- object$spec.name
      ow <- object$which.ns
      mo.call <- object$call 
      mo.call$fit <- FALSE
      mo.call$ncat1 <- object$ncat1
      mo.call$ncat2 <- object$ncat2
      mo <- eval(mo.call)
      maxnlev <- sum(sapply(mo$nlev,max))
      eta.x <- matrix(,object$m,maxnlev)
      which.coef <- c()
      lox <- lapply(ox,length)
      cumint <- c(0,cumsum(lox[-3]))
      cumcoef <- cumsum(lox)
      cumnc <- c(0,nc[1],sum(nc[1:2]))
      which.coef.list <- coef.list <- list()
      b.smooth <- vector("list",length(object$beta.smooth))
      if (se.fit) {
        fit.se <- matrix(0,nrow=mo$m,ncol=maxnlev )
        V <- vcov(object)
        sesmo <- if (object$any.smoother) se.smooth.pblm(object)$se.smooth else NULL
      }
      for (j in 1:length(ox)){
        which.coef <- c( which.coef,ox[[j]]%in%ox[[j]][-(1:nc[j])])
        which.coef.list[[j]] <- ox[[j]]%in%ox[[j]][-(1:nc[j])]
        coef.list[[j]] <- oo[(cumint[j]+1):cumcoef[j]]
      }
      lox2 <- lapply(coef.list,function(x)length(na.omit(x)))
      cumint2 <- c(0,cumsum(lox2[-3]))
      cumcoef2 <- cumsum(lox2)
      mar <- c("mar1","mar2","ass12")
      xnames <- c()
      k <- 0
      subterm=list()
      nterms <- termlab <- list()
      for (j in 1:length(ox)){
        nlevint <- mo$nlev[[j]]
        termlab[[j]] <- attr(terms(mo)[[j]],"term.labels")
        which.xnames <- c("(Intercept)", termlab[[j]])
        aaa <- factor(nlevint, labels = which.xnames)
        aaa2 <- factor(mo$nlev0[[j]], labels = which.xnames)
        asgn <- split(order(nlevint)-1, aaa) 
        asgn$"(Intercept)" <- NULL
        asgn2 <- split(order(mo$nlev0[[j]]), aaa2)
        asgn2$"(Intercept)" <- NULL
        nterms <- length(asgn)
        i1 <- 0
        k2 <- 0
        if (nterms>0) {
          for (i in 1:nterms){
           partial.eta <- 0
           k <- k + 1
           noterms <- !(aaa2%in%levels(aaa)[i+1]) 
           notcoe <- coef.list[[j]][-(1:nc[j])][noterms[-1]]
           notcoe[is.na(notcoe)]<-0 
           coe <- coef.list[[j]][-(1:nc[j])][!noterms[-1]]
           coe[is.na(coe)]<-0 
           x <- mo$Xlist[[j]][,!noterms,drop=TRUE] 
           notx <- if (length(notcoe)>0) notcoe else 0
           mean.x <- 0
           if (sum(noterms[-1])>0){
             Xno <- mo$Xlist[[j]][,noterms,drop=T][,-1]
             mean.x <- if (sum(noterms[-1])>1) 
                           mean(Xno%*%notx)  else                                                                                    
                           mean(notx*Xno)    
           }
           subterm[[k]] <- substr(termlab[[j]][i],1,3)
           if(length(x)>0){
             if (subterm[[k]] == "pvc")
               partial.eta <- mean.x else
               partial.eta <- if (ncol(as.matrix(x))>1) x%*%coe + mean.x else 
                                                        x*coe + mean.x
           } 
           mean.eta <- mean(object$eta[,cumnc[j]+2,drop=F]) 
           is.smooth.x <- ((!is.null(os))&&(length(os)>=j)&&(!is.null(os[[j]]))&&
                          (levels(aaa)[i+1] %in% os[[j]]))
           x.eta <- 0
           if (is.smooth.x){
             i1 <- i1 +1
             if (subterm[[k]] == "pvc"){
               
               dimN <- attr(mo$Xlist[[j]],"dimnames")[[2]]  
               nby  <- sum(dimN %in%  termlab[[j]][i])
               if (is.list(object$by.var[[j]]))
                 fact <- (is.factor(object$by.var[[j]][[i]])) else
                 fact <- (is.factor(object$by.var[[j]])) 
              
               if (fact){
                 for (i2 in i1:(i1+nby-1))
                   x.eta <- x.eta + drop(matrix(object$etalist[[i2]],mo$m,byrow=T)[,(cumnc+1)[j]])
               } else {
                   indeces <- which(levels(aaa)[i+1]==attr(object$fsmooth[[j]],"name"))
                if (length(indeces)>0){
                   x.eta <- attr(object$fsmooth[[j]][[indeces]],"Xor")%*%getSmo(object,j,indeces)
                  } else
                     x.eta <- attr(object$fsmooth[[j]][[1]],"Xor")%*%getSmo(object,j,1) #va indicizzato
               }
             } else {
               x.eta <- drop(matrix(object$etalist[[i1]],mo$m,byrow=T)[,(cumnc+1)[j]])
             }
           } 
          
           eta.x[,k] <- oo[1 + cumint[j]] + partial.eta - mean.eta + x.eta
           if (se.fit){
             if (length(colnames(mo$X[[j]]))>1){
               if (!is.smooth.x){
                 Xm <- sweep(mo$X[[j]], 2L, colMeans(mo$X[[j]]), check.margin=F)
                 ji <- cumint2[j] + asgn[[i]]
                 fit.se[,k] <- rowSums((Xm[,asgn2[[i]],drop=F]%*%
                                        sqrtM(V[ji,ji,drop=F],0.5))^2)
               } else {
                  k2 <- k2 + 1*((ow[[j]])[1])
                 if (subterm[[k]] == "pvc"){
                   dimN <- attr(mo$Xlist[[j]],"dimnames")[[2]]  
                   nby  <- sum(dimN %in%  termlab[[j]][i])
                   for (k3 in k2 : (k2 + nby - 1))
                     fit.se[,k] <- fit.se[,k] + sesmo[[k3]][,(cumnc+1)[j]]^2
                 } else {   
                   fit.se[,k] <- sesmo[[k2]][,(cumnc+1)[j]]^2
                 }   
               }
             }
           } 
         }
       }
       llab <- length(termlab[[j]])   
       xnames <- c(xnames,if (llab>0) paste(mar[j],":", termlab[[j]],sep=""))
       }
    if (se.fit){
        if (ncol(eta.x)>0) {
          fit <- list()
          fit$fitted.values <- eta.x
          colnames(fit$fitted.values) <- xnames
          fit$se.fit <- sqrt(fit.se)
          colnames(fit$se.fit) <- xnames
          return(fit)
        } else {
          fit <- object$eta[,-1, drop=F]
          rownames(fit) <- 1:nrow(fit)
          colnames(fit) <- xnames 
          return(fit)
        }
      } else {  
        if (ncol(eta.x)>0) fit <- eta.x else {
          fit <- object$eta[,-1, drop=F]
          rownames(fit) <- 1:nrow(fit)
          return(fit)
        }
        colnames(fit) <- xnames
      }
      fit <- round(as.matrix(fit),digits=digits)
      rownames(fit) <- 1:nrow(fit)
      attr(fit, "constant") <- colMeans(object$eta[,-1])
      if (ncat1>2)
        attr(fit, "mar1Contrasts") <- coef.list[[1]][2:nc[1]]-coef.list[[1]][1]
      if (ncat2>2)   
        attr(fit, "mar2Contrasts") <- coef.list[[2]][2:nc[2]]-coef.list[[2]][1]
      if ((ncat1>2)|(ncat2>2)){    
        attr(fit, "ass12Contrasts") <- coef.list[[3]][2:nc[3]]-coef.list[[3]][1]
      warning('predictions are with respect the first intercept of mar1, mar2 and
               ass12, respectively. If you want category-dependent predictions, 
               add the corresponding contrasts')
      }
      return(fit)
    } 
  } else {
    if (se.fit) warning("'se.fit=TRUE' is not allowed when 'newdata' is specified")      
    ncat1 <- object$ncat1 
    ncat2 <- object$ncat2
    ncat <- object$ncat
    if (length(object$coef) <= sum(nc)){
      warning("no terms to predict in the model")
      return(NULL)
    } 
    ox <- object$xnames
    oo <- object$coef
    os <- object$spec.name
    ow <- object$which.ns
    ncat1_1 <- ncat1-1 
    ncat2_1 <- ncat2-1 
    ncat_1 <- ncat-1 
    ncat12 <- ncat1+ncat2
    ncat12_1 <- ncat12-1
    ncat_1_2 <- ncat1_1*ncat2_1      
    mo.call <- object$call 
    mo.call$data <- newdata
    mo.call$fit <- FALSE
    mo.call$ncat1 <- object$ncat1
    mo.call$ncat2 <- object$ncat2
    fo1 <- terms.formula(mo.call$fo1)
    fo1[[2]] <- NULL
    mo.call$fo1 <- formula(fo1)
    mo <- eval(mo.call)
    if (object$any.smoother) 
      smo.list <- lapply(1:object$maxB,function(i) mo$Bmat[[i]]%*%object$beta.smooth[[i]])
    eta <- mo$Z%*%object$coef + 
           if (object$any.smoother) Reduce("+",smo.list) else 0
    if (type=="link") {
      fit <- matrix(eta,mo$m,ncat_1,byrow=TRUE)
      grid <- expand.grid(1:nc[2],1:nc[1])[,2:1]
      index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
      colnames(fit) <- c(paste("eta1:",1:nc[1],sep=""),
                         paste("eta2:",1:nc[2],sep=""),
                         paste("eta12:",index.lab ,sep=""))
      return(round(as.data.frame(fit),digits=digits))
    }

    if (type=="response"){
      zerovect <- rep(0,mo$m)
      eta <- cbind(zerovect,matrix(eta,mo$m,ncat_1,byrow=TRUE))
      if (object$plackett){
        p <- mMPORF2(eta,ncat1,ncat2,ncat12,ncat1_1,ncat2_1,ncat12_1,ncat,0,object$type)
      } else {
        pta <- t((object$ta+1e-06)/sum(object$ta+1e-06))
        prob <- matrix(pta,nrow=mo$m,ncol=ncol(pta),byrow=T)
        ris <- mapply(function(a1,a2) 
                        t(find.pi3(a1,a2,mo$M,mo$C,t(mo$C),object$acc2,object$maxit2)),
                        split(eta,1:mo$m), split(log(prob),1:mo$m), SIMPLIFY=FALSE)
        p <- do.call("rbind",ris)   
      }      
      fit <- p 
      grid <- expand.grid(1:ncat2,1:ncat1)[,2:1]
      index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
      colnames(fit) <- paste("p",index.lab,sep="")
      return(round(fit,digits=digits))
    }
    
    if (type=="terms"){
      maxnlev <- sum(sapply(mo$nlev,max))
      eta.x <- matrix(,mo$m,maxnlev)
      which.coef <- c()
      lox <- lapply(ox,length)
      cumint <- c(0,cumsum(lox[-3]))
      cumcoef <- cumsum(lox)
      cumnc <- c(0,nc[1],sum(nc[1:2]))
      which.coef.list <- coef.list <- list()
      b.smooth <- vector("list",length(object$beta.smooth))
      for (j in 1:length(ox)){
        which.coef <- c( which.coef,ox[[j]]%in%ox[[j]][-(1:nc[j])])
        which.coef.list[[j]] <- ox[[j]]%in%ox[[j]][-(1:nc[j])]
        coef.list[[j]] <- oo[(cumint[j]+1):cumcoef[j]]
      }
      lox2 <- lapply(coef.list,function(x)length(na.omit(x)))
      cumint2 <- c(0,cumsum(lox2[-3]))
      cumcoef2 <- cumsum(lox2)
      mar <- c("mar1","mar2","ass12")
      xnames <- c()
      k <- 0
      nterms <- termlab <- list()
      for (j in 1:length(ox)){
        nlevint <- mo$nlev[[j]]
        termlab[[j]] <- attr(terms(mo)[[j]],"term.labels")
        which.xnames <- c("(Intercept)", termlab[[j]])
        aaa <- factor(nlevint, labels = which.xnames)
        aaa2 <- factor(mo$nlev0[[j]], labels = which.xnames)
        asgn <- split(order(nlevint), aaa)
        asgn$"(Intercept)" <- NULL
        nterms <- length(asgn)
        i1 <- 0
        k2 <- 0
        if (nterms>0) {
          for (i in 1:nterms){
            partial.eta <- 0
            k <- k + 1
            noterms <- !(aaa2%in%levels(aaa)[i+1])  
            notcoe <- coef.list[[j]][-(1:nc[j])][noterms[-1]]
            notcoe[is.na(notcoe)]<-0 
            coe <- coef.list[[j]][-(1:nc[j])][!noterms[-1]]
            coe[is.na(coe)]<-0 
            browser()
            x <- mo$Xlist[[j]][,asgn[[i]],drop=TRUE] 
            notx <- if (length(notcoe)>0) notcoe else 0
            mean.x <- 0
            if (sum(noterms[-1])>0){
              mean.x <- if (sum(noterms[-1])>1) 
                mean(mo$Xlist[[j]][,noterms,drop=T][,-1]%*%notx)  else                                                                                    
                  mean(notx*mo$Xlist[[j]][,noterms,drop=T][,-1])    
            }
            
            if(length(x)>0){
              
              if (substr(termlab[[j]][i],1,3)=="pvc")
                partial.eta <- mean.x else
                  partial.eta <- if (ncol(as.matrix(x))>1) x%*%coe + mean.x else x*coe + mean.x
            } 
            mean.eta <- mean(mo$Xlist[[j]]%*%coef.list[[j]]) 
            is.smooth.x <- ((!is.null(os[[j]]))&&(length(os)>=j)&&
                              (colnames(mo$X[[j]])[i+1] %in% os[[j]]))
            x.eta <- 0
            if (is.smooth.x){
              i1 <- i1 +1
              if (substr(termlab[[j]][i],1,3)=="pvc"){
                nby  <- sum(attr(mo$Xlist[[j]],"dimnames")[[2]] %in%  termlab[[j]][i])
                for (i2 in i1:(i1+nby-1))
                  x.eta <- x.eta + drop(matrix(object$etalist[[i2]],mo$m,byrow=T)[,j])
              } else {
                lu <- object$wh.eq[[i1]]
                l1 <- c(1,1+cumsum(lu)[-length(lu)])
                l2 <- cumsum(lu)
                b.smooth.range <- l1[j]:l2[j]
                x.eta <- drop(attributes(mo$fsmooth[[j]][[i1]])$"X"%*%object$beta.smooth[[i1]][b.smooth.range])
              }
            } 
            eta.x[,k] <- oo[1 + cumint[j]] + partial.eta - mean.eta + x.eta
          }
        } 
        xnames <- c(xnames,if (length(termlab[[j]])>0) paste(mar[j],":", termlab[[j]],sep=""))
      }
      if (ncol(eta.x)>0) fit <- eta.x else {
        fit <- eta.x[,, drop=F]
        rownames(fit) <- 1:nrow(fit)
        return(fit)
      }
      colnames(fit) <- xnames
      fit <- round(as.matrix(fit),digits=digits)
      rownames(fit) <- 1:nrow(fit)
      attr(fit,"constant") <- mapply(function(x,b,n)sum(colMeans(x[,-1,drop=F])*na.omit(b[-(1:n)]))+b[1],
                                     mo$X,coef.list,nc,SIMPLIFY=T)
      if (object$ncat1>2)
        attr(fit, "mar1Contrasts") <- coef.list[[1]][2:nc[1]]-coef.list[[1]][1]
      if (object$ncat2>2)   
        attr(fit, "mar2Contrasts") <- coef.list[[2]][2:nc[2]]-coef.list[[2]][1]
      if ((object$ncat1>2)|(object$ncat2>2)){    
        attr(fit, "ass12Contrasts") <- coef.list[[3]][2:nc[3]]-coef.list[[3]][1]
        warning('predictions are made with respect the first intercept of mar1, 
                mar2 and ass12, respectively. For category-dependent predictions,
                add the respective contrasts given as attributes to this object') 
      }  
      return(fit)
    }

    if (type=="count") {
      stop("type='count' is not allowed when newdata is provided")
    }
  }
}



sqrtM <- function(W, power)  with(eigen(W), vectors %*% (values^power * t(vectors))) 

residuals.pblm <- function(object,type=c("working","pearson"), ...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  type <- match.arg(type)
  res <- switch(type,
                "working" = object$w.res,
                "pearson" = t(mapply(function(x,w) if (all(w==0)) w%*%x else sqrtM(w,.5)%*%x,                 
                                  split(object$w.res,1:nrow(object$w.res)),
                                  lapply(as.data.frame(object$W2),
                                         function(x)sym.tri(x,ncol(object$w.res))))))
  grid <- expand.grid(1:(object$ncat1-1),1:(object$ncat2-1))[,2:1]
  index.lab <- drop(t(apply(grid,1,function(x)paste(x,collapse=""))))
  colnames(res) <- c(paste("res1:",1:(object$ncat1-1),sep=""),
                     paste("res2:",1:(object$ncat2-1),sep=""),
                     paste("res12:",index.lab ,sep=""))
  as.data.frame(res)
}

resid.pblm <- residuals.pblm

vcov.pblm <- function(object,...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  V <- solve(object$IMp)
  return(V)
}

coef.pblm <- function(object,digits = max(3, getOption("digits") - 3),...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  print.default(format(drop(object$coef), digits = digits), print.gap = 2, quote = FALSE)
  invisible(object$coef)
}

coefficients.pblm <- function(object,digits = max(3, getOption("digits") - 3),...){
  if (!inherits(object,"pblm")) stop("object is not of class pblm")
  print.default(format(drop(object$coef), digits = digits), print.gap = 2, quote = FALSE)
  invisible(object$coef)
}


Rsq.pblm <- function (object,data, type = c("Cox Snell", "Cragg Uhler", "both")) 
{
  type <- match.arg(type)
  if (class(object)!="pblm") 
    stop("this is design for pblm objects only")
  fo1 <- as.formula(paste(paste(object$fo.list[[1]])[[2]],"~1",sep=""))
  ty <- object$type
  weights <- object$weights
  ncat1 <- object$ncat1
  n <- sum(object$Y*weights)
  m0 <- pblm(fo1,data=data,type=ty,weights=weights,ncat1=ncat1)
  rsq1 <- 1 - exp((2/n) * (logLik(m0)[1] - logLik(object)[1]))
  rsq2 <- rsq1/(1 - exp((2/n) * logLik(m0)[1]))
  if (type == "Cox Snell") 
    return(rsq1)
  if (type == "Cragg Uhler") 
    return(rsq2)
  if (type == "both") 
    return(list(CoxSnell = rsq1, CraggUhler = rsq2))
}

chisq.test.pblm <- function(obj){
  s <- summary(obj)
  Ei <- predict(obj,type="count")
  Xsq <- sum((obj$Y-Ei)^2/(Ei+0.0000001))
  Px <- 1-pchisq(Xsq,s$df.res)
  Gsq <- s$deviance
  Pg <- 1-pchisq(Gsq,s$df.res)
  cat("Goodness-of-fit statistics for pblm objects:\n")
  cat("X^2 =",round(Xsq,digits=3)," df =",s$df.res," p =",Px,"\n")
  cat("G^2 =",round(Gsq,digits=3)," df =",s$df.res," p =",Pg,"\n")
  if (any(obj$Y<5)) warning("Chi-square approximation may be poor")
  return(invisible(list(Xsq=list(Statistic=Xsq,p=Px,df=s$df.res),
                        Gsq=list(Statistic=Gsq,p=Pg,df=s$df.res))))
}

plot2.pblm <- function(object,model="association",v12=NULL,gr=NULL,
                       only.beta=FALSE,ynames=NULL,type="b"){
  if (is.null(ynames)) ynames <- object$ynames
  if (model=="marginal") {
    par(mfrow=c(1,2))
    plot(1:(object$xx1[1]),object$coef[1:(object$xx1[1])],type=type,ylab="log(odds)",
         xlab="intercept index",main=ynames[1])
    plot(1:(object$xx2[1]),object$coef[(sum(object$xx1)+1):(sum(object$xx1)+object$xx2[1])],
         type=type,ylab="log(odds)",xlab="intercept index",main=ynames[2])
  }
  if (model=="association") {
    if (is.null(gr)) {
      g <- expand.grid("Y1"=1:object$xx1[1],"Y2"=1:object$xx2[1])
      if (is.null(ynames)) colnames(g)=object$ynames else colnames(g)=ynames
      
      g$logOR <- object$coef[(sum(object$xx1+object$xx2)+1):(sum(object$xx1+object$xx2)+object$xx12[1])]
      wireframe(logOR ~ Y1*Y2, data = g,zlim=c(min(g$logOR-1),max(g$logOR+1)),
                scales = list(arrows = FALSE, tick.number=c(object$xx1[1],object$xx2[1])), 
                screen = list(z = -130, x = -70),xlab=ynames[1],ylab=ynames[2])
    }  else  {
      if (only.beta) {
        g<-expand.grid("Y1"=1:object$xx1[1],"Y2"=1:object$xx2[1])
         if (is.null(ynames)) colnames(g)=object$ynames else colnames(g)=ynames
        g$beta<- object$coef[rownames(object$coef) %in% paste(object$label,"lor12.",1:object$xx12[1],": ",v12,"1",sep="")]
        colnames(g)=object$ynames
        wireframe(beta ~ Y1*Y2,data = g,zlim=c(min(g$beta-1),max(g$beta+1)),
                  scales = list(arrows = FALSE,tick.number=c(object$xx1[1],object$xx2[1])),
                  screen = list(z = -130, x = -70),xlab=ynames[1],ylab=ynames[2])
      } else {
        g<-expand.grid("Y1"=1:object$xx1[1],"Y2"=1:object$xx2[1],gr)
        if (is.null(ynames)) colnames(g)=c(object$ynames,"gr") else colnames(g)=c(ynames,"gr")
        beta <- object$coef[rownames(object$coef) %in% paste(object$label,"gg.lor12.",1:object$xx12[1],": ",v12,"1",sep="")]
        logOR <- object$coef[(sum(object$xx1+object$xx2)+1):(sum(object$xx1+object$xx2)+object$xx12[1])]
        g$logOR <- c(logOR,logOR+beta)
        
        wireframe(logOR ~ Y1*Y2 , groups=gr,data = g,zlim=c(min(g$logOR)-1,max(g$logOR)+1),
                  scales = list(arrows = FALSE,tick.number=c(object$xx1[1],object$xx2[1])),
                  screen = list(z = -130, x = -70),xlab=ynames[1],ylab=ynames[2])
      }
    }
  }
}



gof<-function(object,probab){ 
  z <- object
  if (!inherits(z,"pblm")) stop("object is not of class pblm")
  ll <- sum(z$Y*log(z$p))
  MSEL  <- sum(z$weights*(rowSums((z$p-probab)^2)))/z$n
  MRSEL  <- sum(z$weights*(rowSums((z$p-probab)^2/probab)))/z$n
  MEL  <- sum(z$weights*( rowSums(probab*log(probab/z$p))))/z$n
  gAIC <- AIC(z,k=c(2,z$gaic.m,log(z$n)))
  rval <- data.frame("CONV"=as.numeric(z$convergence),"MSEL"=MSEL,"MRSEL"=MRSEL,
                     "MEL"=MEL,"AIC"=gAIC[1],"GAIC"=gAIC[2],"BIC"=gAIC[3])
  rval
}




