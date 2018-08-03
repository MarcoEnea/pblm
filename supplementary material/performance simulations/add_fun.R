restrict <- function(x,ncat1,ncat2){
  k<-0
  for (i in 1:ncat1) {
    for (j in 1:ncat2) {
      k<-k+1
      if (which(as.logical(x))==k)  {
        risposta <- c(i,j)
        break
      }
    }
  }
  risposta
}

make.Y<-function(probab,...){
  mm2 <- apply(probab,1,function(x)rmultinom(1,1,x))
  Y<- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
  Y
}


MPORF <- function (a, ncat1, ncat2, type, ...) 
{
  
  ncat12 <- ncat1 + ncat2
  ncat1_1 <- ncat1 - 1
  ncat2_1 <- ncat2 - 1
  ncat12_1 <- ncat12 - 1
  ncat <- ncat1 * ncat2
  
  if (type == "ss") 
    a[, 2:(ncat1_1 + ncat2_1 + 1)] <- -1 * a[, 2:(ncat1_1 + 
                                                    ncat2_1 + 1)]
  P <- A <- matrix(0, ncat1, ncat2)
  k <- 0
  for (i in 1:(ncat1_1)) {
    for (j in 1:(ncat2_1)) {
      k <- k + 1
      i1 <- i + 1
      j1 <- j + 1
      A[i1, j1] <- pblm:::PORF(a[c(i1, ncat1 + j, ncat12_1 + k)])
      P[i, j] <- A[i, j] - A[i, j1] - A[i1, j] + A[i1, 
                                                   j1]
    }
  }
  P[-ncat1, ncat2] <- 1/(1 + exp(-a[2:ncat1])) - if (ncat1 == 
                                                     2) 
    P[-ncat1, -ncat2]
  else (rowSums(P[-ncat1, -ncat2]) + c(0, 1/(1 + exp(-a[2:(ncat1_1)]))))
  P[ncat1, -ncat2] <- 1/(1 + exp(-a[(ncat1 + 1):(ncat12_1)])) - 
    if (ncat2 == 2) 
      P[-ncat1, -ncat2]
  else (colSums(P[-ncat1, -ncat2]) + c(0, 1/(1 + exp(-a[(ncat1 + 
                                                           1):(ncat12 - 2)]))))
  onesP <- 1 - sum(P)
  P[ncat1, ncat2] <- pmax(0, onesP)
  as.vector(t(P))
}

gof<-function(object,probab){         #diagnostics and goodness of fit
  z <- object
  if (!inherits(z,"pblm")) stop("object is not of class pblm")
  ll <- logLik(z,penalized=T)
  MSEL  <- sum(z$weights*(rowSums((z$p-probab)^2)))/z$n
  MRSEL  <- sum(z$weights*(rowSums((z$p-probab)^2/probab)))/z$n
  MEL  <- sum(z$weights*( rowSums(probab*log(probab/z$p))))/z$n
  rval <- data.frame("CONV"=as.numeric(z$convergence),"MSEL"=MSEL,"MRSEL"=MRSEL,
                     "MEL"=MEL,"AIC"=AIC(z),"BIC"=AIC(z,k=log(z$n)))
  rval
}

GOR<-function(x,l=0){
  x <- x + l
  G <- matrix(0,nrow(x)-1,ncol(x)-1)
  for (j in 1:(ncol(x)-1)){
    for (i in 1:(nrow(x)-1)){
      G[i,j]<- log((sum(x[1:i,1:j])*sum(x[(i+1):nrow(x),(j+1):ncol(x)]))/
                       (sum(x[1:i,(j+1):ncol(x)])*sum(x[(i+1):nrow(x),1:j])))
    }    
  }      
  G
}  

