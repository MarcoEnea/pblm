library(pblm)
source("add_fun.R")
#     n: vector of sample sizes
#     m: number of replicates in all the simulations
# ncat1: vector containing the number of levels for the first response 
# ncat2: vector containing the number of levels for the second response 
#   lam: vector of penalty parameters
# param: vector of true regression parameters which the simulations are generated from. 

ncat1 <- 7
ncat2 <- 7
n = 5000
m = 1000  
lam=c(0,0.2,0.5,1,2,5,10,20,50) 
#lam=c(20)
nlam <- length(lam)
df <- c()
ris <- ris_ok <-autov <- nulldist <- mod.false <- vector("list",nlam)
Iee <- Ibb <- Ibe <- Ieb <- Ieeb <-  a.empirical <- vector("list",nlam)

data(bms)
ta <- xtabs(freq ~ fathers + sons, data=bms)
prob <- as.vector(t(ta/sum(ta)+0.0001))
prob <- prob/sum(prob)
C <- pblm:::make.contrasts(7,7)
M <- pblm:::make.marginals(7,7,type="gg")
veroeta <- crossprod(C,log(M%*%prob))

#          b_10            b_11            b_12            
param <- c(veroeta[2:7], rep(-0.1,6),   rep(0.1,6),
#          b_20            b_21            b_23            
           veroeta[8:13], rep(-0.1,6),  rep(0.1,6),
#          b_30            b_31            b_32
           veroeta[14:49], rep(-0.1,36), rep(0,36))
########################################################################
###########################################################################################


############################# Generating Data sets ###################################
set.seed(1234)

# building of design matrix
fat1 <- sample(c(0,1),n,replace=TRUE)
fat2 <- sample(c(0,1),n,replace=TRUE)
Z <- model.matrix(~ fat1 + fat2)

# generating the predictors eta=Xbeta
# this is for the 7x7 case only
eta <- matrix(0,n,49)
for (j in 1:n){
  X1 <- kronecker(t(Z[j,]),diag(1,6))
  X2 <- X1
  X3 <- kronecker(t(Z[j,]),diag(1,36))
  X <- pblm:::Mdiag(list(X1,X2,X3))
  eta[j,] <- c(0,X%*%param) #0 is required for construction
}


# Generating the responses from the linear multiple linear predictors
 probab <- t(apply(eta,1,function(x) MPORF(x,ncat1,ncat2,type="gg")))
 mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
 Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
 dat <- data.frame("Y1"=Y[,1],"Y2"=Y[,2],fat1,fat2)
 dat <- as.data.frame(ftable(fat1+fat2~Y1+Y2,data=dat))
 dat$fat1 <- as.factor(dat$fat1)
 dat$fat2 <- as.factor(dat$fat2)
 dat2 <- multicolumn(Freq~Y2+Y1+fat1+fat2,data=dat)
########################################################################################
fo <- as.formula(paste(attributes(dat2)$"resp"," ~ fat1 + fat2",sep=""))
for (k in 1:nlam){  
  set.seed(k)
  cat("lam=",lam[k],"\n")
  ris[[k]] <- matrix(0, nrow=m,ncol=4)
  mod.false[[k]]<- try(pblm(fo1=fo, data=dat2,
                       proportional=pblm.prop(prop1=c(FALSE,TRUE,TRUE),
                                              prop2=c(FALSE,TRUE,TRUE),
                                              prop12=c(FALSE,TRUE,FALSE)),
                       penalty=pblm.penalty(pnl.type="ARC1", lam3=c(0,0,lam[k])),
                       plackett=TRUE, control=pblm.control(maxit=50,acc=1e-07)))

  for (i in 1:m) {
    cat("i=",i,"\n")
    mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
    Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
    dat <- data.frame("Y1"=Y[,1],"Y2"=Y[,2],fat1,fat2)
    dat <- as.data.frame(ftable(fat1+fat2~Y1+Y2,data=dat))
    dat$fat1 <- as.factor(dat$fat1)
    dat$fat2 <- as.factor(dat$fat2)
    dati2 <- multicolumn(Freq~Y2+Y1+fat1+fat2,data=dat)
     
    # the model under H1
    mod1 <- try(pblm(fo1=fo,
                     data=dati2, plackett=TRUE,
                     proportional=pblm.prop(prop1=c(FALSE,TRUE,TRUE),
                                            prop2=c(FALSE,TRUE,TRUE),
                                            prop12=c(FALSE,TRUE,FALSE)),
                     penalty=pblm.penalty(pnl.type="ARC1",
                                          lam3=c(0,0,lam[k])),
                     control=pblm.control(maxit=50,acc=1e-07)))
    
      # the model under H0             
    mod0 <- try(pblm(fo1=fo,
                     fo2=~fat1+fat2,
                     fo12=~fat1, data=dati2, plackett=TRUE,   
                     proportional=pblm.prop(prop1=c(FALSE,TRUE,TRUE),
                                           prop2=c(FALSE,TRUE,TRUE),
                                           prop12=c(FALSE,TRUE)),
                    control=pblm.control(maxit=50,acc=1e-07)))
                   
    if ((class(mod0) != "try-error")&(class(mod1) != "try-error")){                            
      sum1 <- summary(mod1)
      sum0 <- summary(mod0)
      ris[[k]][i,1] <- sum1$convergence * sum0$convergence
      ris[[k]][i,2] <- 2*(sum1$logLikp - sum0$logLikp)  
      ris[[k]][i,3] <- max(sum1$tol)
      ris[[k]][i,4] <- max(sum0$tol)
    }
  }

  Iee[[k]] <- mod.false[[k]]$IM[54:89,54:89]
  Ibb[[k]] <- mod.false[[k]]$IM[1:53,1:53]
  Ibe[[k]] <- mod.false[[k]]$IM[1:53,54:89]
  Ieb[[k]] <- t(Ibe[[k]])
  Ieeb[[k]] <- (Iee[[k]] - Ieb[[k]]%*%solve(Ibb[[k]])%*%Ibe[[k]])
  # mod.false$P is the same for all simulations of fixed lambda
  autov[[k]] <- eigen(Ieeb[[k]]%*%solve(Ieeb[[k]]+mod.false[[k]]$P[54:89,54:89]))$values  
  nulldist[[k]] <- rowSums(sapply(autov[[k]],function(x)x*rchisq(5000,1)))
  df[k] <- sum(autov[[k]])
  colnames(ris[[k]]) <- c("convergence","LRT","ll1","ll0")
  ris[[k]] <- as.data.frame(ris[[k]])
  ris_ok[[k]] <- ris[[k]][(ris[[k]]$LRT>0&ris[[k]]$convergence==1),]
  a.empirical[[k]] <- sum(ris_ok[[k]]$LRT>qchisq(0.95,df[k]))/nrow(ris_ok[[k]])
  
  
}
  ######################### Histograms ###############################
 
 pdf(file="sim1_n5000_ncat7_ncov2dic.pdf",width=12,height=10)
 par(mfrow=c(3,3),mar=c(4,4,3,1)+0.1)
 for (k in 1:nlam){
   #### simulated values 
   a <- round(a.empirical[[k]],digits=3)
   main.text <- substitute(expression(lambda==la, alpha==al,df==edf,n==en, m==em),
                           list(la=lam[k], al=a,edf=df, en=5000,em=nrow(ris_ok[[k]]), 
                                expression=""))
   hist(ris_ok[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
        main=list(main.text,cex=1.5),ylim=c(0,0.1),
        xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
   lines(seq(0,80,l=500),dchisq(seq(0,80,l=500),df[k]),col=2,lty=2)   
   abline(v=qchisq(0.95,df[k]),col=2,lty=2)
   abline(v=quantile(ris_ok[[k]]$LRT,0.95),col=1)
   #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
 }
 dev.off()
 
save.image(file="sim1_n5000_ncat7_ncov2dic.RData")
