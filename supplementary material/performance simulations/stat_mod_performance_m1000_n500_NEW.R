############################# Section 4, Table 1 ###############################
library(pblm)
source("add_fun.R")
# 3x3 case
# two dicothomous covariate 
# TRUE MODEL: Not Uniform association and Not Proportional Odds Model (NUNPOM)
param<-c(-.6,.6,.3,-.3,  -.6,.6,-.3,.3,  2.6,2.4,2.0,1.7,-.4,.2,-.5,.5)

set.seed(234)
n <- 500  # sample size
m <- 1000  # number of samples
ncat1 <- 3
ncat2 <- 3
ncat <- ncat1 * ncat2
v1 <- runif(n,-1,1)
x <- cbind(rep(1,n),v1)
eta <- matrix(0,n,9) 
for (j in 1:n){
  X1 <- kronecker(t(x[j,]),diag(1,2))
  X2 <- X1
  X3 <- kronecker(t(x[j,]),diag(1,4))
  X <- pblm:::Mdiag(list(X1,X2,X3))
  eta[j,] <- c(0,X%*%param)
}

probab <- t(apply(eta,1,function(x)MPORF(x,ncat1,ncat2, type="gg")))


beta0 <- c("beta101","beta102","beta11",
           "beta201","beta202","beta21",
           "beta301","beta302","beta303","beta304","beta31")

beta1 <- c("beta101","beta102","beta111","beta112",
           "beta201","beta202","beta211","beta212",
           "beta301","beta302","beta303","beta304",
           "beta311","beta312","beta313","beta314")

lambda<-c(0,1,10,100)
ris0<-matrix(0,nrow=m,ncol=5+length(beta0)) 
colnames(ris0)<-c("CONV","MSEL","MRSEL","MEL","AIC",beta0)

ris1 <- vector("list",length(lambda))
for (i in 1:length(lambda)){
  ris1[[i]]<-matrix(0,nrow=m,ncol=6+length(beta1))
  colnames(ris1[[i]])<-c("CONV","MSEL","MRSEL","MEL","AIC","lambda",beta1)
}

for (j in 1:m){ 
  Y <- make.Y(probab)
  dati<-data.frame("Y1"=Y[,1],"Y2"=Y[,2],v1)     
  mod0 <- try(pblm(fo1=cbind(Y1,Y2)~v1,data=dati))
  if (class(mod0)!="try-error"){
    if (mod0$convergence) conv<-TRUE
    coe0 <- mod0$coef
    names(coe0) <- beta0
    values0 <- c(gof(mod0,probab)[c("CONV","MSEL","MRSEL","MEL","AIC")],coe0)
    ris0[j,]<-unlist(values0)
  }     
  
  for (i in 1:length(lambda)) {
    mod1 <- NULL
    conv<-FALSE
    values1 <- NA

    mod1 <- try(pblm(fo1=cbind(Y1,Y2)~v1, data=dati,
                     control = pblm.control(maxit=50),
                     proportional = pblm.prop(prop1=c(F,F),prop2=c(F,F),
                                                 prop12=c(F,F)),
                     penalty = pblm.penalty(pnl.type="ARC1",lam1=c(0,lambda[i]),
                                            lam2=c(0,lambda[i]), 
                                            lam3=c(0,lambda[i]))))                                                                                   
    if ((class(mod1)!="try-error")&&(mod1$convergence)){
      if (mod1$convergence) conv<-TRUE
      coe1 <- mod1$coef
      names(coe1) <- beta1
      ris1[[i]][j,]<-unlist(c(gof(mod1,probab)[c("CONV","MSEL","MRSEL","MEL","AIC")],lambda[i],coe1))
    }     
  }
  cat("j=",j,"\n")
}




ris1_0<-subset(ris1[[1]],ris1[[1]][,1]!=0) 
ris10 <-ris1_0[,-c(1,6)]

ris1_1<-subset(ris1[[2]],ris1[[2]][,1]!=0) 
ris11 <-ris1_1[,-c(1,6)]

ris1_10<-subset(ris1[[3]],ris1[[3]][,1]!=0) 
ris110 <-ris1_10[,-c(1,6)]

ris1_100<-subset(ris1[[4]],ris1[[4]][,1]!=0) 
ris1100 <-ris1_100[,-c(1,6)]


a0<-matrix(0,5,20)
a0[1,]<-colMeans(ris10)
a0[2,]<-colMeans(ris11)
a0[3,]<-colMeans(ris110)
a0[4,]<-colMeans(ris1100)
c0 <- colMeans(ris0[,-1])
co0 <- c(c0[1:6],rep(c0[7],2),c0[8:9],rep(c0[10],2),c0[11:14],rep(c0[15],4))
a0[5,]<-co0


newa0<-matrix(0,5,5)
newa0[,1:4] <- a0[,1:4]
newa0[,5]<-apply(a0[,5:20],1,function(x)sum(abs((x-param)/param)))
colnames(newa0)<-c("MSEL","MRSEL","MEL","AIC","relative bias")
model.name<-c("NPOM(lambda=0)","NPOM(lambda=1)",
              "NPOM(lambda=10)","NPOM(lambda=100)",
              "UPOM")
Ns<-c(nrow(ris10),nrow(ris11),nrow(ris110),nrow(ris1100))
a1<-data.frame("Model"=model.name,round(newa0,4),"Ns"=c(Ns,nrow(ris0)))
a1



library(xtable)
xtable(a1,digits=4)

be <- rbind(round(colMeans(ris10[,-c(1:4)]),digits=2),
            round(colMeans(ris11[,-c(1:4)]),digits=2),
            round(colMeans(ris110[,-c(1:4)]),digits=2),
            round(colMeans(ris1100[,-c(1:4)]),digits=2),
            round(colMeans(ris0[,-c(1:5)]),digits=2))
be <- data.frame(Model=model.name,be)
library(xtable)
xtable(be,digits=2)
sim_m1000_n500<-list(param=param,lambda=lambda,ris0=ris0,
                      ris1=ris1,a1=a1,seed=1234,be=be)
save(sim_m1000_n500,file="sim_m1000_n500_NEW.RData")


ris1 <- sim_m1000_n500$ris1
param <- sim_m1000_n500$param
a1 <- sim_m1000_n500$a1
ris0 <- sim_m1000_n500$ris0
be <- sim_m1000_n500$be
