library(pblm)
source("add_fun.R")
#     n: vector of sample sizes
#     m: number of replicates in all the simulations
# ncat1: vector containing the number of levels for the first response 
# ncat2: vector containing the number of levels for the second response 
#   lam: vector of penalty parameters
# param: vector of true regression parameters which the simulations are generated from. 

ncat1 <- 3
ncat2 <- 3
n = 1000 
m = 1000  
lam=c(0,0.2,0.5,1,2,5,10,20,50) 
nlam <- length(lam)
df <- c()
ris <- ris_ok <-autov <- nulldist <- mod.false <- mod.true <- vector("list",nlam)
Iee <- Ibb <- Ibe <- Ieb <- Ieeb <-  a.empirical <- vector("list",nlam)

####################################### Original parameters ###############################
#-------#  b_10   b_11    b_12    b_21    b_22   b_23    b_30         b_31         b_32  
param<-c(-.5,.5, .2,.2, -.3,.3,  -.1,.6, .3,.3, -.2,.4,  1.5,2,2.5,3, -.2,-.2,-.2,-.2, .5,.5,.5,.5)
###########################################################################################


############################# Generating Data sets ###################################
set.seed(123)

# building of design matrix
fat1 <- sample(c(0,1),n,replace=TRUE)
fat2 <- sample(c(0,1),n,replace=TRUE)
Z <- model.matrix(~ fat1 + fat2)

# generating the predictors eta=Xbeta
# this is for the 3x3 case only
eta <- matrix(0,n,9)
for (j in 1:n){
  X1 <- kronecker(t(Z[j,]),diag(1,2))
  X2 <- X1
  X3 <- kronecker(t(Z[j,]),diag(1,4))
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
                        proportional=pblm.prop(prop1=c(FALSE,TRUE,FALSE),
                                               prop2=c(FALSE,TRUE,FALSE),
                                               prop12=c(FALSE,TRUE,FALSE)),
                        penalty=pblm.penalty(pnl.type="ARC1", lam3=c(lam[k],0,0)),
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
                     proportional=pblm.prop(prop1=c(FALSE,TRUE,FALSE),
                                            prop2=c(FALSE,TRUE,FALSE),
                                            prop12=c(FALSE,TRUE,FALSE)),
                     penalty=pblm.penalty(pnl.type="ARC1",
                                          lam3=c(lam[k],0,0)),
                     control=pblm.control(maxit=50,acc=1e-07)))
    
      # the model under H0             
    mod0 <- try(pblm(fo1=fo, data=dati2, plackett=TRUE,   
                     proportional=pblm.prop(prop1=c(FALSE,TRUE,FALSE),
                                           prop2=c(FALSE,TRUE,FALSE),
                                           prop12=c(FALSE,TRUE,TRUE)),
                     penalty=pblm.penalty(pnl.type="ARC1",
                                          lam3=c(lam[k],0,0)),
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


   nulldist[[k]] <- rchisq(5000,3)
  colnames(ris[[k]]) <- c("convergence","LRT","ll1","ll0")
  ris[[k]] <- as.data.frame(ris[[k]])
  ris_ok[[k]] <- ris[[k]][(ris[[k]]$LRT>0&ris[[k]]$convergence==1),]
  a.empirical[[k]] <- sum(ris_ok[[k]]$LRT>qchisq(0.95,3))/nrow(ris_ok[[k]])
  
  
}
  ######################### Histograms ###############################
load(file="sim2_n1000_ncat3_ncov2dic.RData") 
 
pdf(file="sim2_n1000_ncat3_ncov2dic_chisq3df.pdf",width=12,height=10)
par(mfrow=c(3,3),mar=c(4,4,3,1)+0.1)
 

for (k in 1:nlam){
  #### simulated values 
  #### simulated values 
  a <- round(a.empirical[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                          list(la=lam[k], al=a, em=nrow(ris_ok[[k]]), 
                               expression=""))
  hist(ris_ok[[k]]$LRT,freq=F,breaks=20,ylim=c(0,0.25),cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),xlim=c(0,max(ris_ok[[k]]$LRT)+5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),3),col=2,lty=2)   
  abline(v=qchisq(0.95,3),col=2,lty=2)
  abline(v=quantile(ris_ok[[k]]$LRT,0.95),col=1)
  legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}
dev.off()


save.image(file="sim2_n1000_ncat3_ncov2dic.RData")
