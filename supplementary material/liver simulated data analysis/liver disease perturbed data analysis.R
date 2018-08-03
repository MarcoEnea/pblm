##########################################################################
############ Simulated Liver Disease Patients Data Analysis ##############
##########################################################################
## Author: Marco Enea
##
## Warning:
## This example is based on data simulated from the original Liver Disease 
## Patients Data. Accordingly, all results in this example might not
## correspond to those of the original analysis reported on the paper. 
## Further, by changing the set.seed(), model selection might be different, 
## while the code refers to the original analysis. 
##########################################################################
library(pblm)
da <- read.csv2("pertliver.csv",header=T)


############################ generating simulated data ###################
load("mod7_predicted_probability.RData")
source("add_fun.R")
set.seed(13) 
mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
Y <- t(apply(mm2,2,function(x) restrict(x,ncat1=3,ncat2=3)))
da$STAGE <- Y[,1]
da$STIFF <- Y[,2]
da <- subset(da,!(STAGE==1&da$STIFF==3)) # this record is not present into the original data
ftable(STIFF~STAGE,da)
##########################################################################

# type "?pblm" for details on the package.

# Trying to estimate a NUPPOM model on set P^0_1 (specified in the paper) 
# does not work as the original liver data don't contain a patient profile with 
# STAGE=1 and STIFF=3. 

mod1 <- try(pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, data=da, 
             proportional=pblm.prop(prop2=c(F,T,T,T,F)),center=TRUE))

class(mod1)=="try-error" 

# In this case, applying a penalization (e.g. ARC1) on the association 
# intercepts can be useful.

# Step1: grid search of the minimun lambda value for which the estimating 
# algorithm does not fail. 

# In this example we use grid search typical values. Using a finer scale, 
# would take more time    

lam.min <- c(0.1,0.2,0.5,1,2,5,10,20,50)

# This takes a while.
lmin <- NA
for (i in 1:length(lam.min)){
  mod1 <- try(pblm(fo1=cbind(STAGE,STIFF) ~ SEX + ALT + AGE + PLT, 
                   data=da, 
                   proportional=pblm.prop(prop2=c(F,T,T,T,F)),
                   penalty=pblm.penalty(lam3=lam.min[i],pnl.type="ARC1"),
                   center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  if ((class(mod1)!="try-error"&&mod1$convergence==TRUE)){
     lmin <- lam.min[i]
    cat("minimum lambda =",lam.min[i],"\n") 
    break
  }
}



# Step 2. Searching the optimal lambda on a sequence of values starting 
# from lam_min  
aic1 <- AIC(mod1)
lam <- c(10,20,50,100,200) # it would be better a finer scale...
pos.min <- 1
for (i in 1:length(lam)){
  mod_i <- try(pblm(fo1=cbind(STAGE,STIFF) ~ SEX + ALT + AGE + PLT, 
                   data=da, 
                   proportional=pblm.prop(prop2=c(F,T,T,T,F)),
                   penalty=pblm.penalty(lam3=c(lam[i]),pnl.type="ARC1"),
                   center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  if ((class(mod1)!="try-error"&&mod1$convergence==TRUE)) {
    if (AIC(mod_i) < aic1) {
      aic1 <- AIC(mod1) 
      mod1 <- mod_i
      pos.min <- i
    }
  }
}
cat(AIC(mod1)," and lambda =",lam[pos.min],"\n") 


lam.opt <- lam[pos.min]
sm1 <- summary(mod1)

# Model 2 is a NUPOM fitted on set P^0_2 
mod2 <- pblm(fo1=cbind(STAGE,STIFF) ~ SEX + ALT + AGE + PLT, 
             data=da, 
             proportional=pblm.prop(prop2=c(F,T,T,T,T)),#<--- the default
             penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
             center=TRUE,control=pblm.control(maxit=50,acc=1e-06)) 
AIC(mod2)
sm2 <- summary(mod2)


# LR_p test for the comparison mod1 vs mod2: using the theoretical 
# chi-square distribution. 
df_12 <- 1
LRP_12 <- 2*(logLik(mod1,penalized=T)-logLik(mod2,penalized=T))
print(pvalue_12 <- 1-pchisq(LRP_12,df_12))  
# in this case is preferable mod2
########################################################################

########################################################################
# Simulating the LR_p test, for the comparison mod1 vs mod2, to check if 
# the theoretical chi-square distribution with 1 degrees of freedom well 
# approximate the simulated one.
#
######### WARNIG: this simulation takes about one hour ################# 
#
# H0: NUPOM (mod2)
# H1: NUPPOM (mod1)
ncat1 <- 3
ncat2 <- 3
n <- nrow(da)
m <- 1500 
ris_12 <- matrix(0, nrow=m,ncol=4)
source("add_fun.R")
probab <- mod2$p
set.seed(123)
for (i in 1:m) {
  cat("i=",i,"\n")
  mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
  Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
  da2 <- da
  da2$STAGE <- Y[,1]
  da2$STIFF <- Y[,2]

  # the model under H1
  mod1b <- try(pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, 
                    data=da2, 
                    proportional=pblm.prop(prop2=c(F,T,T,T,F)),
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                    center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  
    
  # the model under H0             
  mod2b <- try(pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, 
                    data=da2, 
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                    center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
    
    if ((class(mod1b) != "try-error")&(class(mod2b) != "try-error")){                            
      sum1b <- summary(mod1b)
      sum2b <- summary(mod2b)
      ris_12[i,1] <- sum1b$convergence * sum2b$convergence
      ris_12[i,2] <- 2*(logLik(mod1b,penalized=T)-logLik(mod2b,penalized=T))  
      ris_12[i,3] <- max(sum1b$tol)
      ris_12[i,4] <- max(sum2b$tol)
    }
  }
  
  colnames(ris_12) <- c("convergence","LRT","tol1","tol2")
  ris_12 <- as.data.frame(ris_12)
  ris_12ok <- subset(ris_12, LRT > 0 & (convergence==1)) 
 
  ########################### Figure SM-1 ###################################
  pdf("LR_p_mod1_vs_mod2_SM-1.pdf",width = 12,height=8)
  par(mfrow=c(1,1))
  
  # Histogram
  a <- round(sum(ris_12ok$LRT>qchisq(0.95,df_12))/nrow(ris_12ok), digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                          list(la=0.5, al=a, em=nrow(ris_12ok), 
                               expression=""))
  hist(ris_12ok$LRT,freq=F,breaks=20,ylim=c(0,1.5),cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),xlim=c(0,max(ris_12ok$LRT)),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,12,l=500),dchisq(seq(0,12,l=500),df_12),col=2,lty=2)   
  abline(v=qchisq(0.95,df_12),col=2,lty=2)
  abline(v=quantile(ris_12ok$LRT,0.95),col=1)
  dev.off()
  # QQ-plot
  nulldist <- rchisq(5000,1)
  qqplot(nulldist,ris_12ok$LRT,xlab=list("theoretical",cex=1.5),
         ylab=list("simulated",cex=1.5),cex.axis=1.5,cex=1.5,
         xlim=c(0,15),ylim=c(0,15),main=list(main.text,cex=1.5))
  abline(0,1)
  dev.off()

############################ end of Figure SM-1 #############################
    

############################################################################

#############################################################################

mod3 <- pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, data=da, 
             proportional=pblm.prop(prop2=c(F,T,T,T,F), prop12=c(T,T,T,T,T)),
             center=TRUE)
AIC(mod3)
sm3 <- summary(mod3)

### mod1 vs mod3
print(df_13 <- round(sm1$df.tot-sm3$df.tot,digits=0)) # approximately 1
print(LRP_13 <- 2*(logLik(mod1,penalized=T)-logLik(mod3,penalized=F)))
print(1-pchisq(LRP_13,df_13))  # mod1 is kept 
###############################################################################

###############################################################################

mod4 <- pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, data=da, 
             proportional=pblm.prop(prop12=c(T,T,T,T,T)),
             center=TRUE)
AIC(mod4)
sm4 <- summary(mod4)

### mod1 vs mod4
print(df_14 <- round(sm1$df.tot-sm4$df.tot,digits=0)) # approximately 1
print(LRP_14 <- 2*(logLik(mod1,penalized=T)-logLik(mod4,penalized=F)))
print(1-pchisq(LRP_14,df_14))  # mod1 is kept 
################################################################################

################################################################################

mod5 <- pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, data=da, 
             penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
             proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T,T)),
             center=TRUE,control=pblm.control(maxit=50,acc=1e-06))
AIC(mod5)
sm5 <- summary(mod5)

### mod1 vs mod5
df_15 <- 3
print(LRP_15 <- 2*(logLik(mod1,penalized=T)-logLik(mod5,penalized=T)))
print(1-pchisq(LRP_15,df_15))  # mod5 is kept 

###################### LR_P simulation for the comparison mod5 vs mod1 ####################
#H0: NUPOM (mod5)
#H1: NUPPOM (mod1)
ncat1 <- 3
ncat2 <- 3
n <- nrow(da)
m <- 1500
ris_15 <- matrix(0, nrow=m,ncol=4)
source("add_fun.R")
probab <- mod5$p
set.seed(123)
for (i in 1:m) {
  cat("i=",i,"\n")
  mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
  Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
  da2 <- da
  da2$STAGE <- Y[,1]
  da2$STIFF <- Y[,2]
  
  # the model under H1
  mod1b <- try(pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, data=da2, 
                    proportional=pblm.prop(prop2=c(F,T,T,T,F)),
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),center=TRUE,
                    control=pblm.control(maxit=50,acc=1e-06)))
  
  
  # the model under H0             
  mod5b <- try(pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, data=da2, 
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                    proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T,T)),
                    center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  
  if ((class(mod1b) != "try-error")&(class(mod5b) != "try-error")){                            
    sum1b <- summary(mod1b)
    sum5b <- summary(mod5b)
    ris_15[i,1] <- sum1b$convergence * sum5b$convergence
    ris_15[i,2] <- 2*(logLik(mod1b,penalized=T)-logLik(mod5b,penalized=T))  
    ris_15[i,3] <- max(sum1b$tol)
    ris_15[i,4] <- max(sum5b$tol)
  }
}

colnames(ris_15) <- c("convergence","LRT","tol1","tol5")
nulldist <- rchisq(500,1)
ris_15 <- as.data.frame(ris_15)
ris_15ok <- subset(ris_15, LRT > 0 & convergence==1)


pdf("LR_p_mod1_vs_mod5_liver.pdf",width = 12,height=8)
par(mfrow=c(1,2))

a <- round(sum(ris_15ok$LRT>qchisq(0.95,df_15))/nrow(ris_15ok),digits=3)
main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                        list(la=0.5, al=a, em=nrow(ris_15ok), 
                             expression=""))
hist(ris_15ok$LRT,freq=F,breaks=20,ylim=c(0,0.25),cex.axis=1.5,cex=1.5,
     main=list(main.text,cex=1.5),xlim=c(0,max(ris_15ok$LRT)+5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_15),col=2,lty=2)   
abline(v=qchisq(0.95,df_15),col=2,lty=2)
abline(v=quantile(ris_15ok$LRT,0.95),col=1)
legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)

# QQ-plot
nulldist <- rchisq(5000,df_15)
qqplot(nulldist,ris_15ok$LRT,xlab=list("theoretical",cex=1.5),
       ylab=list("simulated",cex=1.5),cex.axis=1.5,
       xlim=c(0,15),ylim=c(0,15),main=list(main.text,cex=1.5))
abline(0,1)
dev.off()

###############################################################################

################################################################################

mod6 <- pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, 
             fo12=~ALT + AGE, 
             data=da, 
             penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
             proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T)),
             center=TRUE)
AIC(mod6)
sm6 <- summary(mod6)

### mod5 vs mod6
df_56 <- 1 # approximately 1
print(LRP_56 <- 2*(logLik(mod5,penalized=T)-logLik(mod6,penalized=T)))
print(1-pchisq(LRP_56,df_56))  # mod6 is kept 

###################### LR_P simulation for the comparison mod6 vs mod5 ####################

#H0: NUPOM (mod6)
#H1: NUPPOM (mod5)
ncat1 <- 3
ncat2 <- 3
n <- nrow(da)
m <- 1500
ris_56 <- matrix(0, nrow=m,ncol=4)
source("add_fun.R")
probab <- mod6$p
set.seed(123)
for (i in 1:m) {
  cat("i=",i,"\n")
  #generating simulated responses under H0
  mm2 <- apply(probab,1,function(x) rmultinom(1,1,x)) 
  Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
  da2 <- da
  da2$STAGE <- Y[,1]
  da2$STIFF <- Y[,2]
  
  # the model under H1
  mod5b <- try(pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, data=da2, 
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                    proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T,T)),
                    center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  
  # the model under H0             
  mod6b <- try(pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, 
               fo12=~ALT + AGE, 
               data=da2, 
               penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
               proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T)),
               center=TRUE,control=pblm.control(maxit=50,acc=1e-06)))
  
  if ((class(mod5b) != "try-error")&(class(mod6b) != "try-error")){                            
    sum5b <- summary(mod5b)
    sum6b <- summary(mod6b)
    ris_56[i,1] <- sum5b$convergence * sum6b$convergence
    ris_56[i,2] <- 2*(logLik(mod5b,penalized=T)-logLik(mod6b,penalized=T))  
    ris_56[i,3] <- max(sum5b$tol)
    ris_56[i,4] <- max(sum6b$tol)
  }
}

colnames(ris_56) <- c("convergence","LRT","tol5","tol6")
nulldist <- rchisq(500,1)
ris_56 <- as.data.frame(ris_56)
ris_56ok <- subset(ris_56, LRT > 0 & convergence==1)


pdf("LR_p_mod5_vs_mod6_liver.pdf",width = 12,height=8)
par(mfrow=c(1,2))

a <- round(sum(ris_56ok$LRT>qchisq(0.95,df_56))/nrow(ris_56ok),digits=3)
main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                        list(la=0.5, al=a, em=nrow(ris_56ok), 
                             expression=""))
hist(ris_56ok$LRT,freq=F,breaks=20,ylim=c(0,0.25),cex.axis=1.5,cex=1.5,
     main=list(main.text,cex=1.5),xlim=c(0,max(ris_56ok$LRT)+5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_56),col=2,lty=2)   
abline(v=qchisq(0.95,df_56),col=2,lty=2)
abline(v=quantile(ris_56ok$LRT,0.95),col=1)
legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)

# QQ-plot
nulldist <- rchisq(5000,df_56)
qqplot(nulldist,ris_56ok$LRT,xlab=list("theoretical",cex=1.5),
       ylab=list("simulated",cex=1.5),cex.axis=1.5,
       xlim=c(0,15),ylim=c(0,15),main=list(main.text,cex=1.5))
abline(0,1)
dev.off()


################################################################################

################################################################################

mod7 <- pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, 
             fo12=~ALT, 
             data=da, 
             penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
             proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T)),
             center=TRUE)
AIC(mod7)
sm7 <- summary(mod7)

### 1. mod6 vs mod7 and ...
df_67 <- 1 
print(LRP_67 <- 2*(logLik(mod6,penalized=T)-logLik(mod7,penalized=T)))
print(1-pchisq(LRP_67,df_67))  # mod7 is kept 

### 2. ... mod1 vs mod7
df_17 <- 5
print(LRP_17 <- 2*(logLik(mod1,penalized=T)-logLik(mod7,penalized=T)))
print(1-pchisq(LRP_17,df_17))  # mod7 is kept 

###################### LR_P simulation for the comparison mod7 vs mod1 ####################

#H0: NUPPOM (mod7)
#H1: NUPPOM (mod6)
ncat1 <- 3
ncat2 <- 3
n <- nrow(da)
m <- 1500
ris_67 <- matrix(0, nrow=m,ncol=4)
ris_17 <- matrix(0, nrow=m,ncol=4)

source("add_fun.R")
probab <- mod7$p
set.seed(123)
for (i in 1:m) {
  cat("i=",i,"\n")
  mm2 <- apply(probab,1,function(x) rmultinom(1,1,x))
  Y <- t(apply(mm2,2,function(x) restrict(x,ncat1,ncat2)))
  da2 <- da
  da2$STAGE <- Y[,1]
  da2$STIFF <- Y[,2]
  
  # the model under H1 (for testing mod6 vs mod7 )
  mod6b <- try(pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, 
                    fo12=~ALT + AGE, 
                    data=da2, 
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                    proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T,T)),
                    center=TRUE))
  
  # the model under H1 (for testing mod1 vs mod7 )
  mod1b <- try(pblm(fo1=cbind(STAGE,STIFF)~ SEX + ALT + AGE + PLT, data=da2, 
                    proportional=pblm.prop(prop2=c(F,T,T,T,F)),
                    penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),center=TRUE))
  
  # the model under H0             
  mod7b <- try(pblm(fo1=cbind(STAGE,STIFF)~ ALT + AGE + PLT, 
                fo12=~ALT, 
                data=da2, 
                penalty=pblm.penalty(lam3=lam.opt,pnl.type="ARC1"),
                proportional=pblm.prop(prop2=c(F,T,T,F), prop12=c(F,T)),
                center=TRUE))
  
  if ((class(mod6b) != "try-error")&(class(mod7b) != "try-error")){                            
    sum6b <- summary(mod6b)
    sum7b <- summary(mod7b)
    ris_67[i,1] <- sum6b$convergence * sum7b$convergence
    ris_67[i,2] <- 2*(logLik(mod6b,penalized=T)-logLik(mod7b,penalized=T))  
    ris_67[i,3] <- max(sum6b$tol)
    ris_67[i,4] <- max(sum7b$tol)
  }
  if ((class(mod1b) != "try-error")&(class(mod7b) != "try-error")){                            
    sum1b <- summary(mod1b)
    sum7b <- summary(mod7b)
    ris_17[i,1] <- sum1b$convergence * sum7b$convergence
    ris_17[i,2] <- 2*(logLik(mod1b,penalized=T)-logLik(mod7b,penalized=T))  
    ris_17[i,3] <- max(sum1b$tol)
    ris_17[i,4] <- max(sum7b$tol)
  }
}

colnames(ris_67) <-  c("convergence","LRT","tol6","tol7")
colnames(ris_17) <-  c("convergence","LRT","tol1","tol7")
ris_67 <- as.data.frame(ris_67)
ris_67ok <- subset(ris_67, LRT > 0 & convergence==1)
ris_17 <- as.data.frame(ris_17)
ris_17ok <- subset(ris_17, LRT > 0 & convergence==1)


pdf("LR_p_mod6_vs_mod7_liver.pdf",width = 12,height=8)
par(mfrow=c(1,2))

a <- round(sum(ris_67ok$LRT>qchisq(0.95,df_67))/nrow(ris_67ok),digits=3)
main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                        list(la=0.5, al=a, em=nrow(ris_67ok), 
                             expression=""))
hist(ris_67ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
     main=list(main.text,cex=1.5),xlim=c(0,max(ris_67ok$LRT)+5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_67),col=2,lty=2)   
abline(v=qchisq(0.95,df_67),col=2,lty=2)
abline(v=quantile(ris_67ok$LRT,0.95),col=1)
legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)

# QQ-plot
nulldist <- rchisq(5000,df_67)
qqplot(nulldist,ris_67ok$LRT,xlab=list("theoretical",cex=1.5),
       ylab=list("simulated",cex=1.5),cex.axis=1.5,
       xlim=c(0,15),ylim=c(0,15),main=list(main.text,cex=1.5))
abline(0,1)
dev.off()


###### mod1 vs mod7
pdf("LR_p_mod1_vs_mod7_liver.pdf",width = 12,height=8)
par(mfrow=c(1,2))

a <- round(sum(ris_17ok$LRT>qchisq(0.95,df_17))/nrow(ris_17ok),digits=3)
main.text <- substitute(expression(lambda==la, alpha==al, m==em),
                        list(la=0.5, al=a, em=nrow(ris_17ok), 
                             expression=""))
hist(ris_17ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
     main=list(main.text,cex=1.5),xlim=c(0,max(ris_17ok$LRT)+5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_17),col=2,lty=2)   
abline(v=qchisq(0.95,df_17),col=2,lty=2)
abline(v=quantile(ris_17ok$LRT,0.95),col=1)
legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)

# QQ-plot
nulldist <- rchisq(5000,df_17)
qqplot(nulldist,ris_17ok$LRT,xlab=list("theoretical",cex=1.5),
       ylab=list("simulated",cex=1.5),cex.axis=1.5,
       xlim=c(0,15),ylim=c(0,15),main=list(main.text,cex=1.5))
abline(0,1)
dev.off()

#############################  Figure 8 ##################################
pdf("LR_p_liver.pdf",width = 12,height=8)
par(mfrow=c(2,3),mar= c(5, 4, 4, 2) +0.1)

## mod1 vs mod2
a <- round(sum(ris_12ok$LRT>qchisq(0.95,df_12))/nrow(ris_12ok), digits=3)
main.text <- substitute(expression("mod1 vs mod2",df==d, alpha==al, "m*"==em),
                        list(d=1, al=a, em=nrow(ris_12ok), 
                             expression=""))
hist(ris_12ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,lwd=2,
     main=list(main.text,cex=1.5),#xlim=c(0,max(ris_12ok$LRT)),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,12,l=500),dchisq(seq(0,12,l=500),df_12),col=2,lty=2,lwd=2)   
abline(v=qchisq(0.95,df_12),col=2,lty=2,lwd=2)
abline(v=quantile(ris_12ok$LRT,0.95),col=1,lwd=2)

## mod1 vs mod5
a <- round(sum(ris_15ok$LRT>qchisq(0.95,df_15))/nrow(ris_15ok),digits=3)
main.text <- substitute(expression("mod1 vs mod5",df==d, alpha==al, "m*"==em),
                        list(d=3, al=a, em=nrow(ris_15ok), 
                             expression=""))
hist(ris_15ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,lwd=2,
     main=list(main.text,cex=1.5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_15),col=2,lty=2,lwd=2)   
abline(v=qchisq(0.95,df_15),col=2,lty=2,lwd=2)
abline(v=quantile(ris_15ok$LRT,0.95),col=1,lwd=2)

## mod5 vs mod6
a <- round(sum(ris_56ok$LRT>qchisq(0.95,df_56))/nrow(ris_56ok),digits=3)
main.text <- substitute(expression("mod5 vs mod6", df==d, alpha==al, "m*"==em),
                        list(d=1, al=a, em=nrow(ris_56ok), 
                             expression=""))
hist(ris_56ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,lwd=2,
     main=list(main.text,cex=1.5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_56),col=2,lty=2,lwd=2)   
abline(v=qchisq(0.95,df_56),col=2,lty=2,lwd=2)
abline(v=quantile(ris_56ok$LRT,0.95),col=1,lwd=2)

## mod6 vs mod7
a <- round(sum(ris_67ok$LRT>qchisq(0.95,df_67))/nrow(ris_67ok),digits=3)
main.text <- substitute(expression("mod6 vs mod7",df==d, alpha==al, "m*"==em),
                        list(d=1, al=a, em=nrow(ris_67ok), 
                             expression=""))
hist(ris_67ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,lwd=2,
     main=list(main.text,cex=1.5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_67),col=2,lty=2,lwd=2)   
abline(v=qchisq(0.95,df_67),col=2,lty=2,lwd=2)
abline(v=quantile(ris_67ok$LRT,0.95),col=1,lwd=2)

## mod1 vs mod7
a <- round(sum(ris_17ok$LRT>qchisq(0.95,df_17))/nrow(ris_17ok),digits=3)
main.text <- substitute(expression("mod1 vs mod7",df==d, alpha==al, "m*"==em),
                        list(d=5, al=a, em=nrow(ris_17ok), 
                             expression=""))
hist(ris_17ok$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,lwd=2,
     main=list(main.text,cex=1.5),
     xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df_17),col=2,lty=2,lwd=2)   
abline(v=qchisq(0.95,df_17),col=2,lty=2,lwd=2)
abline(v=quantile(ris_17ok$LRT,0.95),col=1,lwd=2)

plot.new()
legend("center",c("theoretical","simulated"),lty=2:1,col=2:1,cex=2,lwd=2)
dev.off()



#################################################################################








 

