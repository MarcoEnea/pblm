############# data from the analysis by Lapp et al. (1998)  ##############
# Cross-tabulation of a sample of British males according to seven 
# occupational status categories of fathers and their sons
library(pblm)
data(bms)
source("add_fun.R")
dat <- bms

############################# Table 2 ####################################
print(CT <- xtabs(freq ~ fathers + sons, data=dat))
##########################################################################

######################### The saturated model ############################
# the logLik of the saturated model ca be calculated on the 
# observed data directly
a<-(dat$freq)/sum(dat$freq)
fullog <- sum(dat$freq*log(ifelse(a==0,1,a)))
df_full <- 48 #i.e. ncat1*ncat2 - 1
aic_full <- -2*(fullog-df_full)

# a model based alternative with a very small adjustement factor for zeros
m_full <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq, 
               control=pblm.control(zero.adj = 1e-10))
AIC(m_full)


##########################################################################

###################### Figure 6, top-left plot  ##########################
# Plotting the observed association structure
require(lattice)
g <- expand.grid(1:6,1:6)
colnames(g)=c("sons","fathers")
g$logGOR <- as.vector(t(GOR(CT,.5)))
wireframe(logGOR~sons*fathers, data=g, zlim=c(0,max(g$logGOR+1)), 
          scales=list(arrows=FALSE), screen = list(z=-130, x=-60))
##########################################################################


############################### Models ###################################

# independence
m0.0 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             proportional = pblm.prop(prop12=TRUE),
             penalty=pblm.penalty(pnl.type="ridge",lam3=c(1e9)))
AIC(m0.0)
G2_00 <- 2*(fullog-logLik(m0.0))
df00 <- round((df_full-m0.0$df.fix),digit=0)
1-pchisq(G2_00,df00)

# uniform association
m0.1 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             proportional = pblm.prop(prop12=TRUE))
AIC(m0.1)
G2_01 <- 2*(fullog-logLik(m0.1))
df01 <- round((df_full-m0.1$df.fix),digit=0)
1-pchisq(G2_01,df01)

# third degree polynomial surface with approximately integer scores
m1.0 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
              penalty=pblm.penalty(pnl.type="ARC2",lam3=c(1e7), 
                                   lam4=c(1e7), s3=c(4), s4=c(4)))

# second degree polynomial surface with approximately integer scores
m1.1 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(1e7), lam4=c(1e7), 
                                  s3=c(3), s4=c(3)))
AIC(m1.1)
G2_1 <- 2*(fullog-logLik(m1.1))
df1 <- round((df_full-m1.1$df.fix),digit=0)
1-pchisq(G2_1,df1)

# AIC-based lambda selection for a surface tending to a second 
# degree polynomial surface with non integer scores. 
# The option auto.select=TRUE allows to estimate a model with the optimal 
# lambda. The underlying algorithm makes use of the R function optim().
m1.2 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(20), 
                                  lam4=c(20), s3=c(3), s4=c(3)), 
             control=pblm.control(auto.select=TRUE))

# AIC-based lambda selection for a surface tending to a third 
# degree polynomial surface with non integer scores.
m1.3 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(20), 
                                  lam4=c(20), s3=c(4), s4=c(4)), 
             control=pblm.control(auto.select=TRUE))


# First degree polynomial surface with approximately integer scores.
m1.4 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(1e7), 
                                  lam4=c(1e7), s3=c(2), s4=c(2)))
AIC(m1.4)
G2_4 <- 2*(fullog-logLik(m1.4))
df4 <- round((df_full-m1.4$df.fix),digit=0)
1-pchisq(G2_4,df4)

# AIC-based lambda selection for a surface tending to a third 
# degree polynomial surface with non integer scores.
m1.5 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(20), 
                                  lam4=c(20), s3=c(2), s4=c(2)), 
             control=pblm.control(auto.select=TRUE))

# BIC-based lambda selection
m1.6 <- pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(20), 
                                  lam4=c(20), s3=c(4), s4=c(4)), 
             control=pblm.control(auto.select=TRUE,gaic.m=log(sum(CT))))

AIC(m1.0, m1.1, m1.2, m1.3, m1.4, m1.5, m1.6)

g$logOR <- m1.2$coef[13:48]

wireframe(logOR ~ sons*fathers, data = g,zlim=c(min(g$logOR-1),max(g$logOR+1)), 
          scales = list(arrows = FALSE), screen = list(z = -130, x = -70))
#####


# Searching the optimal lambda into a sequence of values
# WARNING: running this example take a while
lam <- c(seq(0.1,0.9,by=0.1),seq(1,99,by=1),seq(100,900,by=100),
         seq(1000,10000,by=1000))
aic1 <- aic2 <-aic3 <- aic4 <-c()
for (i in 1:length(lam)){
  m.aic1 <- try(pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
                     penalty=pblm.penalty(pnl.type="ARC2",
                     lam3=lam[i], lam4=lam[i], s3=c(1), s4=c(1)),
                     control=pblm.control(maxit=50,acc=1e-06)))
  aic1[i] <- ifelse(class(m.aic1)!="try-error", m.aic1$GAIC, NA)

  m.aic2 <- try(pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
                     penalty=pblm.penalty(pnl.type="ARC2",
                     lam3=lam[i], lam4=lam[i], s3=c(2), s4=c(2))))
  aic2[i] <- ifelse(class(m.aic2)!="try-error", m.aic2$GAIC, NA)
  
  m.aic3 <- try(pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq, 
                     penalty=pblm.penalty(pnl.type="ARC2",
                     lam3=lam[i], lam4=lam[i], s3=c(3), s4=c(3))))
  aic3[i] <- ifelse(class(m.aic3)!="try-error", m.aic3$GAIC, NA)
  
  m.aic4 <- try(pblm(fo1=cbind(fathers,sons)~1, data=dat, weights=freq,
                     penalty=pblm.penalty(pnl.type="ARC2",
                     lam3=lam[i], lam4=lam[i], s3=c(4), s4=c(4))))
  aic4[i] <- ifelse(class(m.aic4)!="try-error", m.aic4$GAIC, NA)
}

min(aic1)
which.min(aic1)
lam[which.min(aic1)]

min(aic2)
which.min(aic2)
lam[which.min(aic2)]

min(aic3)
which.min(aic3)
lam[which.min(aic3)]

min(aic4)
which.min(aic4)
lam[which.min(aic4)]
############################# Figure 5 ###################################
par(mfrow=c(1,2))
plot(log(lam),aic4,type="l",ylim=c(22220,22400),ylab="AIC",
     xlab=expression(log(lambda)),lwd=2)
lines(log(lam),aic2,col=2,lty=2,lwd=2)
lines(log(lam),aic1,col=3,lty=3,lwd=2)
lines(log(lam),aic3,col=4,lty=4,lwd=2)
abline(h=m0.1$GAIC,lwd=2,lty=5,col="grey")

plot(log(lam),aic4,type="l",ylim=c(22235,22260),ylab="AIC",
     xlab=expression(log(lambda)),lwd=2)
lines(log(lam),aic2,col=2,lty=2,lwd=2)
lines(log(lam),aic1,col=3,lty=3,lwd=2)
lines(log(lam),aic3,col=4,lty=4,lwd=2)

#########################################################################

#Dale and Goodman models. See also Table 6 in Lapp et al.(1998)

# Goodman's model 
#estimation od m2.0 can fail
m2.0 <- pblm(fo1=cbind(fathers,sons)~1, type="ll", weights=freq, RC.fo=~Row*Col, data=dat)
m2.1 <- pblm(fo1=cbind(fathers,sons)~1, type="ll", weights=freq, RC.fo=~Row+Col, data=dat)
m2.2 <- pblm(fo1=cbind(fathers,sons)~1, type="ll", weights=freq, RC.fo=~Row, data=dat)
m2.3 <- pblm(fo1=cbind(fathers,sons)~1, type="ll", weights=freq, RC.fo=~Col, data=dat)
m2.4 <- pblm(fo1=cbind(fathers,sons)~1, type="ll", weights=freq, RC.fo=~1, data=dat)

# Goodman's model deviances 
c(2*(fullog-logLik(m2.0)),round((df_full-m2.0$df.fix),digit=0)) 
c(2*(fullog-logLik(m2.1)),round((df_full-m2.1$df.fix),digit=0)) 
c(2*(fullog-logLik(m2.2)),round((df_full-m2.2$df.fix),digit=0)) 
c(2*(fullog-logLik(m2.3)),round((df_full-m2.3$df.fix),digit=0)) 
c(2*(fullog-logLik(m2.4)),round((df_full-m2.4$df.fix),digit=0))


# Dale's model
#estimation od m3.0 can fail
m3.0 <- pblm(fo1=cbind(fathers,sons)~1, RC.fo=~Row*Col, weights=freq, data=dat)
m3.1 <- pblm(fo1=cbind(fathers,sons)~1, RC.fo=~Row+Col, weights=freq, data=dat)
m3.2 <- pblm(fo1=cbind(fathers,sons)~1, RC.fo=~Row, weights=freq, data=dat)
m3.3 <- pblm(fo1=cbind(fathers,sons)~1, RC.fo=~Col, weights=freq, data=dat)
m3.4 <- pblm(fo1=cbind(fathers,sons)~1, RC.fo=~1, weights=freq, data=dat)

# Dale's model Deviances (compare with Lapp et al.(1998))
c(2*(fullog-logLik(m3.1)),round((df_full-m3.1$df.fix),digit=0)) 
c(2*(fullog-logLik(m3.2)),round((df_full-m3.2$df.fix),digit=0)) 
c(2*(fullog-logLik(m3.3)),round((df_full-m3.3$df.fix),digit=0)) 
c(2*(fullog-logLik(m3.4)),round((df_full-m3.4$df.fix),digit=0))


AIC(m1.0,m1.1,m1.2,m1.3,m1.4,m2.0,m2.1,m2.2,m2.3,m2.4,m3.0,m3.1,m3.2,m3.3,m3.4)


######################### end of the analysis #############################


########## (Fast) Alternative way to use pblm with grouped data  ###########
da <- multicolumn(freq~fathers+sons,dat) 
attributes(da)
fo <- as.formula(paste(attributes(da)$"resp","~1",sep=""))

m0 <- pblm(fo, data=da, proportional=pblm.prop(prop12=T))


m1 <- pblm(fo, data=da, penalty=pblm.penalty(pnl.type="ARC2",
                                             lam3=c(1e7), lam4=c(1e7), s3=c(3), s4=c(3)))

m1.3 <- pblm(fo1=fo, data=da, 
             control=pblm.control(auto.select=TRUE), 
             penalty=pblm.penalty(pnl.type="ARC2",
                                  lam3=c(20), lam4=c(20), 
                                  s3=c(3), s4=c(3)))


# Searching the optimal lambda into a sequence of values
lam <- c(seq(0.1,0.9,by=0.1),seq(1,99,by=1),seq(100,900,by=100),
         seq(1000,10000,by=1000))
aic1 <- aic2 <-aic3 <- aic4 <-c()
for (i in 1:length(lam)){
  m.aic1 <- try(pblm(fo1=fo, data=da, penalty=pblm.penalty(pnl.type="ARC2",
                                                           lam3=lam[i], lam4=lam[i], s3=c(1), s4=c(1)),
                     control=pblm.control(maxit=50,acc=1e-06)))
  aic1[i] <- ifelse(class(m.aic1)!="try-error", m.aic1$GAIC, NA)
  
  m.aic2 <- try(pblm(fo1=fo, data=da, penalty=pblm.penalty(pnl.type="ARC2",
                                                           lam3=lam[i], lam4=lam[i], s3=c(2), s4=c(2))))
  aic2[i] <- ifelse(class(m.aic2)!="try-error", m.aic2$GAIC, NA)
  
  m.aic3 <- try(pblm(fo1=fo, data=da, penalty=pblm.penalty(pnl.type="ARC2",
                                                           lam3=lam[i], lam4=lam[i], s3=c(3), s4=c(3))))
  aic3[i] <- ifelse(class(m.aic3)!="try-error", m.aic3$GAIC, NA)
  
  m.aic4 <- try(pblm(fo1=fo, data=da, penalty=pblm.penalty(pnl.type="ARC2",
                                                           lam3=lam[i], lam4=lam[i], s3=c(4), s4=c(4))))
  aic4[i] <- ifelse(class(m.aic4)!="try-error", m.aic4$GAIC, NA)
}

min(aic1)
which.min(aic1)
lam[which.min(aic1)]

min(aic2)
which.min(aic2)
lam[which.min(aic2)]

min(aic3)
which.min(aic3)
lam[which.min(aic3)]

min(aic4)
which.min(aic4)
lam[which.min(aic4)]
############################# Figure 5 (Alternative way) ###################################
par(mfrow=c(1,2))
plot(log(lam),aic1,type="l",ylim=c(320,500),ylab="AIC",
     xlab=expression(log(lambda)),lwd=2)
lines(log(lam),aic2,col=2,lty=2,lwd=2)
lines(log(lam),aic3,col=3,lty=3,lwd=2)
lines(log(lam),aic4,col=4,lty=4,lwd=2)
abline(h=m0$GAIC,lwd=2,lty=5,col="grey")

plot(log(lam),aic1,type="l",ylim=c(320,345),ylab="AIC",
     xlab=expression(log(lambda)),lwd=2)
lines(log(lam),aic2,col=2,lty=2,lwd=2)
lines(log(lam),aic3,col=3,lty=3,lwd=2)
lines(log(lam),aic4,col=4,lty=4,lwd=2)
#########################################################################

#########################################################################


