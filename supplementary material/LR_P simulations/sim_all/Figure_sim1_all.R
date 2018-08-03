############################### Creating pdf ############################# 
pdf(file="sim1_all.pdf",width=18,height=10)

load(file="sim_n200_ncat3_ncov2dic.RData")
ri200 <- ris_ok
df200 <- df
a.empirical200 <- a.empirical
load(file="sim_n500_ncat5_ncov2dic.RData")
ri500 <- ris_ok
df500 <- df
a.empirical500 <- a.empirical
load(file="sim_n1000_ncat5_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- df
a.empirical1000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10))  
par(mfrow=c(3,4),mar=c(5,4,3,2)+0.1)
for (k in wlam){
   #### simulated values 
   a <- round(a.empirical200[[k]],digits=3)
   main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==3,m==em),
                           list(la=lam[k], al=a, def=round(df200[k],digits=2),
                                en=200,em=nrow(ri200[[k]]), 
                                expression=""))
   hist(ri200[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
        main=list(main.text,cex=1.5),
        xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
   lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df200[k]),col=2,lty=2)   
   abline(v=qchisq(0.95,df200[k]),col=2,lty=2)
   abline(v=quantile(ri200[[k]]$LRT,0.95),col=1)
   #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical500[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==5,m==em),
                          list(la=lam[k], al=a, def=round(df500[k],digits=2),
                               en=500,em=nrow(ri500[[k]]), 
                               expression=""))
  hist(ri500[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df500[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df500[k]),col=2,lty=2)
  abline(v=quantile(ri500[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical1000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==7,m==em),
                          list(la=lam[k], al=a, def=round(df1000[k],digits=2),
                               en=1000,em=nrow(ri1000[[k]]), 
                               expression=""))
  hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df1000[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df1000[k]),col=2,lty=2)
  abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}
dev.off()




#############################################################
############################### Creating pdf ############################# 
pdf(file="sim1_all_nsup.pdf",width=18,height=10)

load(file="sim_n500_ncat3_ncov2dic.RData")
ri500 <- ris_ok
df500 <- df
a.empirical500 <- a.empirical
load(file="sim_n1000_ncat5_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- df
a.empirical1000 <- a.empirical
load(file="sim_n2000_ncat5_ncov2dic.RData")
ri2000 <- ris_ok
df2000 <- df
a.empirical2000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10))  
par(mfrow=c(3,4),mar=c(5,4,3,2)+0.1)
for (k in wlam){
  #### simulated values 
  a <- round(a.empirical500[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==3,m==em),
                          list(la=lam[k], al=a, def=round(df500[k],digits=2),
                               en=500,em=nrow(ri500[[k]]), 
                               expression=""))
  hist(ri500[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df500[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df500[k]),col=2,lty=2)
  abline(v=quantile(ri500[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical1000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==5,m==em),
                          list(la=lam[k], al=a, def=round(df1000[k],digits=2),
                               en=1000,em=nrow(ri1000[[k]]), 
                               expression=""))
  hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=1000),dchisq(seq(0,22,l=1000),df1000[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df1000[k]),col=2,lty=2)
  abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical2000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,ncat==7,m==em),
                          list(la=lam[k], al=a, def=round(df2000[k],digits=2),
                               en=2000,em=nrow(ri2000[[k]]), 
                               expression=""))
  hist(ri2000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df2000[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df2000[k]),col=2,lty=2)
  abline(v=quantile(ri2000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}
dev.off()



