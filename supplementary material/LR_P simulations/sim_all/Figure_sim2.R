############################### Creating pdf ############################# 
pdf(file="sim2_all.pdf",width=18,height=10)

load(file="sim2_n200_ncat3_ncov2dic.RData")
ri200 <- ris_ok
df200 <- 3
a.empirical200 <- a.empirical
load(file="sim2_n500_ncat5_ncov2dic.RData")
ri500 <- ris_ok
df500 <- 15
a.empirical500 <- a.empirical
load(file="sim2_n1000_ncat7_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- 35
for (k in 1:9) {
  ri1000[[k]] <-ri1000[[k]][ri1000[[k]]$LRT<100,] # removing outliers
  a.empirical[[k]] <- sum(ri1000[[k]]$LRT>qchisq(0.95,df1000))/nrow(ri1000[[k]])
  }
a.empirical1000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10))  
par(mfrow=c(3,4),mar=c(5,4,3,2)+0.1)
for (k in wlam){
   #### simulated values 
   a <- round(a.empirical200[[k]],digits=3)
   main.text <- substitute(expression(lambda==la, alpha==al, df==def, ncat==3,n==en,m==em),
                           list(la=lam[k], al=a, def=round(df200,digits=2),
                                en=200,em=nrow(ri200[[k]]), 
                                expression=""))
   hist(ri200[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
        main=list(main.text,cex=1.5),ylim=c(0,0.24),
        xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
   lines(seq(0,max(ri200[[k]]$LRT),l=500),dchisq(seq(0,max(ri200[[k]]$LRT),l=500),df200),col=2,lty=2)   
   abline(v=qchisq(0.95,df200),col=2,lty=2)
   abline(v=quantile(ri200[[k]]$LRT,0.95),col=1)
   #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical500[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, ncat==5,n==en,m==em),
                          list(la=lam[k], al=a, def=round(df500,digits=2),
                               en=500,em=nrow(ri500[[k]]), 
                               expression=""))
  hist(ri500[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.08),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri500[[k]]$LRT),l=500),dchisq(seq(0,max(ri500[[k]]$LRT),l=500),df500),col=2,lty=2)   
  abline(v=qchisq(0.95,df500),col=2,lty=2)
  abline(v=quantile(ri500[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical1000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def,ncat==7, 
                                     n==en,m==em),
                          list(la=lam[k], al=a, def=round(df1000,digits=2),
                               en=1000,em=nrow(ri1000[[k]]), expression=""))
  hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.06),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri1000[[k]]$LRT),l=500),dchisq(seq(0,max(ri1000[[k]]$LRT),l=500),df1000),col=2,lty=2)   
  abline(v=qchisq(0.95,df1000),col=2,lty=2)
  abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

dev.off()



############################### Creating pdf ############################# 
pdf(file="sim2_all_nsup.pdf",width=18,height=10)

load(file="sim2_n500_ncat3_ncov2dic.RData")
ri500 <- ris_ok
df500 <- 3
a.empirical500 <- a.empirical

load(file="sim2_n1000_ncat5_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- 15
a.empirical1000 <- a.empirical

load(file="sim2_n2000_ncat7_ncov2dic.RData")
ri2000 <- ris_ok
df2000 <- 35
for (k in 1:9) {
  ri2000[[k]] <-ri2000[[k]][ri2000[[k]]$LRT<100,] # removing outliers
  a.empirical[[k]] <- sum(ri2000[[k]]$LRT>qchisq(0.95,df2000))/nrow(ri2000[[k]])
}
a.empirical2000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10))  
par(mfrow=c(3,4),mar=c(5,4,3,2)+0.1)
for (k in wlam){
  #### simulated values 
  a <- round(a.empirical500[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, ncat==3,n==en,m==em),
                          list(la=lam[k], al=a, def=round(df500,digits=2),
                               en=500,em=nrow(ri500[[k]]), 
                               expression=""))
  hist(ri500[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.24),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri500[[k]]$LRT),l=500),dchisq(seq(0,max(ri500[[k]]$LRT),l=500),df500),col=2,lty=2)   
  abline(v=qchisq(0.95,df500),col=2,lty=2)
  abline(v=quantile(ri500[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical1000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, ncat==5,n==en,m==em),
                          list(la=lam[k], al=a, def=round(df1000,digits=2),
                               en=1000,em=nrow(ri1000[[k]]), 
                               expression=""))
  hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.08),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri1000[[k]]$LRT),l=1000),dchisq(seq(0,max(ri1000[[k]]$LRT),l=1000),df1000),col=2,lty=2)   
  abline(v=qchisq(0.95,df1000),col=2,lty=2)
  abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical2000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def,ncat==7, 
                                     n==en,m==em),
                          list(la=lam[k], al=a, def=round(df2000,digits=2),
                               en=2000,em=nrow(ri2000[[k]]), expression=""))
  hist(ri2000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.06),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri2000[[k]]$LRT),l=500),dchisq(seq(0,max(ri2000[[k]]$LRT),l=500),df2000),col=2,lty=2)   
  abline(v=qchisq(0.95,df2000),col=2,lty=2)
  abline(v=quantile(ri2000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

dev.off()







