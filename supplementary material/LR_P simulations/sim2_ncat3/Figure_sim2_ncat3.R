############################### Creating pdf ############################# 
pdf(file="sim2_ncat3_ncov2dic_all.pdf",width=18,height=10)

load(file="sim2_n200_ncat3_ncov2dic.RData")
ri200 <- ris_ok
df200 <- 3
a.empirical200 <- a.empirical
load(file="sim2_n500_ncat3_ncov2dic.RData")
ri500 <- ris_ok
df500 <- 3
a.empirical500 <- a.empirical
load(file="sim2_n1000_ncat3_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- 3
a.empirical1000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10,50))  
par(mfrow=c(3,5),mar=c(5,4,3,2)+0.1)
for (k in wlam){
   #### simulated values 
   a <- round(a.empirical200[[k]],digits=3)
   main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                           list(la=lam[k], al=a, def=round(df200,digits=2),
                                en=200,em=nrow(ri200[[k]]), 
                                expression=""))
   hist(ri200[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
        main=list(main.text,cex=1.5),ylim=c(0,0.25),
        xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
   lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df200),col=2,lty=2)   
   abline(v=qchisq(0.95,df200),col=2,lty=2)
   abline(v=quantile(ri200[[k]]$LRT,0.95),col=1)
   #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical500[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                          list(la=lam[k], al=a, def=round(df500,digits=2),
                               en=500,em=nrow(ri500[[k]]), 
                               expression=""))
  hist(ri500[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.25),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df500),col=2,lty=2)   
  abline(v=qchisq(0.95,df500),col=2,lty=2)
  abline(v=quantile(ri500[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical1000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                          list(la=lam[k], al=a, def=round(df1000,digits=2),
                               en=1000,em=nrow(ri1000[[k]]), 
                               expression=""))
  hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),ylim=c(0,0.25),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,22,l=500),dchisq(seq(0,22,l=500),df1000),col=2,lty=2)   
  abline(v=qchisq(0.95,df1000),col=2,lty=2)
  abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

dev.off()



