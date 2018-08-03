############################### Creating pdf ############################# 
pdf(file="sim1_ncat7_ncov2dic_all.pdf",width=18,height=10)

load(file="sim_n1000_ncat7_ncov2dic.RData")
ri1000 <- ris_ok
df1000 <- df
a.empirical1000 <- a.empirical
load(file="sim1_n2000_ncat7_ncov2dic.RData")
ri2000 <- ris_ok
df2000 <- df
a.empirical2000 <- a.empirical
load(file="sim1_n5000_ncat7_ncov2dic.RData")
ri5000 <- ris_ok
df5000 <- df
a.empirical5000 <- a.empirical

wlam <- which(lam %in% c(0,1,5,10,50))  
par(mfrow=c(3,5),mar=c(5,4,3,2)+0.1)
for (k in wlam){
   #### simulated values 
   a <- round(a.empirical1000[[k]],digits=3)
   main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                           list(la=lam[k], al=a, def=round(df1000[k],digits=2),
                                en=1000,em=nrow(ri1000[[k]]), 
                                expression=""))
   hist(ri1000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
        main=list(main.text,cex=1.5),
        xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
   lines(seq(0,max(ri1000[[k]]$LRT),l=500),dchisq(seq(0,max(ri1000[[k]]$LRT),l=500),df1000[k]),col=2,lty=2)   
   abline(v=qchisq(0.95,df1000[k]),col=2,lty=2)
   abline(v=quantile(ri1000[[k]]$LRT,0.95),col=1)
   #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

for (k in wlam){
  #### simulated values 
  a <- round(a.empirical2000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                          list(la=lam[k], al=a, def=round(df2000[k],digits=2),
                               en=2000,em=nrow(ri2000[[k]]), 
                               expression=""))
  hist(ri2000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri2000[[k]]$LRT),l=500),dchisq(seq(0,max(ri2000[[k]]$LRT),l=500),df2000[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df2000[k]),col=2,lty=2)
  abline(v=quantile(ri2000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}


for (k in wlam){
  #### simulated values 
  a <- round(a.empirical5000[[k]],digits=3)
  main.text <- substitute(expression(lambda==la, alpha==al, df==def, n==en,m==em),
                          list(la=lam[k], al=a, def=round(df5000[k],digits=2),
                               en=5000,em=nrow(ri5000[[k]]), 
                               expression=""))
  hist(ri5000[[k]]$LRT,freq=F,breaks=20,cex.axis=1.5,cex=1.5,
       main=list(main.text,cex=1.5),
       xlab=list(expression(LR[P]),cex=1.5),ylab=list("Density",cex=1.5))
  lines(seq(0,max(ri5000[[k]]$LRT),l=500),dchisq(seq(0,max(ri5000[[k]]$LRT),l=500),df5000[k]),col=2,lty=2)   
  abline(v=qchisq(0.95,df5000[k]),col=2,lty=2)
  abline(v=quantile(ri5000[[k]]$LRT,0.95),col=1)
  #legend("right",c("theoretical","simulated"),lty=2:1,col=2:1,cex=1.5)
}

dev.off()
