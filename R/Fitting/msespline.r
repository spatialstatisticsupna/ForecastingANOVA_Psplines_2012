
######################## calcular MSE 

g1<-diag(Z%*%(G-G%*%t(Z)%*%V.inv%*%Z%*%G)%*%t(Z))
g2<-diag((X-Z%*%G%*%t(Z)%*%V.inv%*%X)%*%ginv(XVX)%*%t(X-Z%*%G%*%t(Z)%*%V.inv%*%X))


d.v1<-Z%*%(dG1%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV1%*%V.inv )
d.v2<-Z%*%(dG2%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV2%*%V.inv )
d.v3<-Z%*%(dG3%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV3%*%V.inv )
d.v4<-Z%*%(dG4%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV4%*%V.inv )
d.v5<-Z%*%(dG5%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV5%*%V.inv )
d.v6<-Z%*%(dG6%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV6%*%V.inv )


g3<-c(rep(0,nrow(Z)))

for (i in 1:nrow(Z)){
d.v12<-rbind(d.v1[i,],d.v2[i,],d.v3[i,],d.v4[i,],d.v5[i,],d.v6[i,])
g3[i]<-sum(diag(d.v12%*%V%*%t(d.v12)%*%G.inv))
}


mse<-g1+g2+2*g3

mse1<-g1+g2

######################## Dibujar SMRs 

betas2<-a

alphas2<-b


b.ij2<-X%*%betas2 + Z%*%alphas2

radjS<-exp(b.ij2)

yadjS<-radjS*esp




Devian<-2*sum(y*log(y/yadjS)-(y-yadjS))
Devian

matrizH<-X%*%XVX.inv%*%t(X)%*%V.inv+(Z%*%G%*%t(Z)%*%P)
Df<-sum(diag(matrizH))
Df

AIC<-Devian+2*Df
AIC

BIC<-Devian+log(dim(Z)[1])*Df
BIC

AICyBIC<-cbind(Devian,Df,AIC,BIC)

IlinfS<-exp(b.ij2+qnorm(0.025)*sqrt(mse))

IlsupS<-exp(b.ij2+qnorm(0.975)*sqrt(mse))
    
smrs<-prost$casos/prost$esp


smr<-matrix(smrs,nrow=50,ncol=34)
IlinfS95<-matrix(IlinfS,nrow=50,ncol=34)
IlsupS95<-matrix(IlsupS,nrow=50,ncol=34)
risks<-matrix(radjS,nrow=50,ncol=34)
mse<-matrix(mse,nrow=50,ncol=34)
yadjS<-matrix(yadjS,nrow=50,ncol=34)



dump("smr",file="dumpsdata/smr")
dump("IlinfS95",file="dumpsdata/IlinfS95")
dump("IlsupS95",file="dumpsdata/IlsupS95")
dump("mse",file="dumpsdata/mse")
dump("risks",file="dumpsdata/risks")
dump("a",file="dumpsdata/a")
dump("b",file="dumpsdata/b")

