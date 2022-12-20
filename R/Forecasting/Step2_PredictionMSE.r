
##########################################################
############ CALCULO DEL ERROR CUADRATICO MEDIO ##########
##########################################################

library(splines)
library(mgcv)
library(nlme)
library(MASS)
library(mvtnorm)
library(maptools)
library(RColorBrewer)
library(splancs)
library(spdep)

setwd("C:/APLICACIONES/Jaione/ForecastingANOVAPspline/Prediccion")


source("dumpsdataPred/riskspredict")


source("dumpsdata/dG1")
source("dumpsdata/dG2")
source("dumpsdata/dG3")

source("dumpsdata/dG4")
source("dumpsdata/dG5")
source("dumpsdata/dG6")


source("dumpsdata/V.inv")
source("dumpsdata/XVX.inv")
source("dumpsdata/P")
source("dumpsdata/V")

source("dumpsdata/Tr2t")
source("dumpsdata/Tr2st")



source("dumpsdata/F1")
source("dumpsdata/F2")
source("dumpsdata/F3")
source("dumpsdata/F4")
source("dumpsdata/F5")
source("dumpsdata/F6")

source("dumpsdata/F7")
source("dumpsdata/F8")
source("dumpsdata/F9")
source("dumpsdata/F10")
source("dumpsdata/F11")

source("dumpsdata/Z")


##### MATRIZ DE COVARIANZAS PARA PARTE ESPACIAL

diagonalG1amp<-c(1/F1,1/F2,1/F3)

dim1<-length(diagonalG1amp)

G1.amp<-diag(diagonalG1amp,ncol=dim1,nrow=dim1)



##### MATRIZ DE COVARIANZAS AMPLIADA PARA PARTE TEMPORAL

if(dim(t(as.matrix(E3)))[1]==1){
TE3<-t(Tr2t)%*%E3} else

if(dim(t(as.matrix(E3)))[1]>1){
TE3<-t(Tr2t)%*%t(E3)
}

if(dim(t(as.matrix(E3)))[1]==1){
E3T<-t(E3)%*%Tr2t} else

if(dim(t(as.matrix(E3)))[1]>1){
E3T<-E3%*%Tr2t
}


diagonalG2amp<-c(1/F4)

dim1<-length(diagonalG2amp)

G2amp<-diag(diagonalG2amp,ncol=dim1,nrow=dim1)


GTE3<-G2amp%*%TE3
E3TG<-E3T%*%G2amp

I.F3<-diag(1,nrow=dimF3[1],ncol=dimF3[2])

G2.amp<-matrix(0,nrow=dim(G2amp)[1]+dim(E3TG)[1],ncol=dim(G2amp)[2]+dim(GTE3)[2])
G2.amp[1:dim(G2amp)[1],1:dim(G2amp)[2]]<-G2amp
G2.amp[1:dim(G2amp)[1],dim(G2amp)[2]+1:dim(GTE3)[2]]<--GTE3
G2.amp[dim(G2amp)[1]+1:dim(E3TG)[1],1:dim(G2amp)[2]]<--E3TG
G2.amp[dim(G2amp)[1]+1:dim(E3TG)[1],dim(G2amp)[2]+1:dim(GTE3)[2]]<-I.F3+I.F3%*%E3T%*%G2amp%*%TE3%*%I.F3








##### MATRIZ DE COVARIANZAS AMPLIADA PARA INTERACCION

if(dim(t(as.matrix(E3)))[1]==1){
TE3I<-t(Tr2st)%*%kronecker(E3,Is)} else

if(dim(t(as.matrix(E3)))[1]>1){
TE3I<-t(Tr2st)%*%kronecker(t(E3),Is)
}

if(dim(t(as.matrix(E3)))[1]==1){
E3IT<-kronecker(t(E3),Is)%*%Tr2st} else

if(dim(t(as.matrix(E3)))[1]>1){
E3IT<-kronecker(E3,Is)%*%Tr2st
}


diagonalG3amp<-c(1/F5,1/F6,1/F7,1/F8,1/F9,1/F10,1/F11)

dim1<-length(diagonalG3amp)

G3amp<-diag(diagonalG3amp,ncol=dim1,nrow=dim1)


GTE3I<-G3amp%*%TE3I
E3ITG<-E3IT%*%G3amp

I.F3<-diag(1,nrow=dimF3[1],ncol=dimF3[2])


I.F<-kronecker(I.F3,Is)

G3.amp<-matrix(0,nrow=dim(G3amp)[1]+dim(E3ITG)[1],ncol=dim(G3amp)[2]+dim(GTE3I)[2])
G3.amp[1:dim(G3amp)[1],1:dim(G3amp)[2]]<-G3amp
G3.amp[1:dim(G3amp)[1],dim(G3amp)[2]+1:dim(GTE3I)[2]]<--GTE3I
G3.amp[dim(G3amp)[1]+1:dim(E3ITG)[1],1:dim(G3amp)[2]]<--E3ITG
G3.amp[dim(G3amp)[1]+1:dim(E3ITG)[1],dim(G3amp)[2]+1:dim(GTE3I)[2]]<-I.F+I.F%*%E3IT%*%G3amp%*%TE3I%*%I.F



#### MATRIZ DE VARIANZAS Y COVARIANZAS COMPLETA


G.amp<-matrix(0,nrow=dim(G1.amp)[1]+dim(G2.amp)[1]+dim(G3.amp)[1],ncol=dim(G1.amp)[2]+dim(G2.amp)[2]+dim(G3.amp)[2])
G.amp[1:dim(G1.amp)[1],1:dim(G1.amp)[2]]<-G1.amp
G.amp[(dim(G1.amp)[1]+1):(dim(G1.amp)[1]+dim(G2.amp)[1]),(dim(G1.amp)[2]+1):(dim(G1.amp)[2]+dim(G2.amp)[2])]<-G2.amp
G.amp[(dim(G1.amp)[1]+dim(G2.amp)[1]+1):(dim(G1.amp)[1]+dim(G2.amp)[1]+dim(G3.amp)[1]),(dim(G1.amp)[2]+dim(G2.amp)[2]+1):(dim(G1.amp)[2]+dim(G2.amp)[2]
+dim(G3.amp)[2])]<-G3.amp


Zst.amp<-cbind(Z.b31,Z.b41,Z.b51)+Z.b61%*%(-E3IT)

Z.amp<-matrix(0,dim(Z.b10)[1],ncol=dim(Z)[2])
Z.amp[1:dim(Z.b10)[1],1:dim(Z)[2]]<-cbind(Z.b10,Z.b21+Z.b22%*%(E3%*%Tr2t),Zst.amp)


Qpred<-(X0 - Z.amp%*%G%*%t(Z)%*%V.inv%*%X)



M1<-G1.amp
M2<-G2.amp[,1:dim(G2amp)[2]]
M3<-G3.amp[,1:dim(G3amp)[2]]


d1<-dim(M1)[1]+dim(M2)[1]+dim(M3)[1]
d2<-dim(M1)[2]+dim(M2)[2]+dim(M3)[2]

M<-diag(0,nrow=d1,ncol=d2)
M[1:dim(M1)[1],1:dim(M1)[2]]<-M1
M[dim(M1)[1]+1:dim(M2)[1],dim(M1)[2]+1:dim(M2)[2]]<-M2
M[dim(M1)[1]+dim(M2)[1]+1:dim(M3)[1],dim(M1)[2]+dim(M2)[2]+1:dim(M3)[2]]<-M3


Z.estr<-Z0

g1.estr<-diag(Z.estr%*%(G.amp-M%*%t(Z)%*%V.inv%*%Z%*%t(M))%*%t(Z.estr))

g1<-diag(Z.amp%*%(G-G%*%t(Z)%*%V.inv%*%Z%*%G)%*%t(Z.amp))

g2<-diag(Qpred%*%XVX.inv %*%t(Qpred))




source("dumpsdata/dV1")
source("dumpsdata/dV2")
source("dumpsdata/dV3")
source("dumpsdata/dV3")
source("dumpsdata/dV4")
source("dumpsdata/dV5")
source("dumpsdata/dV6")

source("dumpsdata/dG1")
source("dumpsdata/dG2")
source("dumpsdata/dG3")
source("dumpsdata/dG3")
source("dumpsdata/dG4")
source("dumpsdata/dG5")
source("dumpsdata/dG6")


d.v1<-Z.amp%*%(dG1%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV1%*%V.inv )
d.v2<-Z.amp%*%(dG2%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV2%*%V.inv )
d.v3<-Z.amp%*%(dG3%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV3%*%V.inv )
d.v4<-Z.amp%*%(dG4%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV4%*%V.inv )
d.v5<-Z.amp%*%(dG5%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV5%*%V.inv )
d.v6<-Z.amp%*%(dG6%*%t(Z)%*%V.inv-G%*%t(Z)%*%V.inv%*%dV6%*%V.inv )


PdV1<-P %*% dV1
PdV2<-P %*% dV2
PdV3<-P %*% dV3
PdV4<-P %*% dV4
PdV5<-P %*% dV5
PdV6<-P %*% dV6


G.s<-matrix(0, 6, 6)


G.s[1, 1]<- (1/2) * (sum(diag(PdV1 %*% PdV1)))
G.s[1, 2]<- (1/2) * (sum(diag(PdV1 %*% PdV2)))
G.s[1, 3]<- (1/2) * (sum(diag(PdV1 %*% PdV3)))
G.s[1, 4]<- (1/2) * (sum(diag(PdV1 %*% PdV4)))
G.s[1, 5]<- (1/2) * (sum(diag(PdV1 %*% PdV5)))
G.s[1, 6]<- (1/2) * (sum(diag(PdV1 %*% PdV6)))


G.s[2, 1]<- (1/2) * (sum(diag(PdV2 %*% PdV1)))
G.s[2, 2]<- (1/2) * (sum(diag(PdV2 %*% PdV2)))
G.s[2, 3]<- (1/2) * (sum(diag(PdV2 %*% PdV3)))
G.s[2, 4]<- (1/2) * (sum(diag(PdV2 %*% PdV4)))
G.s[2, 5]<- (1/2) * (sum(diag(PdV2 %*% PdV5)))
G.s[2, 6]<- (1/2) * (sum(diag(PdV2 %*% PdV6)))


G.s[3, 1]<- (1/2) * (sum(diag(PdV3 %*% PdV1)))
G.s[3, 2]<- (1/2) * (sum(diag(PdV3 %*% PdV2)))
G.s[3, 3]<- (1/2) * (sum(diag(PdV3 %*% PdV3)))
G.s[3, 4]<- (1/2) * (sum(diag(PdV3 %*% PdV4)))
G.s[3, 5]<- (1/2) * (sum(diag(PdV3 %*% PdV5)))
G.s[3, 6]<- (1/2) * (sum(diag(PdV3 %*% PdV6)))


G.s[4, 1]<- (1/2) * (sum(diag(PdV4 %*% PdV1)))
G.s[4, 2]<- (1/2) * (sum(diag(PdV4 %*% PdV2)))
G.s[4, 3]<- (1/2) * (sum(diag(PdV4 %*% PdV3)))
G.s[4, 4]<- (1/2) * (sum(diag(PdV4 %*% PdV4)))
G.s[4, 5]<- (1/2) * (sum(diag(PdV4 %*% PdV5)))
G.s[4, 6]<- (1/2) * (sum(diag(PdV4 %*% PdV6)))


G.s[5, 1]<- (1/2) * (sum(diag(PdV5 %*% PdV1)))
G.s[5, 2]<- (1/2) * (sum(diag(PdV5 %*% PdV2)))
G.s[5, 3]<- (1/2) * (sum(diag(PdV5 %*% PdV3)))
G.s[5, 4]<- (1/2) * (sum(diag(PdV5 %*% PdV4)))
G.s[5, 5]<- (1/2) * (sum(diag(PdV5 %*% PdV5)))
G.s[5, 6]<- (1/2) * (sum(diag(PdV5 %*% PdV6)))



G.s[6, 1]<- (1/2) * (sum(diag(PdV6 %*% PdV1)))
G.s[6, 2]<- (1/2) * (sum(diag(PdV6 %*% PdV2)))
G.s[6, 3]<- (1/2) * (sum(diag(PdV6 %*% PdV3)))
G.s[6, 4]<- (1/2) * (sum(diag(PdV6 %*% PdV4)))
G.s[6, 5]<- (1/2) * (sum(diag(PdV6 %*% PdV5)))
G.s[6, 6]<- (1/2) * (sum(diag(PdV6 %*% PdV6)))

G.inv<-solve(G.s)

g3<-c(rep(0,nrow(Z.amp)))
for (i in 1:nrow(Z.amp)){
d.v12<-rbind(d.v1[i,],d.v2[i,],d.v3[i,],d.v4[i,],d.v5[i,],d.v6[i,])
g3[i]<-sum(diag(d.v12%*%V%*%t(d.v12)%*%G.inv))
}

dump("g1.estr","dumpsdataPred/g1.estr")
dump("g1","dumpsdataPred/g1")
dump("g2","dumpsdataPred/g2")
dump("g3","dumpsdataPred/g3")


##CALCULAMOS EL ECM DE LOS LOG-RIESGOS

mse.amp<-g1.estr+g2+2*g3


#### ORGANIZAMOS LOS ECM EN FORMA MATRICIAL DE FORMA QUE CADA FILA
##CORRESPONDA UNA PROVINCIA DISTINTA

mse.array<-NULL
for(i in 1:p) {
lim1<-50*(i-1)+1
lim2<-50*i
mse.array<-cbind(mse.array,mse.amp[c(lim1:lim2)])
 }


#### AÑADIMOS ADEMÁS LOS ECM RIESGOS OBTENIDOS EN EL AJUSTE. LA MATRIZ
## msecom  GUARDARÁ LOS VALORES DEL ECM DE LOS RIESGOS AJUSTADOS Y PREDICHOS.
source("dumpsdata/mse")
msecom<-cbind(mse,mse.array)


#### GUARDAMOS LA MATRIZ DE LOS ECM EN LA CARPETA dumpsdataPred.
dump("msecom","dumpsdataPred/msecom")
