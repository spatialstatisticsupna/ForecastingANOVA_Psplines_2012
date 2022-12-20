################################################################################
########### CALCULO DE LAS MATRICES PARA EL INTERVALO (1973,2008) ##############

rm(list=ls(all=TRUE))

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

################################################################################
###### CONSTRUCCIÓN DE LAS MATRICES PARA EL MODELO DE PREDICCIÓN ###############

#######Ejecutar las definiciones para construir un producto tensorial por filas
## y las bases de B-splines.

Rten<-function(X1,X2){
one1<-matrix(1,1,ncol(X1))
one2<-matrix(1,1,ncol(X2))
kronecker(X1,one2)*kronecker(one1,X2)
}

bspline<-function(x, xl, xr, ndx, bdeg){
dx <- (xr-xl)/ndx
 knots <- seq(xl-bdeg*dx, xr+bdeg*dx, by=dx)
 B <- spline.des(knots,x,bdeg+1,0*x)$design
 B
 }

### LEEMOS LOS DATOS

dat<- read.table("Data/ProstH.txt",header=TRUE)
ProstH<-dat


bdeg<-3  ## SELECCIONAMOS BASES DE ORDEN 3
pord<-2  ## SELECCIONAMOS PENALIZACIONES DE ORDEN 2

#### DEFINICIÓN DE LAS VARIABLES EXPLICATIVAS

xSC<- readShapePoly("ZonificacionEspanaProv/esp_prov.shp",IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))

x1<-coordinates(xSC)[,1]
x1<-(x1-min(x1))/1000

x2<-coordinates(xSC)[,2]
x2<-(x2-min(x2))/1000


p<-9 ### SELECCIONAMOS EL NÚMERO DE AÑOS EN LOS QUE QUEREMOS PREDECIR

x3<-unique(seq(min(ProstH$agno),max(ProstH$agno)+p)) ### AMPLIAMOS EL RANGO DEL TIEMPO PARA PREDECIR
x3<-(x3-min(x3))/1000

#### SELECCION DE NODOS

##NOTA: ndx indica el número de intervalos internos en el rango de las variables
## por tanto el número de nodos es ndx+1. Los valores de ndx1, ndx2 y ndx3 deben ser
## los mismos que los utilizados en el ajuste. Mirar archivo Inicio.r de la carpeta ajuste


ndx1<-10 #9#min(ceiling(length(x1)/4),40)#8
dis1 <- (max(x1)-min(x1))/ndx1
x1l<-min(x1)-dis1*0.05
x1r<-max(x1)+dis1*0.05
dx1 <- (x1r-x1l)/ndx1
knots1<-seq(x1l-bdeg*dx1, x1r+bdeg*dx1, by=dx1)


ndx2<-10 #9#min(ceiling(length(x2)/4),40) #x
dis2 <- (max(x2)-min(x2))/ndx2
x2l<-min(x2)-dis2*0.05
x2r<-max(x2)+dis2*0.05
dx2 <- (x2r-x2l)/ndx2
knots2<-seq(x2l-bdeg*dx2, x2r+bdeg*dx2, by=dx2)



x3na<-unique(ProstH$agno)
x3na<-(x3na-min(x3na))/1000
ndx3<-8 #min(ceiling(length(x3)/4),40)#2#
dis3 <- (max(x3na)-min(x3na))/ndx3
x3l<-min(x3na)-dis3*0.05
x3r<-max(x3na)+dis3*0.05 ### frontera hasta la no ampliada
dx3 <- (x3r-x3l)/ndx3

x3ra<-max(x3)+dis3*0.05 ### frontera para la ampliada
knots3<-seq(x3l-bdeg*dx3, x3ra+bdeg*dx3, by=dx3)



#### CONSTRUCCIÓN DE LAS BASES DE BSPLINES

B1<-spline.des(knots1,x1,bdeg+1,0*x1)$design ##Base marginal longitud
B2<-spline.des(knots2,x2,bdeg+1,0*x2)$design ##Base marginal latitud
B3<-spline.des(knots3,x3,bdeg+1,0*x3,outer.ok = TRUE)$design   ##Base marginal tiempo

B<-kronecker(B3,Rten(B2,B1))


#### CONSTRUCCIÓN DE LAS MATRICES DE PENALIZACION

##Penalización marginal longitud

m1=ncol(B1)
D1=diff(diag(m1),differences=pord)
P1.svd=svd(t(D1)%*%D1)
U1s=(P1.svd$u)[,1:(m1-pord)] # le quito dos columnas por los dos autovalores que son 0
U1n=(P1.svd$u)[,(m1-pord+1):m1] # dos columnas de los dos autovalores que son 0
d1=(P1.svd$d)[1:(m1-pord)]
Delta1=diag(d1)

##Penalización marginal latitud

m2=ncol(B2)
D2=diff(diag(m2),differences=pord)
P2.svd=svd(t(D2)%*%D2)
U2s=(P2.svd$u)[,1:(m2-pord)] # le quito dos columnas por los dos autovalores que son 0
U2n=(P2.svd$u)[,(m2-pord+1):m2] # dos columnas de los dos autovalores que son 0
d2=(P2.svd$d)[1:(m2-pord)]
Delta2=diag(d2)


Bs<-Rten(B2,B1)

##Penalización marginal tiempo

m3=ncol(B3)
D3=diff(diag(m3),differences=pord)
P3.svd=svd(t(D3)%*%D3)
U3s=(P3.svd$u)[,1:(m3-pord)] # le quito dos columnas por los dos autovalores que son 0
U3n=(P3.svd$u)[,(m3-pord+1):m3] # dos columnas de los dos autovalores que son 0
d3=(P3.svd$d)[1:(m3-pord)]
Delta3=diag(d3)

source("dumpsdata/rowsD3Ajuste")
source("dumpsdata/colsD3Ajuste")

D3antigua<-D3[1:rowsD3Ajuste,1:colsD3Ajuste]
E3<-D3[(rowsD3Ajuste+1):dim(D3)[1],1:colsD3Ajuste]
MatF3<-D3[(rowsD3Ajuste+1):dim(D3)[1],(colsD3Ajuste+1):dim(D3)[2]]

dimF3<-dim(MatF3)

##CONSTRUCCIÓN DE LAS MATRICES DE EFECTOS FIJOS Y EFECTOS ALEATORIOS
#CORRESPONDIENTES A LA REPARAMETRIZACIÓN COMO MODELO MIXTOS DEL
#MODELO DE P-SPLINES AMPLIADO



############ MATRIZ X
#source("C:/APLICACIONES/Jaione/ForecastingSANOVA/Ajuste/dumpsdata/X")


X1<-cbind(1,x1)
X2<-cbind(1,x2)

Xs<-Rten(X2,X1)
xs<-Xs[,4]
One.s<-Xs[,1]


nao<-length(unique(ProstH$agno))

Xt<-cbind(1,x3[1:nao])
Xt0<-cbind(1,x3[(nao+1):length(x3)])

One.t<-Xt[,1]
One.t0<-Xt0[,1]

xt<-x3[1:nao]
xt0<-x3[(nao+1):length(x3)]

xhat<-cbind(x1,x2,xs)

X<-cbind(kronecker(One.t,Xs),kronecker(xt,One.s),kronecker(xt,xhat)) ### ANTIGUA
X0<-cbind(kronecker(One.t0,Xs),kronecker(xt0,One.s),kronecker(xt0,xhat)) ### AMPLIADA


################ MATRIZ Z AMPLIADA

#### Z AMPLIADA PARTE ESPACIAL

Z1<-B1%*%U1s
Z2<-B2%*%U2s

Z.s<-cbind(Rten(Z2,X1),Rten(X2,Z1),Rten(Z2,Z1))
Z.b10<-kronecker(One.t0,Z.s)

#### Z AMPLIADA PARTE TEMPORAL

source("dumpsdata/Delta3")
source("dumpsdata/rowsB3Ajuste")
source("dumpsdata/colsB3Ajuste")
source("dumpsdata/U3s")

numrow<-dim(B3)[1]

Z.t<-B3[1:rowsB3Ajuste,1:colsB3Ajuste]%*%U3s
Z.b2<-kronecker(Z.t,One.s)


Z.t1<-B3[(rowsB3Ajuste+1):numrow,1:colsB3Ajuste]%*%U3s
Z.b21<-kronecker(Z.t1,One.s)

Z.t2<-B3[(rowsB3Ajuste+1):numrow,(colsB3Ajuste+1):dim(B3)[2]]%*%solve(MatF3)

Z.b22<-kronecker(Z.t2,One.s)

#### Z AMPLIADA PARTE INTERACCION


Z.b31<-kronecker(xt0,Z.s)
Z.b41<-kronecker(Z.t1,xhat)
Z.b51<-kronecker(Z.t1,Z.s)

Z.b61<-kronecker(Z.t2,Bs)



Z0<-cbind(Z.b10,Z.b21,Z.b22,Z.b31,Z.b41,Z.b51,Z.b61)




##CALCULAMOS LAS PREDICCIONES DE LOS RIESGOS

source("dumpsdata/Tr2t")
source("dumpsdata/Tr2st")
source("dumpsdata/a")
source("dumpsdata/b")
source("dumpsdata/G")

#### Seleccion de los coeficientes aleatorios de cada componente

                                                     
bspa<-b[1:dim(Z.b10)[2]]

bt<-b[(dim(Z.b10)[2]+1):(dim(Z.b10)[2]+dim(Z.b21)[2])]

bst<-b[(dim(Z.b10)[2]+dim(Z.b21)[2]+1):dim(G)[2]]


up.ampt<-cbind(E3%*%Tr2t%*%bt)

Bs<-Rten(B2,B1)
Is1<-diag(1,nrow=dim(B1)[2],ncol=dim(B1)[2])
Is2<-diag(1,nrow=dim(B2)[2],ncol=dim(B2)[2])
Is<-kronecker(Is2,Is1)



if(dim(t(as.matrix(E3)))[1]==1){
E3IT<-kronecker(t(E3),Is)%*%Tr2st} else

if(dim(t(as.matrix(E3)))[1]>1){
E3IT<-kronecker(E3,Is)%*%Tr2st
}

up.ampst<-cbind(-E3IT%*%bst)




risksp.amp<-exp(X0%*%a + Z.b10%*%bspa+ Z.b21%*%bt 
+ Z.b22%*%up.ampt+ cbind(Z.b31,Z.b41,Z.b51)%*%bst +  Z.b61%*%up.ampst)


risks.array<-NULL
for(i in 1:p) {
lim1<-50*(i-1)+1
lim2<-50*i
risks.array<-cbind(risks.array,risksp.amp[c(lim1:lim2)])
 }
 
source("dumpsdata/risks")

risks.com<-cbind(risks,risks.array)
riskspredict<-risks.com

dump("riskspredict","dumpsdataPred/riskspredict") 
 
 
 
 
 
pdf("Plot_RisksSoloTrendPrediction.pdf",onefile=TRUE,width=23, height=13)
for(i in 1:50){

plot(c(1975:(2008+p)),riskspredict[i,],type="l")
}

dev.off()
 
