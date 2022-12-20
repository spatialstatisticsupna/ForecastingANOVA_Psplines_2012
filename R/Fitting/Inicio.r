


################################################################################
################# CALCULO DE LOS CASOS ESPERADOS ###############################

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

setwd("C:/APLICACIONES/Jaione/ForecastingANOVAPspline/Ajuste") 



################################################################################
####SELECCIONAMOS DEL ARCHIVO MORTALIDAD LOS DATOS DE CANCER DE PROSTATA
####PARA 50 PROVINCIAS DE ESPAÑA. 


dat1<-read.table(file="Data/Mortalidad.txt", header = TRUE, sep = ";")
dat1$agno<-rep(0,dim(dat1)[1])
dat1$agno<-dat1$ano+1900
dat<-dat1[dat1$tumorlb=="PROSTATA" & dat1$sex==1 & dat1$prov<=50,]


dat<-dat[,c(3:4,6:7,9:11)]  ###SELECCIONAMOS SOLO LAS COLUMNAS QUE NOS INTERESAN

#### CALCULAMOS EL NÚMERO DE CASOS ESPERADOS PARA CADA PROVINCIA Y AÑO. 

dat<-dat[order(dat$edadgr), ] 
dat<-dat[order(dat$prov,dat$agno), ] 

Dj<-tapply(dat$casos, dat$edadgr,sum)
Nj<-tapply( dat$pob, dat$edadgr,sum)
Rj<-Dj/Nj
Rj.amp<-c(rep(Rj,34*50))
dat$Rj<-Rj.amp
dat$esp<-(dat$pob)*(dat$Rj)

#### GUARDAMOS LOS DATOS 

ProstataH<-aggregate(dat[,c(5,9)],list(prov= dat$prov,agno= dat$agno),sum)
write.table(ProstataH, file = "Data/ProstataH.txt", sep = "\t", col.names =TRUE) 

################################################################################
####EN EL ARCHIVO prostataH QUE ACABAMOS DE GENERAR AÑADIMOS LA LONGITUD  
####Y LATITUD DEL CENTROIDE DE CADA UNA DE LAS AREAS. ESTA INFORMACIÓN SE EXTRAE
####DEL ARCHIVO esp_prov.shp DE LA CARPERTA ZonificacionEspanaProv 


dat<-read.table(file="Data/ProstataH.txt", header = TRUE, sep = "\t")


xSC<- readShapePoly("ZonificacionEspanaProv/esp_prov.shp",IDvar="NAME", 
proj4string=CRS("+proj=longlat +ellps=clrk66"))


### DEL CONJUNTO DE DATOS xSC EXTRAEMOS LA LONGITUD Y LATITUD DE CADA UNA DE
### LAS PROVINCIAS


x<-coordinates(xSC)[,1]
y<-coordinates(xSC)[,2]

x<-rep(x,34)
y<-rep(y,34)

dat$long<-x
dat$lat<-y

## GUARDAMOS EL CONJUNTO DE DATOS QUE CONTIENE, LOS CASOS OBSERVADOS Y ESPERADOS
# EN CADA PROVINCIA Y AÑO. TAMBIÉN LA LONGITUD Y LA LATITUD EL CENTROIDE DE CADA
# PROVINCIA.

ProstH<-dat[order(dat$agno), ] 
write.table(ProstH,"Data/ProstH.txt")

################################################################################
############## CONSTRUCCIÓN DE LAS MATRICES PARA EL MODELO #####################

#######Ejecutar las definiciones para construir un producto tensorial por filas
## y las bases de B-splines. 

bspline<-function(x, xl, xr, ndx, bdeg){
 dx <- (xr-xl)/ndx
 knots <- seq(xl-bdeg*dx, xr+bdeg*dx, by=dx)
 B <- spline.des(knots,x,bdeg+1,0*x)$design
 B
 }

####### Producto tensorial por filas

Rten<-function(X1,X2){
one1<-matrix(1,1,ncol(X1))
one2<-matrix(1,1,ncol(X2))
kronecker(X1,one2)*kronecker(one1,X2)
} 
 

### LEEMOS LOS DATOS

prost<- read.table("Data/ProstH.txt",header=TRUE)###Leemos los datos de Cancer

bdeg<-3  ## SELECCIONAMOS BASES DE ORDEN 3  
pord<-2  ## SELECCIONAMOS PENALIZACIONES DE ORDEN 2 

#### DEFINICIÓN DE LAS VARIABLES EXPLICATIVAS

x1<-coordinates(xSC)[,1]
x1<-(x1-min(x1))/1000


x2<-coordinates(xSC)[,2]
x2<-(x2-min(x2))/1000



x3<-c(1975:2008)
x3<-(x3-min(x3))/1000    

#x3<-((x3-mean(x3))/sqrt(var(x3))) 
#x3<-(x3-mean(x3))/max(x3)

#### SELECCION DE NODOS

##NOTA: ndx indica el número de intervalos internos en el rango de las variables
## por tanto el número de nodos es ndx+1


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


ndx3<-8 #min(ceiling(length(x3)/4),40)#2#
dis3 <- (max(x3)-min(x3))/ndx3
x3l<-min(x3)-dis3*0.05
x3r<-max(x3)+dis3*0.05
dx3 <- (x3r-x3l)/ndx3
knots3<-seq(x3l-bdeg*dx3, x3r+bdeg*dx3, by=dx3)


#### CONSTRUCCIÓN DE LAS BASES DE BSPLINES 

B1<-spline.des(knots1,x1,bdeg+1,0*x1)$design  ##Base marginal longitud
B2<-spline.des(knots2,x2,bdeg+1,0*x2)$design  ##Base marginal latitud
B3<-spline.des(knots3,x3,bdeg+1,0*x3)$design  ##Base marginal tiempo

rowsB3Ajuste<-dim(B3)[1]
colsB3Ajuste<-dim(B3)[2]


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


###### Penalizacion 3

m3=ncol(B3)
D3=diff(diag(m3),differences=pord)
P3.svd=svd(t(D3)%*%D3)
U3s=(P3.svd$u)[,1:(m3-pord)] # le quito dos columnas por los dos autovalores que son 0
U3n=(P3.svd$u)[,(m3-pord+1):m3] # dos columnas de los dos autovalores que son 0
d3=(P3.svd$d)[1:(m3-pord)]
Delta3=diag(d3)


Tr2t<-U3s



Tr2st<-as.matrix(cbind(
kronecker(kronecker(U3n[,2],U2s),U1n),
kronecker(kronecker(U3n[,2],U2n),U1s),
kronecker(kronecker(U3n[,2],U2s),U1s),
kronecker(kronecker(U3s,U2n[,1]),U1n[,2]),
kronecker(kronecker(U3s,U2n[,2]),U1n[,1]),
kronecker(kronecker(U3s,U2n[,2]),U1n[,2]),
kronecker(kronecker(U3s,U2s),U1n),
kronecker(kronecker(U3s,U2n),U1s),
kronecker(kronecker(U3s,U2s),U1s)))


rowsD3Ajuste<-dim(D3)[1]
colsD3Ajuste<-dim(D3)[2]



##CONSTRUCCIÓN DE LAS MATRICES DE EFECTOS FIJOS Y EFECTOS ALEATORIOS 
#CORRESPONDIENTES A LA REPARAMETRIZACIÓN COMO MODELO MIXTOS DEL 
#MODELO DE P-SPLINES TIPO ANOVA.

############ MATRIZ X


X1<-cbind(1,x1)
X2<-cbind(1,x2)
Xt<-cbind(1,x3)

Xs<-Rten(X2,X1)

One.t<-Xt[,1]
One.s<-Xs[,1]

xs<-Xs[,4]
xt<-x3
xhat<-cbind(x1,x2,xs) 



X<-cbind(kronecker(One.t,Xs),kronecker(xt,One.s),kronecker(xt,xhat))

############ MATRIZ Z


Z1<-B1%*%U1s
Z2<-B2%*%U2s
Z3<-B3%*%U3s

Z.s<-cbind(Rten(Z2,X1),Rten(X2,Z1),Rten(Z2,Z1))
Z.t<-Z3

One.t<-Xt[,1]
One.s<-Xs[,1]


Z.b1<-kronecker(One.t,Z.s)
Z.b2<-kronecker(Z.t,One.s) 
 
Z.b3<-kronecker(xt,Z.s)
Z.b4<-kronecker(Z.t,xhat)
Z.b5<-kronecker(Z.t,Z.s)

Z<-cbind(Z.b1,Z.b2,Z.b3,Z.b4,Z.b5)

k1<-m1
k2<-m2
k3<-m3



###GUARDAMOS LAS MATRICES GENERADAS EN LA CARPERTA LLAMADA dumpsdata.

dump("k1",file="dumpsdata/k1")
dump("k2",file="dumpsdata/k2")
dump("k3",file="dumpsdata/k3")


dump("X",file="dumpsdata/X")
dump("Z",file="dumpsdata/Z")

dump("Delta1",file="dumpsdata/Delta1")
dump("Delta2",file="dumpsdata/Delta2")
dump("Delta3",file="dumpsdata/Delta3")

dump("B1",file="dumpsdata/B1")
dump("B2",file="dumpsdata/B2")
dump("B3",file="dumpsdata/B3")


dump("U3n","dumpsdata/U3n")
dump("U3s","dumpsdata/U3s")
dump("Tr2t","dumpsdata/Tr2t")
dump("Tr2st","dumpsdata/Tr2st")

dump("rowsD3Ajuste","dumpsdata/rowsD3Ajuste")
dump("colsD3Ajuste","dumpsdata/colsD3Ajuste")

dump("rowsB3Ajuste","dumpsdata/rowsB3Ajuste")
dump("colsB3Ajuste","dumpsdata/colsB3Ajuste")

