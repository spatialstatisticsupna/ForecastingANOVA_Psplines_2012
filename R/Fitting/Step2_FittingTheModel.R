

################################################################################
###### Ajustemos el modelo de P-splines


rm(list=ls())


library(splines)
library(mgcv)
library(nlme)
library(MASS)
library(mvtnorm)
library(maptools)
library(RColorBrewer)
library(splancs)
library(spdep)



setwd(".../Fitting") 


################################################################################
####SELECCIONAMOS LOS DATOS DE CANCER DE PROSTATA PARA 50 PROVINCIAS DE ESPAÑA 

prost<- read.table("Data/ProstH.txt",header=TRUE)###Leemos los datos de Cancer 

y<-prost$casos

prov<-prost$prov

l<-length(prov)

n<-length(y)

esp<-prost$esp

offset1<-log(esp)

### LEEMOS LAS MATRICES QUE HEMOS GENERADO EN EL ARCHIVO Inicio.r Y QUE HAN
# SIDO GUARDADAS EN LA CARPETA dumpsdata.  


source("dumpsdata/Z")
source("dumpsdata/X")

source("dumpsdata/Delta1")
source("dumpsdata/Delta2")
source("dumpsdata/Delta3")


source("dumpsdata/k1")
source("dumpsdata/k2")
source("dumpsdata/k3")



### DEFINIMOS VALORES INICIALES PARA LOS PARAMETROS DEL MODELO 

lambda1<-1
lambda2<-2
lambda3<-5

tau1<-6
tau2<-0.8
tau3<-50

a<-as.vector(array(0,dim(X)[2]))
b<-as.vector(array(0,dim(Z)[2]))

param<-c(lambda1,lambda2,lambda3,tau1,tau2,tau3)

oldparam<-c(a,b,lambda1,lambda2,lambda3,tau1,tau2,tau3)

conv<-1

counter<-0


### EJECUTAMOS EL PQL. ESTE ARCHIVO PERMITE EL AJUSTE DEL MODELO.  
### TARDA UN RATO EN AJUSTARSE, VARIAS HORAS. 

memory.size(4000)
date()
source("PQLs.r")   
date()



### CALCULAMOS LOS INTERVALOS DE CONFIANZA


source("msespline.r")


### GENERAMOS LAS FIGURAS DE LAS TENDENCIAS Y LOS MAPAS DE LOS RIESGOS 

source("Plots.r")



#### DIBUJAR LOS DISTITOS COMPONENTES DEL MODELO


source("PlotsComponenetsColor.r")


### GUARDAMOS LOS RESULTADOS DEL MODELO. CAMBIAR SI SE DESEA
## LA FECHA DEL AJUSTE EN Date PARA QUE LOS DISTINTOS AJUSTES QUEDEN GUARDADOS


#####GUARDAMOS ELEMENTOS QUE VAMOS A NECESITAR EN LA PREDICCION
dump("a","dumpsdata/a")
dump("b","dumpsdata/b")
dump("G","dumpsdata/G")


#####GUARDAMOS MAS ELEMENTOS QUE VAMOS A NECESITAR EN LA PREDICCION
dump("dG1","dumpsdata/dG1")
dump("dG2","dumpsdata/dG2")
dump("dG3","dumpsdata/dG3")
dump("dG3","dumpsdata/dG3")
dump("dG4","dumpsdata/dG4")
dump("dG5","dumpsdata/dG5")
dump("dG6","dumpsdata/dG6")



dump("dV1","dumpsdata/dV1")
dump("dV2","dumpsdata/dV2")
dump("dV3","dumpsdata/dV3")
dump("dV3","dumpsdata/dV3")
dump("dV4","dumpsdata/dV4")
dump("dV5","dumpsdata/dV5")
dump("dV6","dumpsdata/dV6")



dump("V","dumpsdata/V")
dump("V.inv","dumpsdata/V.inv")
dump("XVX.inv","dumpsdata/XVX.inv")

dump("P","dumpsdata/P")



dump("F1","dumpsdata/F1")
dump("F2","dumpsdata/F2")
dump("F3","dumpsdata/F3")
dump("F4","dumpsdata/F4")
dump("F5","dumpsdata/F5")
dump("F6","dumpsdata/F6")

dump("F7","dumpsdata/F7")  
dump("F8","dumpsdata/F8")  
dump("F9","dumpsdata/F9")  
dump("F10","dumpsdata/F10")  
dump("F11","dumpsdata/F11")  
 


save.image("ANOVATypePspline_Estimation_22_01_13.RDATA")







