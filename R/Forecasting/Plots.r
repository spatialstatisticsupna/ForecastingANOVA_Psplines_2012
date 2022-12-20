

################################################################################
################################################################################
### SINTAXIS QUE PERMITE DIBUJAR LAS TENDENCIAS Y LOS MAPAS DE LOS RIESGOS #########


source("dumpsdataPred/msecom")
source("dumpsdataPred/riskspredict")

p<-3

totalyears<-length(unique(seq(min(dat$agno),max(dat$agno)+p)))


msecom<-msecom[,1:totalyears]
b.ij2<-log(riskspredict[,1:totalyears])

IlinfS<-exp(b.ij2+qnorm(0.025)*sqrt(msecom))
IlsupS<-exp(b.ij2+qnorm(0.975)*sqrt(msecom))
risk<-riskspredict[,1:totalyears]



Year<-c(1975:(1975+totalyears-1))

#### ORDENAMOS LAS PROVICIAS POR COMUNIDADES AUTONOMAS DE NORTE A SUR.

ord<-c(15,27,32,36,33,39,26,1,20,48,31,22,50,44,9,34,24,49,47,42,
37,5,40,25,17,8,43,12,46,3,30,7,28,10,6,19,45,16,13,2,21,41,11,14,29,23,18,4,35,38)

nombres<-c("Alava","Albacete","Alicante", "Almeria ",
"Avila", "Badajoz","Baleares", "Barcelona", "Burgos", "Caceres", "Cadiz",
"Castellon", "CiudadReal", "Cordoba", "Coruña", "Cuenca", "Girona",
"Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen",
"Leon", "Lleida", "Rioja", "Lugo", "Madrid", "Malaga", "Murcia",
"Navarra", "Ourense", "Asturias", "Palencia","Las Palmas", "Pontevedra", "Salamanca",
 "Santa Cruz de Tenerife","Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
"Valencia", "Valladolid", "Vizcaya", "Zamora", "Zaragoza")


source("dumpsdata/smr")





plot.smr<-function(Year,risk,smr,IlinfS,IlsupS,nombres,aa,bb){
#par(mfrow=c(2,3),cex.axis=1.3)

for(j  in aa:bb){
ii<-ord[j]


plot(c(1975:2008),risk[ii,1:34],type="l",xlab="",ylab="",yaxt="n",xaxt="n",lty=1,lwd=1.2,ylim=c(min(smr),max(smr)),xlim=c(1975,2012))
xx1 <- c(c(1975:2008, rev(c(1975:2008))))
yy <- c(IlinfS[ii,1:34], rev(IlsupS[ii,1:34]))
polygon(xx1, yy, col="#D1D1D1", border = NA)

xx2 <- c(c(2008:(2008+p)), rev(c(2008:(2008+p))))
yy2 <- c(IlinfS[ii,34:(34+p)], rev(IlsupS[ii,34:(34+p)]))
polygon(xx2, yy2, col="#949494", border = NA)


lines(c(1975:2008),risk[ii,1:34],type="l",cex=2)
lines(c(2008:(2008+p)),risk[ii,34:(34+p)],type="l",cex=2)

axis(1,at=c(1975:(2008+p)), labels=FALSE, tick = TRUE,las=2,cex.axis=2.25)

text(c(min(Year),1992,2008), par("usr")[3]-0.03,
     labels = c("1975", "1992", "2008"),
     srt = 1, pos = 1, xpd = TRUE,cex=3)


#axis(1,at=c(min(Year),1992,2008),las=2,cex.axis=2.25)
axis(2,at=c(0.5,1,1.5),las=3,cex.axis=3)

lines(Year[1:34],smr[ii,],type="o",col=rgb(102,102,102,max = 255),lwd=1.3)
lines(Year[1:34],smr[ii,],type="l",col=rgb(102,102,102,max = 255),cex=1)
lines(Year,IlinfS[ii,1:(34+p)],type="l",main=nombres[ii],cex=3,lwd=1.6)
lines(Year,IlsupS[ii,1:(34+p)],type="l",main=nombres[ii],cex=3,lwd=1.6)
legend(1970,1.9,legend=nombres[ii],cex=4,bty = "n")

abline(h=1,lty=2,lwd=1.5)
abline(v=2008,lty=2,lwd=1.5)
}
}


plot.smr.col<-function(Year,risk,smr,IlinfS,IlsupS,nombres,aa,bb){
#par(mfrow=c(2,3),cex.axis=1.3)

for(j  in aa:bb){
ii<-ord[j]


plot(c(1975:2008),risk[ii,1:34],type="l",xlab="",ylab="",yaxt="n",xaxt="n",lty=1,lwd=1.2,ylim=c(min(smr),max(smr)),xlim=c(1975,2012))
xx1 <- c(c(1975:2008, rev(c(1975:2008))))
yy <- c(IlinfS[ii,1:34], rev(IlsupS[ii,1:34]))
polygon(xx1, yy, col="#D1D1D1", border = NA)

xx2 <- c(c(2008:(2008+p)), rev(c(2008:(2008+p))))
yy2 <- c(IlinfS[ii,34:(34+p)], rev(IlsupS[ii,34:(34+p)]))
polygon(xx2, yy2, col="#3A5FCD", border = NA)
                      

lines(c(1975:2008),risk[ii,1:34],type="l",cex=2)
lines(c(2008:(2008+p)),risk[ii,34:(34+p)],type="l",cex=2)

axis(1,at=c(1975:(2008+p)), labels=FALSE, tick = TRUE,las=2,cex.axis=2.25)

text(c(min(Year),1992,2008), par("usr")[3]-0.03,
     labels = c("1975", "1992", "2008"),
     srt = 1, pos = 1, xpd = TRUE,cex=3)


#axis(1,at=c(min(Year),1992,2008),las=2,cex.axis=2.25)
axis(2,at=c(0.5,1,1.5),las=3,cex.axis=3)


lines(Year[1:34],smr[ii,],type="o",col=rgb(102,102,102,max = 255),lwd=1.3)
lines(Year[1:34],smr[ii,],type="l",col=rgb(102,102,102,max = 255),cex=1)
lines(Year,IlinfS[ii,1:(34+p)],type="l",main=nombres[ii],cex=3,lwd=1.6)
lines(Year,IlsupS[ii,1:(34+p)],type="l",main=nombres[ii],cex=3,lwd=1.6)
legend(1975,1.9,legend=nombres[ii],cex=4,bty = "n")

abline(h=1,lty=2,lwd=1.5)
abline(v=2008,lty=2,lwd=1.5)
}
}




#### ESTA SINTAXIS PERMITE GUARDAR DIRECTAMENTE LAS FIGURAS EN UN PDF.

pdf("Plot_RisksTrendPrediction.pdf",onefile=TRUE,width=23, height=13)
par(mfrow=c(2,3),cex.axis=1.3)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=1,bb=6)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=7,bb=12)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=13,bb=18)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=19,bb=24)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=25,bb=30)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=31,bb=36)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=37,bb=42)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=43,bb=48)
plot.smr(Year, risk,smr,IlinfS,IlsupS,nombres,aa=49,bb=50)
dev.off()


#### FIGURAS EN COLOR
#### ESTA SINTAXIS PERMITE GUARDAR DIRECTAMENTE LAS FIGURAS EN UN PDF.

pdf("Plot_RisksTrendPredictionCol.pdf",onefile=TRUE,width=20, height=13)
par(mfrow=c(1,1),cex.axis=1.3)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=1,bb=6)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=7,bb=12)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=13,bb=18)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=19,bb=24)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=25,bb=30)#
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=31,bb=36)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=37,bb=42)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=43,bb=48)
plot.smr.col(Year, risk,smr,IlinfS,IlsupS,nombres,aa=49,bb=50)
dev.off()



################################################################################
################################################################################
########### SINTAXIS QUE PERMITE DIBUJAR LOS MAPAS DE LOS RIESGOS ##############







xName <- readShapePoly("ZonificacionEspanaProv/esp_prov.shp",IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))

xNameTras <- xName[xName$ESP_PROV_I !=51 & xName$ESP_PROV_I !=52,]

#Ïndice de los polinomos a trasladar
idx<-c(35, 38)


#plot(xName[idx,])

#Traslación que aplicamos a long/lat
#shft<-c(6, 7) #NUEVO
#shft<-c(5, 8) #NUEVO
shft<-c(18, 8) ### Ponerlas debajo de Baleares

#Bucle sobre las provincias a trasladar
for(prov.index in idx)
{
    prov<-slot(((xName [prov.index,])), "polygons")
    prov.polygons <-slot(prov[[1]], "Polygons")

    new.polygons <-as.list(rep(NA, length(prov.polygons)))
    for(pol.index in 1:length(prov.polygons))
    {
        xx<-  slot(prov.polygons[pol.index][[1]], "coords")
        xx<-xx+matrix(shft, byrow=TRUE, ncol=2, nrow=nrow(xx))
        #points(xx)
        new.polygons[pol.index]<-Polygon(xx)
    }

    xNameTras@polygons[prov.index]<-Polygons(new.polygons, ID=slot(prov[[1]], "ID"))
}


#Aumentar un poco los márgenes del bounding box del mapa total
bbox.xNameTras <- bbox(xNameTras[(1:50),])                   #
bbox.xNameTras[c(1,2)] <- bbox.xNameTras[c(1,2)] - 0.1       # Bounding box del mapa total
bbox.xNameTras[c(3,4)] <- bbox.xNameTras[c(3,4)] + 0.1       #
xNameTras@bbox<-bbox.xNameTras

### Determinar bounding box de las islas
bbox.islas <- bbox(xNameTras[(1:50)[idx],])

### Aumentar un poco los márgenes del rectángulo
bbox.islas[c(1,2)] <- bbox.islas[c(1,2)] - 0.1
bbox.islas[c(3,4)] <- bbox.islas[c(3,4)] + 0.1

### Construir rectángulo para las islas.
rect.islas <- rbind(c(bbox.islas[1],bbox.islas[2]),c(bbox.islas[1],bbox.islas[4]),c(bbox.islas[3],bbox.islas[4]),c(bbox.islas[3],bbox.islas[2]),c(bbox.islas[1],bbox.islas[2]))

### Construir objetos Polygon y Polygons con el rectángulo para las islas.
polygon.islas <- Polygon(rect.islas)
polygons.islas <- Polygons(list(polygon.islas), "Cuadro")

### Incluir el nuevo polígono en xNameTras
num.rect <- length(xNameTras@polygons)+1
xNameTras@polygons[[num.rect]] <- polygons.islas

### Determinar el orden en que se dibujará el nuevo polígono, así se dibujará primero.
xNameTras@plotOrder <- c(as.integer(num.rect),xNameTras@plotOrder)


xtodas <- xNameTras

# Incluir un nuevo nivel para el rectángulo de canarias
levels(xtodas$NAME)<- cbind(xtodas$NAME, " ")

# Añadir la nueva fila al data frame correspondiente al rectángulo de canarias, con valores no válidos
df <- xtodas@data
df[num.rect,] <- c(num.rect, as.factor(NA), as.factor(" "), rep(NA, dim(df)[2]-3))
xtodas@data <- df


xSC<-xtodas


risks<-riskspredict

xSC$risk1<-c(risks[,1],NA)
xSC$risk2<-c(risks[,2],NA)
xSC$risk3<-c(risks[,3],NA)
xSC$risk4<-c(risks[,4],NA)
xSC$risk5<-c(risks[,5],NA)
xSC$risk6<-c(risks[,6],NA)
xSC$risk7<-c(risks[,7],NA)
xSC$risk8<-c(risks[,8],NA)
xSC$risk9<-c(risks[,9],NA)
xSC$risk10<-c(risks[,10],NA)
xSC$risk11<-c(risks[,11],NA)
xSC$risk12<-c(risks[,12],NA)
xSC$risk13<-c(risks[,13],NA)
xSC$risk14<-c(risks[,14],NA)
xSC$risk15<-c(risks[,15],NA)
xSC$risk16<-c(risks[,16],NA)
xSC$risk17<-c(risks[,17],NA)
xSC$risk18<-c(risks[,18],NA)
xSC$risk19<-c(risks[,19],NA)
xSC$risk20<-c(risks[,20],NA)
xSC$risk21<-c(risks[,21],NA)
xSC$risk22<-c(risks[,22],NA)
xSC$risk23<-c(risks[,23],NA)
xSC$risk24<-c(risks[,24],NA)
xSC$risk25<-c(risks[,25],NA)
xSC$risk26<-c(risks[,26],NA)
xSC$risk27<-c(risks[,27],NA)
xSC$risk28<-c(risks[,29],NA)
xSC$risk29<-c(risks[,29],NA)
xSC$risk30<-c(risks[,30],NA)
xSC$risk31<-c(risks[,31],NA)
xSC$risk32<-c(risks[,32],NA)
xSC$risk33<-c(risks[,33],NA)
xSC$risk34<-c(risks[,34],NA)
xSC$risk35<-c(risks[,35],NA)
xSC$risk36<-c(risks[,36],NA)
xSC$risk37<-c(risks[,37],NA)
xSC$risk38<-c(risks[,38],NA)
xSC$risk39<-c(risks[,39],NA)




#### ESTA SINTAXIS PERMITE GUARDAR DIRECTAMENTE LOS MAPAS EN UN PDF.

pdf("Plots_RisksMapsPredict.pdf")
spplot(xSC,zcol=c("risk1","risk2","risk3","risk4","risk5","risk6","risk7","risk8","risk9","risk10",
"risk11","risk12","risk13","risk14","risk15","risk16","risk17","risk18","risk19","risk20",
"risk21","risk22","risk23","risk24","risk25","risk26","risk27","risk28","risk29","risk30","risk31","risk32",
"risk33","risk34","risk35","risk36","risk37"),
names.attr=c("1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985",
"1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998",
"1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011"),
as.table=TRUE,col.regions=brewer.pal(5, "YlOrRd"),layout=c(6,7),at=c(.60,.80,1.00,1.20,1.4,1.6))
dev.off()
