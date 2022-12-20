################################################################################



source("dumpsdata/smr")
source("dumpsdata/IlinfS95")
source("dumpsdata/IlsupS95")
source("dumpsdata/mse")
source("dumpsdata/risks")

################################################################################
################################################################################
########### SINTAXIS QUE PERMITE DIBUJAR LAS TENDENCIAS DE LOS RIESGOS #########



times<-c(1975:2008)



nombres<-c("Alava","Albacete","Alicante", "Almeria ",
"Avila", "Badajoz","Baleares", "Barcelona", "Burgos", "Caceres", "Cadiz",
"Castellon", "CiudadReal", "Cordoba", "Coruña", "Cuenca", "Girona",
"Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen",
"Leon", "Lleida", "Rioja", "Lugo", "Madrid", "Malaga", "Murcia",
"Navarra", "Ourense", "Asturias", "Palencia","Las Palmas", "Pontevedra", "Salamanca",
 "Santa Cruz de Tenerife","Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
"Valencia", "Valladolid", "Vizcaya", "Zamora", "Zaragoza")

#### ORDENAMOS LAS PROVICIAS POR COMUNIDADES AUTONOMAS DE NORTE A SUR.

orden.prov<-c(15,27,32,36,33,39,26,1,20,48,31,22,50,44,25,17,8,43,24,34,9,49,47,42,
              37,5,40,28,10,6,19,45,16,13,2, 12,46,3,30,7,21,41,11,14,29,23,18,4,35,38)



plot.smr<-function(time,smr,risk,int.infb95,int.supb95,nombres,aa,bb){
par(mfrow=c(1,1),cex.axis=2)
for(ii in aa:bb){
#plot(time,smrs[ii,],type="l",main=nombres[ii],ylab="",xaxt="n",lty=1,lwd=1,ylim=c(min(smrs),max(smrs)))
plot(time,smr[orden.prov[ii],],type="l",ylab="",xlab="",xaxt="n",yaxt="n",lty=1,lwd=1,ylim=c(min(smr),max(smr)))
legend(1975,2,legend=nombres[orden.prov[ii]],cex = 6,bty = "n")


axis(1,at=c(1975:2011), labels=FALSE, tick = TRUE,las=2,cex.axis=2.25,srt=45, adj=0)
#axis(1,at=c(min(times),1992,2008),las=2,cex.axis=2.25 ,srt=45, adj=0)
axis(2,at=c(0.5,1,1.5),las=3,cex.axis=5)

text(c(min(times),1992,2008), par("usr")[3]-0.03,
     labels = c("1975", "1992", "2008"),
     srt = 1, pos = 1, xpd = TRUE,cex=5)


lines(time,risk[orden.prov[ii],],type="l",col=4,cex=1)
lines(time,int.infb95[orden.prov[ii],],type="l",col=4,cex=2)
lines(time,int.supb95[orden.prov[ii],],type="l",col=4,cex=2)
xx <- c(time, rev(time))
yy <- c(int.infb95[orden.prov[ii],], rev(int.supb95[orden.prov[ii],]))
polygon(xx, yy, col=rgb(0.1,0.1,0.1,0.25), border = NA)
abline(h=1,lty=2,lwd=1.5)
}
}



int.infb95<-IlinfS95

int.supb95<-IlsupS95

#### ESTA SINTAXIS PERMITE GUARDAR DIRECTAMENTE LAS FIGURAS EN UN PDF.

pdf("Plot_RisksTrends.pdf",onefile=TRUE,width=20, height=15)
for(i in 1:50)
plot.smr(times,smr,risks,int.infb95,int.supb95,nombres,aa=i  ,bb= i)
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




source("dumpsdata/risks")

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


#### ESTA SINTAXIS PERMITE GUARDAR DIRECTAMENTE LOS MAPAS EN UN PDF.
pdf("Plot_RisksMaps.pdf")

spplot(xSC,zcol=c("risk1","risk2","risk3","risk4","risk5","risk6","risk7","risk8","risk9","risk10",
"risk11","risk12","risk13","risk14","risk15","risk16","risk17","risk18","risk19","risk20",
"risk21","risk22","risk23","risk24","risk25","risk26","risk27","risk28","risk29","risk30","risk31","risk32","risk33","risk34"),
names.attr=c("1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985",
"1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998",
"1999","2000","2001","2002","2003","2004","2005","2006","2007","2008"),
as.table=TRUE,col.regions=brewer.pal(5, "YlOrRd"),layout=c(1,1),at=c(.60,.80,1.00,1.20,1.4,1.6))

dev.off()

