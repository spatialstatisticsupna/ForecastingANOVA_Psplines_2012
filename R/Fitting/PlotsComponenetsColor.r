

Z.b1<-kronecker(One.t,Z.s)
Z.b2<-kronecker(Z.t,One.s)  
Z.b3<-kronecker(xt,Z.s)
Z.b4<-kronecker(Z.t,xhat)
Z.b5<-kronecker(Z.t,Z.s)

Z<-cbind(Z.b1,Z.b2,Z.b3,Z.b4,Z.b5)


source("dumpsdata/a")
source("dumpsdata/b")
source("dumpsdata/risks")



times<-c(1975:2008)



nombres<-c("Alava","Albacete","Alicante", "Almeria ",
"Avila", "Badajoz","Baleares", "Barcelona", "Burgos", "Caceres", "Cadiz",
"Castellon", "CiudadReal", "Cordoba", "Coruña", "Cuenca", "Girona",
"Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen",
"Leon", "Lleida", "Rioja", "Lugo", "Madrid", "Malaga", "Murcia",
"Navarra", "Ourense", "Asturias", "Palencia","Las Palmas", "Pontevedra", "Salamanca",
 "Santa Cruz de Tenerife","Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
"Valencia", "Valladolid", "Vizcaya", "Zamora", "Zaragoza")


orden.prov<-c(15,27,32,36,33,39,26,1,20,48,31,22,50,44,25,17,8,43,24,34,9,49,47,42,
              37,5,40,28,10,6,19,45,16,13,2, 12,46,3,30,7,21,41,11,14,29,23,18,4,35,38)

#################### f(x1,x2)


Xfx12<-X[,1:4]
Zfx12<-Z.b1
                                                     
                                                     
b.ij2x12<-Xfx12%*%a[1:4] + Zfx12%*%b[1:165]

radjSx12<-exp(b.ij2x12)

risksx12<-matrix(radjSx12,nrow=50,ncol=34)



#################### f(t)

Xft<-X[,5]
Zft<-Z.b2

b.ij2t<-cbind(Xft)*a[5] + Zft%*%b[166:174]

radjSt<-exp(b.ij2t)

riskst<-matrix(radjSt,nrow=50,ncol=34)


#################### f(x1,x2,xt)


Xfx12t<-X[,6:8]
Zfx12t<-cbind(Z.b3,Z.b4,Z.b5)

b.ij2x12t<-Xfx12t%*%a[6:8] + Zfx12t%*%b[175:1851]

radjSx12t<-exp(b.ij2x12t)

risksx12t<-matrix(radjSx12t,nrow=50,ncol=34)


risks<-matrix(radjS,nrow=50,ncol=34)



source("dumpsdata/smr")


######### RIESGOS PARA 6 COMUNIDADES f(t), f(x1,x2), f(x1,x2,t) y f-ANOVA

plot.smr<-function(time,smrs,riskT,risk,risk1,risk2,nombres,aa,bb){
for(ii in aa:bb){
#plot(time,smrs[ii,],type="l",main=nombres[ii],ylab="",xaxt="n",lty=1,lwd=1,ylim=c(min(smrs),max(smrs)))
plot(time,smrs[orden.prov[ii],],type="l",ylab="",xlab="",xaxt="n",lty=1,lwd=1,ylim=c(min(smrs),max(smrs)))
legend(1972,2,legend=nombres[orden.prov[ii]],cex = 3,bty = "n")
axis(1,las=3,cex.axis=2,xaxp=c(min(time), max(time), 11))
lines(time,riskT[orden.prov[ii],],type="l",lwd=3,lty = 1)
lines(time,risk[orden.prov[ii],],lty = 5, lwd=3,col="#3A5FCD")
lines(time,risk1[orden.prov[ii],], lwd=3,col="#CD0000")
lines(time,risk2[orden.prov[ii],],lty = 6,lwd=3,col="#228B22")
legend(1996,2,legend=c("Total risk",expression(italic(f[t](x[t]))),expression(italic(f[s](x[1],x[2]))), expression(italic(f[(st)](x[1],x[2],x[t])))),
cex =2,lty = c(1,5,1,6),bty = "n",lwd=3,col=c("#000000","#3A5FCD","#CD0000","#228B22"))
abline(h=1,lty=2,lwd=1.5)
}
}





pdf("PlotsComponenetsColor.pdf",onefile=TRUE,width=20, height=15)
par(mfrow=c(2,3))
for(i in c(1,11,17,28,41,49)){
plot.smr(times,smr,risks,riskst,risksx12,risksx12t,nombres,aa=i,bb= i) }
dev.off()


