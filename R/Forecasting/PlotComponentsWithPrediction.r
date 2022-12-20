


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


risksp.amp<-exp(X0%*%a + Z.b10%*%bspa+ Z.b21%*%bt
+ Z.b22%*%up.ampt+ cbind(Z.b31,Z.b41,Z.b51)%*%bst +  Z.b61%*%up.ampst)

X0<-cbind(kronecker(One.t0,Xs),kronecker(xt0,One.s),kronecker(xt0,xhat)) ### AMPLIADA

#################### f(t)




Xft<-X[,5]
Zft<-kronecker(Z.t,One.s)

b.ij2t<-cbind(Xft)*a[5] + Zft%*%b[166:174]

radjSt<-exp(b.ij2t)

riskst<-matrix(radjSt,nrow=50,ncol=34)

risksp.ampt<-matrix(exp(cbind(kronecker(xt0,One.s))*a[5]+Z.b21%*%bt+ Z.b22%*%up.ampt),nrow=50,ncol=9)

comptemp<-cbind(riskst,risksp.ampt)


#################### f(x1,x2)


Xfx12<-X[,1:4]

Z.b1<-kronecker(One.t,Z.s)
Zfx12<-Z.b1


risksspa<-matrix(exp(kronecker(One.t,Xs)%*%a[1:4] + Z.b1%*%bspa),nrow=50,ncol=34)


risksp.amps<-matrix(exp(kronecker(One.t0,Xs)%*%a[1:4]+Z.b10%*%bspa),nrow=50,ncol=9)

compspa<-cbind(risksspa,risksp.amps)

#################### f(x1,x2,xt)


Z.b3<-kronecker(xt,Z.s)
Z.b4<-kronecker(Z.t,xhat)
Z.b5<-kronecker(Z.t,Z.s)


Xfx12t<-X[,6:8]
Zfx12t<-cbind(Z.b3,Z.b4,Z.b5)

b.ij2x12t<-Xfx12t%*%a[6:8] + Zfx12t%*%b[175:1851]

radjSx12t<-exp(b.ij2x12t)

risksx12t<-matrix(radjSx12t,nrow=50,ncol=34)


risksp.ampst<-matrix(exp(kronecker(xt0,xhat)%*%a[6:8]+cbind(Z.b31,Z.b41,Z.b51)%*%bst +  Z.b61%*%up.ampst),nrow=50,ncol=9)

compst<-cbind(risksx12t,risksp.ampst)






risksp.amp<-matrix(exp(X0%*%a + Z.b10%*%bspa+ Z.b21%*%bt
+ Z.b22%*%up.ampt+ cbind(Z.b31,Z.b41,Z.b51)%*%bst +  Z.b61%*%up.ampst),nrow=50,ncol=9)



riskstot<-cbind(risks,risksp.amp)
#####COMPONENTE TEMPORAL




source("dumpsdata/smr")


######### RIESGOS PARA 6 COMUNIDADES f(t), f(x1,x2), f(x1,x2,t) y f-ANOVA

######### RIESGOS PARA 6 COMUNIDADES f(t), f(x1,x2), f(x1,x2,t) y f-ANOVA

times<-c(1975:2011)
plot.smr<-function(times,smr,comptemp,compspa,compst,riskstot,nombres,aa,bb){
for(ii in aa:bb){
#plot(time,smrs[ii,],type="l",main=nombres[ii],ylab="",xaxt="n",lty=1,lwd=1,ylim=c(min(smrs),max(smrs)))
plot(times,riskstot[orden.prov[ii],1:37],type="l",ylab="",xlab="",xaxt="n",lty=1,lwd=1,ylim=c(min(smr),max(smr)))
 lines(times[1:34],smr[orden.prov[ii],],type="l",lwd=1,lty = 1)
legend(1972,2,legend=nombres[orden.prov[ii]],cex = 3,bty = "n")
axis(1,las=3,cex.axis=2,xaxp=c(min(times), max(times), 12))

lines(times,riskstot[orden.prov[ii],1:37],lty = 1, lwd=4,col="#000000")
lines(times,comptemp[orden.prov[ii],1:37],lty = 5, lwd=3,col="#3A5FCD")
lines(times,compspa[orden.prov[ii],1:37], lty = 1, lwd=1,col="#CD0000")
lines(times,compst[orden.prov[ii],1:37],lty = 6   ,lwd=3,col="#228B22")


#pares=seq(1,37, by = 2)
#points(times[pares],riskstot[orden.prov[ii],pares], pch=1,cex=3,col="#000000")
#points(times[pares],comptemp[orden.prov[ii],pares], pch=2,cex=3,col="#3A5FCD")
#points(times[pares],compspa[orden.prov[ii],pares],  pch=3,cex=3,col="#CD0000")
#points(times[pares],compst[orden.prov[ii],pares],  pch=4,cex=3,col="#228B22")


legend(1996,2,legend=c("Total risk",expression(italic(f[t](x[t]))),expression(italic(f[s](x[1],x[2]))), expression(italic(f[(st)](x[1],x[2],x[t])))),
cex =2,lty = c(1,5,1,6),bty = "n",lwd=c(4,3,1,3),col=c("#000000","#3A5FCD","#CD0000","#228B22"))
abline(h=1,lty=2,lwd=1.5)
abline(v=2008,lty=2,lwd=1.5)
}
}

pdf("PlotsComponenetsColorPruebaTodas.pdf",onefile=TRUE,width=20, height=15)
par(mfrow=c(2,3))
for(i in 1:50 ){
plot.smr(times,smr,comptemp,compspa,compst,riskstot,nombres,aa=i,bb= i) }
dev.off()





pdf("PlotsComponenetsColorPrueba2.pdf",onefile=TRUE,width=20, height=15)
par(mfrow=c(2,3))
for(i in c(1,11,17,28,41,49)){
plot.smr(times,smr,comptemp,compspa,compst,riskstot,nombres,aa=i,bb= i) }
dev.off()  


selec=c(which( nombres[orden.prov]==c("Lugo")),
which( nombres[orden.prov]==c("Zaragoza")),
which( nombres[orden.prov]==c("Tarragona")),
which( nombres[orden.prov]==c("Valladolid")),
which( nombres[orden.prov]==c("Valencia")),
which( nombres[orden.prov]==c("Malaga")))

pdf("PlotsComponenetsWithPrediction4.pdf",onefile=TRUE,width=20, height=15)
par(mfrow=c(2,3))
for(i in selec){
plot.smr(times,smr,comptemp,compspa,compst,riskstot,nombres,aa=i,bb= i) }
dev.off()


