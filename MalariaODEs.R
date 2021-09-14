require(deSolve)

SEA<-function(t,y,p){
  PSU = y[1]
  PSD = y[2]
  #PS = y[3]
  PRU = y[3]
  PRD = y[4]
  #PR = y[6]
  U = y[5]
  D = y[6]
  with(as.list(p),{
    N = y[5]+y[6]
    m = (PSU+PSD+PRU+PRD)/N
    #b = -A*cos(L*sin(2*pi*t)-(2*pi*t))+A
    dPSU.dt = -(b*(PSU/N)*ws*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PSD/N)*ws*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muSU)*(PSU)
    #dPSU.dt = -(b*(PSU/N)*ws*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PSD/N)*ws*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) + (b*(PSU/N)*ws*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muSU*PSU)
    dPSD.dt = -(b*(PSD/N)*ws*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PSU/N)*ws*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muSD)*(PSD)
    #dPSD.dt = -(b*(PSD/N)*ws*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PSU/N)*ws*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) + (b*(PSD/N)*ws*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K))- (muSD*PSD)
    #dPS.dt = dPSU.dt + dPSD.dt
    dPRU.dt = -(b*(PRU/N)*wr*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PRD/N)*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muRU)*(PRU)
    #dPRU.dt = -(b*(PRU/N)*wr*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PRD/N)*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) + (b*(PRU/N)*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muRU*PRU)
    dPRD.dt = -(b*(PRD/N)*wr*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PRU/N)*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muRD)*(PRD)
    #dPRD.dt = -(b*(PRD/N)*wr*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (b*(PRU/N)*D*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) + (b*(PRD/N)*U*(1/K)*(((PSU+PSD+PRU+PRD)/N)-K)) - (muRD*PRD)
    #dPR.dt = dPRU.dt + dPRD.dt
    dU.dt = -b*((1-exp(-m))*(((PSD+PSU+((PRD+PRU)*(wr))))/(PSD+PSU+PRD+PRU))*(U)*d) + (D/Tau)
    dD.dt = b*((1-exp(-m))*(((PSD+PSU+((PRD+PRU)*(wr))))/(PSD+PSU+PRD+PRU))*(U)*d) - (D/Tau)
    return(list(c(dPSU.dt,dPSD.dt,dPRU.dt,dPRD.dt,dU.dt,dD.dt)))
  })
}

b = 0.2239


L=1
A = 0.70
per = 1
ws = 1
wr = 1-0.6
#muSU = 1
muSU=0.05
#muSD = 20
muSD=0.99
#muRD = 1
muRD = 0.05
#muRU = 1
muRU = 0.05
Tau = 20*365
d = 0.3
#m = 3
K=10
p2 = list(b=b,ws=ws,wr=wr,muSU=muSU,muSD=muSD,muRD=muRD,muRU=muRU,Tau=Tau,d=d,m=m)

#t2 = c(1,5,10,20)
t2 = seq(from=0,to=500,by=1/365) 


N0 = c(500,500,0,1,3000,3000)
out2 = ode(y=N0,times=t2,func=SEA,parms=p2)

plot(out2[,1],(out2[,2]),type="l", xlab="Time", ylab="Population", main = "Peak Transmission Level=0.40")
lines(out2[,1],(out2[,3]),col="RED")
legend(x=1500,y=20000,legend=c("PSU","PSD"),bty="n",lwd = c(2,1),col=c("BLACK","RED"))

plot(out2[,1],(out2[,4]),type="l",ylim=c(0,max(out2[,4])), xlab = "Time", ylab="Population", main = "Peak Transmission Level=0.40")
lines(out2[,1],(out2[,5]),col="RED")
legend(x=2000,y=7000,legend=c("PRU","PRD"),bty="n",lwd = c(2,1),col=c("BLACK","RED"))

#plot(out2[,1],(out2[,6]),type="l",ylim=c(0,12500), xlab="Time", ylab="Population", main = "Peak Transmission Level=0.40")
#lines(out2[,1],(out2[,7]),col="RED")
#legend(x=0,y=8000,legend=c("U","D"),bty="n",lwd=c(2,1),col=c("BLACK","RED"))

out3eq = vector(mode = "numeric", length = 21)
out3eq2 = vector(mode = "numeric", length = 21)
out3eqrat = vector(mode= "numeric", length = 21)
out4eq = vector(mode = "numeric", length = 13)
out4eq2 = vector(mode="numeric", length=13)
out4eqrat = vector(mode="numeric", length = 13)
dvalues = c(0.10,0.15,0.20,0.25,0.27,0.30,0.32,0.35,0.37,0.40,0.45,0.50,0.55,0.60,0.70,0.80,0.90)
bvalues = c(0.10,0.15,0.20,0.25,0.27,0.30,0.32,0.35,0.37,0.40,0.45,0.50,0.55,0.60)
N0 = c(500,500,0,1,3000,3000)
wrvalues = c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1)
Avalues = c(0.10,0.15,0.20,0.25,0.27,0.30,0.32,0.35,0.37,0.40,0.45,0.50,0.55,0.60,0.70,0.80,0.90)
Lvalues= c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1)
for (i in 1:21){
  #wr = wrvalues[i]
  #A=Avalues[i]
  L=Lvalues[i]
  p3 = list(b=b,ws=ws,wr=wr,muSU=muSU,muSD=muSD,muRD=muRD,muRU=muRU,Tau=Tau,d=d)
  out3 = ode(y=N0,times=t2,func=SEA,parms=p3)
  #out3eq[i] = out3[3999,2] + out3[3999,3]+out3[3999,4] + out3[3999,5]
  out3eq2[i] = (out3[3999,4]+out3[3999,5])/(out3[3999,2] + out3[3999,3] + out3[3999,4] + out3[3999,5])
  #out3eqrat[i] = out3[3999,6]/out3[3999,7]
}
#print(out3eqrat)
plot(Lvalues,out3eq2,xlab="Relative Length of Dry Season",ylab="Frequency of Resistance", main="d=0.30, A=0.70")
plot(Lvalues,out3eq,xlab="Relative Length of Dry Season",ylab="Total Parasite Population", main="d=0.30, A=0.70")
#print(out3eq2)
#plot(bvalues,out3eqrat,xlab="b",ylab="U/D")
#plot(wrvalues,out3eq,xlab="1-S",ylab="Total Parasite Population",main="d=0.30, b=0.15")
#plot(Lvalues,out3eq2,xlab="Relative Length of Dry Season",ylab="Frequency of Resistance",main="d=0.30, Peak Transmission Level=0.05")
#for (i in 1:13){
  #b = bvalues[i]
  #d = dvalues[i]
#  wr = wrvalues[i]
#  p4 = list(b=b,ws=ws,wr=wr,muSU=muSU,muSD=muSD,muRD=muRD,muRU=muRU,Tau=Tau,d=d)
#  out4 = ode(y=N0,times=t2,func=SEA,parms=p4)
#  out4eq[i] = out4[3999,2] + out4[3999,3] + out4[3999,4] + out4[3999,5]
#  out4eq2[i] = (out4[3999,4]+out4[3999,5])/(out4[3999,2] + out4[3999,3] + out4[3999,4] + out4[3999,5])
#  out4eqrat[i] = out4[3999,6]/out4[3999,7]
#}
#print(out4eqrat)
#plot(wrvalues,out4eq)
#plot(wrvalues,out4eq2)
#plot(wrvalues,out4eqrat)

#mat = matrix(data=0,nrow=17,ncol=17)
#mat2 = matrix(data=0,nrow=17,ncol=17)
#mat3 = matrix(data=0,nrow=17,ncol=17)
#out3eq = vector(mode = "numeric", length = 17)
#out3eqrat = vector(mode= "numeric", length = 17)
#out4eq = vector(mode = "numeric", length = 17)
#out4eqrat = vector(mode="numeric", length = 17)
#dvalues = c(0.10,0.15,0.20,0.25,0.27,0.30,0.32,0.35,0.37,0.40,0.45,0.50,0.55,0.60,0.70,0.80,0.90)
#Avalues = c(0.10,0.15,0.20,0.25,0.27,0.30,0.32,0.35,0.37,0.40,0.45,0.50,0.55,0.60,0.70,0.80,0.90)
#N0 = c(500,500,0,1,3000,3000)
#for (i in 1:17){
#  for (j in 1:17){
#    d = dvalues[i]
#    A = Avalues[j]
#    p2 = list(b=b,ws=ws,wr=wr,muSU=muSU,muSD=muSD,muRD=muRD,muRU=muRU,Tau=Tau,d=d,m=m)
#    out5 = ode(y=N0,times=t2,func=SEA,parms=p2)
#    PTotal = out5[3999,2] + out5[3999,3] + out5[3999,4] + out5[3999,5]
#    mat[i,j] = (out5[3999,4] + out5[3999,5])/PTotal
#    mat2[i,j] = out5[3999,2] + out5[3999,3] + out5[3999,4] + out5[3999,5]
    #mat3[i,j] = out5[3999,6]/out5[3999,7]
#  }
#}
#print(mat)

#library(tcltk)
#library(plot3D)

#install.packages(c("rgl", "car"))

#library("car")


#surf3D(bvalues,dvalues,mat)
#lambda2 = -b(1-exp(-m))(((PSD+PSU+((PRD+PRU)*(wr))))/(PSD+PSU+PRD+PRU))d - (1/Tau)

#condition for D=U: 2b(1-exp(-m))(((PSD+PSU+((PRD+PRU)*(wr))))/(PSD+PSU+PRD+PRU))dTau=1

#U/D~(PSD+PSU+((PRD+PRU)*(wr))))/(PSD+PSU+PRD+PRU)

#(6b) < muSU+muSD ; (PSD+PSU)/N=carrying capacity

FLUC<-function(t,A,L){
  L=1
  A = 0.40
  b = -A*cos(L*sin(t*2*pi)-(t*2*pi))+A
  return(b)
}

integrate(FLUC,0,1)