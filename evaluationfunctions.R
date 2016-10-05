# Description of function

plot_forecast<-function(rnr,hnr,syear,smonth,flood_values,ptype="greyshade"){
ndays<-c(1,31,28,31,30,31,30,31,31,30,31,30,fpath='../data/netcdf/')
fd<-c("01/01/","01/02/","01/03/","01/04/","01/05/","01/06/","01/07/","01/08/","01/09/","01/10/","01/11/","01/12/")
ttext<-c(paste(fd,syear,sep=""))
ci<-which(flood_values[,3]==paste(rnr,'_',hnr,sep=''))
xtt<-cumsum(ndays)[smonth:12]
xtt<-xtt-xtt[1]+1
dlabel<-ttext[smonth:12]

nc<-open.nc(paste(fpath,"seasonal_forecast_database_",rnr,'.',hnr,".nc",sep=''),main=flood_values[ci,3])
fyears<-var.get.nc(nc,"Year")
fmonths<-var.get.nc(nc,"Month")
scenario<-var.get.nc(nc,"Scenario")
leadtime<-var.get.nc(nc,"LeadTime")
qobs<-var.get.nc(nc,"qobs")
qsim<-var.get.nc(nc,"vf")

yi=which(syear==fyears)
mi=which(smonth==fmonths)
si=which(syear!=scenario)

scenario_si<-scenario[si]

qobs_sel<-qobs[yi,mi,]
qsim_sel<-qsim[yi,mi,si,]
ymax=max(c(qobs_sel,qsim_sel),na.rm=TRUE)
ymax2=max(flood_values[ci,4:9])
ymax=max(ymax,ymax2)

maxind<-which(max(qsim_sel,na.rm=TRUE)==qsim_sel,arr.ind=TRUE)[1]
qsim_max<-na.omit(qsim_sel[maxind,])
nlead<-length(qsim_max)

windows(12,8)
plot(qobs_sel,col=1,ylim=c(0,ymax),xlim=c(1,nlead),xlab="Date",type='l',ylab="Q (m3/s)",main=flood_values[ci,2],xaxt='n')
axis(1,at=xtt,labels=dlabel)

if(ptype=="spaghetti"){
for(i in 1 : dim(qsim_sel)[1]){
lines(qsim_sel[i,],col="grey")
}
}

if(ptype=="greyshade"){
q25<-na.omit(apply(qsim_sel,2,quantile,0.25,na.rm=TRUE))
q75<-na.omit(apply(qsim_sel,2,quantile,0.75,na.rm=TRUE))
ii<-length(q25)
q0025<-na.omit(apply(qsim_sel,2,quantile,0.025,na.rm=TRUE))
q0975<-na.omit(apply(qsim_sel,2,quantile,0.975,na.rm=TRUE))
ii<-length(q25)


polygon(c(1:ii,ii:1),c(q0975[1:ii],q0025[ii:1]),col="grey80")
polygon(c(1:ii,ii:1),c(q75[1:ii],q25[ii:1]),col="grey30")
}


lines(qobs_sel,col="blue")
lines(qsim_max,col="cyan")

abline(h=flood_values[ci,4],col="yellow")
abline(h=flood_values[ci,5],col="orange")
abline(h=flood_values[ci,6],col="red")

abline(h=flood_values[ci,7],col="yellow",lty=2)
abline(h=flood_values[ci,8],col="orange",lty=2)
abline(h=flood_values[ci,9],col="red",lty=2)

nover<-c(1:3)
for(j in 1 : 3){
ttest<-qsim_sel>flood_values[ci,6+j]
nover_t<-rowSums(ttest,na.rm=TRUE)
nover_t[nover_t>=1]=1
nover[j]=mean(nover_t)
text(10,flood_values[ci,6+j],round(nover[j],2))
}

if(ptype=="greyshade"){
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"50% interval","95% interval","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)",
"Mean flood (hbv)","5 years flood (hbv)", "50 years flood (hbv)"),
lty=c(1,1,1,1,1,1,1,2,2,2),col=c("blue","cyan",0,0,"yellow","orange","red","yellow","orange","red"),
fill=c(0,0,"grey30","grey80",0,0,0,0,0,0),border=c(0,0,1,1,0,0,0,0,0,0))
}

if(ptype=="spaghetti"){
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"Scenarios","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)",
"Mean flood (hbv)","5 years flood (hbv)", "50 years flood (hbv)"),
lty=c(1,1,1,1,1,1,2,2,2),col=c("blue","cyan","grey","grey","yellow","orange","red","yellow","orange","red"))

}

}





analyze_forecast<-function(rnr,hnr,smonth){

ci<-which(flood_values[,3]==paste(rnr,'_',hnr,sep=''))
ndays<-c(31,28,31,30,31,30,31)
lmax<-sum(ndays[smonth:7])

nc<-open.nc(paste("seasonal_forecast_database_",rnr,'.',hnr,".nc",sep=''),main=flood_values[ci,3])
fyears<-var.get.nc(nc,"Year")
fmonths<-var.get.nc(nc,"Month")
scenario<-var.get.nc(nc,"Scenario")
leadtime<-var.get.nc(nc,"LeadTime")
qobs<-var.get.nc(nc,"qobs")
qsim<-var.get.nc(nc,"vf")

#yi=which(syear==fyears)
mi=which(smonth==fmonths)
#si=which(syear!=scenario)
maxvalues_sim<-matrix(ncol=length(scenario)-1,nrow=length(fyears))
maxvalues_obs<-rep(NA,length(fyears))
ens.ref<-matrix(ncol=(length(scenario)-1),nrow=length(fyears))
pp_qmean<-rep(NA,length(fyears))
pp_q5<-rep(NA,length(fyears))
pp_q50<-rep(NA,length(fyears))
bb_qmean<-rep(NA,length(fyears))
bb_q5<-rep(NA,length(fyears))
bb_q50<-rep(NA,length(fyears))

Pthresholds<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
False_alarm_rate<-matrix(ncol=length(Pthresholds),nrow=3)
Hit_rate<-False_alarm_rate
csi<-False_alarm_rate

for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
qobs_sel<-qobs[i,mi,1:lmax]
qsim_sel<-qsim[i,mi,si,1:lmax]
maxvalues_sim[i,]<-apply(qsim_sel,1,max,na.rm=TRUE)
if ( any(is.na(qobs_sel))){
maxvalues_obs[i]<-NA
bb_qmean[i]<-NA
bb_q5[i]<- NA
bb_q50[i]<- NA
pp_qmean[i]<-NA
pp_q5[i]<- NA
pp_q50[i]<- NA
}
else {
maxvalues_obs[i]<-max(qobs_sel,na.rm=TRUE)
bb_qmean[i]<-as.integer(maxvalues_obs[i]>flood_values[ci,7])
bb_q5[i]<- as.integer(maxvalues_obs[i]>flood_values[ci,8])
bb_q50[i]<- as.integer(maxvalues_obs[i]>flood_values[ci,9])
pp_qmean[i]<-mean(maxvalues_sim[i,]>flood_values[ci,7])
pp_q5[i]<- mean(maxvalues_sim[i,]>flood_values[ci,8])
pp_q50[i]<- mean(maxvalues_sim[i,]>flood_values[ci,9])
}

}



for(j in 1 : 3){
outcome<-as.integer(maxvalues_obs>flood_values[ci,(3+j)])
for(k in 1 : length(Pthresholds)){
hits<-sum((apply(maxvalues_sim,1,quantile,Pthresholds[k])>flood_values[ci,(6+j)])& outcome==1,na.rm=TRUE) 
false<-sum((apply(maxvalues_sim,1,quantile,Pthresholds[k])>flood_values[ci,(6+j)])& outcome==0,na.rm=TRUE) 
misses<-sum((apply(maxvalues_sim,1,quantile,Pthresholds[k])<flood_values[ci,(6+j)])& outcome==1,na.rm=TRUE)
noevents<-sum((apply(maxvalues_sim,1,quantile,Pthresholds[k])<flood_values[ci,(6+j)])& outcome==0,na.rm=TRUE)
Hit_rate[j,k]<-hits/(hits+misses)
False_alarm_rate[j,k]<-false/(false+noevents)
csi[j,k]<-hits/(false+hits+misses)
if((false+noevents)==0)False_alarm_rate[j,k]<-0
if((hits+misses)==0)Hit_rate[j,k]<-0
if((false+hits+misses)==0)csi[j,k]<-0
}
}

windows(12,12)
par(mfrow=c(2,2))
tt<-c("Mean flood","Q5","Q50")
for(j in 1: 3){
plot(False_alarm_rate[j,],Hit_rate[j,],xlim=c(0,1),ylim=c(0,1),xlab="False alarm rate",ylab="Hit rate",main=tt[j])
for(l in 1 : length(Pthresholds)) text(False_alarm_rate[j,l],(Hit_rate[j,l])+0.03,Pthresholds[l],cex=0.7)
abline(0,1)
}

windows(12,12)
par(mfrow=c(2,2))
tt<-c("Mean flood","Q5","Q50")
for(j in 1: 3){
plot(Pthresholds,csi[j,],xlim=c(0,1),ylim=c(0,1),ylab="CSI",xlab="Probability level",main=tt[j])
#for(l in 1 : length(Pthresholds)) text(False_alarm_rate[j,l],(Hit_rate[j,l])+0.03,Pthresholds[l],cex=0.7)
#abline(0,1)
}




for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
ens.ref[i,]<-maxvalues_obs[si]
}

windows(12,8)
boxplot(t(maxvalues_sim),xaxt="n",ylab="Flood (m3/s)",xlab="Year")
axis(1,at=c(1:length(fyears)),labels=as.character(fyears))
points(maxvalues_obs,col="blue",pch=15)

qsim_95<-apply(maxvalues_sim,1,quantile,0.95)
qsim_5<-apply(maxvalues_sim,1,quantile,0.05)
qsim_50<-apply(maxvalues_sim,1,quantile,0.5)
qsim_max<-apply(maxvalues_sim,1,max,na.rm=TRUE)

ssmax<-max(qsim_max,maxvalues_obs,na.rm=TRUE)
ssmin<-min(qsim_5,maxvalues_obs,na.rm=TRUE)
windows(8,8)
plot(maxvalues_obs,qsim_50,xlim=c(ssmin,ssmax),ylim=c(ssmin,ssmax),xlab="Observed floods (m3/s)",ylab="Scenario (m3/s)")
points(maxvalues_obs,qsim_max,col="red")
points(maxvalues_obs,qsim_95,col="orange")
points(maxvalues_obs,qsim_5,col="green")
legend('bottomright',legend=c(paste("Median scenario cor:",round(cor(maxvalues_obs,qsim_50,use="pairwise.complete.obs"),2)),
paste("Maximum scenario cor:",round(cor(maxvalues_obs,qsim_max,use="pairwise.complete.obs"),2))
,paste("95% scenariocor:",round(cor(maxvalues_obs,qsim_95,use="pairwise.complete.obs"),2))
,paste("5% scenario cor:",round(cor(maxvalues_obs,qsim_5,use="pairwise.complete.obs"),2)))
,col=c(1,"red","orange","green"),pch=c(1,1,1,1))
abline(0,1)

windows(6,4)
ReliabilityDiagram(na.omit(pp_qmean), na.omit(bb_qmean),plot=TRUE)
windows(6,4)
ReliabilityDiagram(na.omit(pp_q5), na.omit(bb_q5),plot=TRUE)
windows(6,4)
ReliabilityDiagram(na.omit(pp_q50), na.omit(bb_q50),plot=TRUE)

windows(6,4)
out<-list()
out$rankhist<-Rankhist(maxvalues_sim[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)])
PlotRankhist(out$rankhist)
out$sim_max<-maxvalues_sim
out$obs_max<-maxvalues_obs
out$crps<-EnsCrps(maxvalues_sim[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)])
out$brier_qmean<-EnsBrier(maxvalues_sim[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,7])
out$brier_q5<-EnsBrier(maxvalues_sim[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,8])
out$brier_q50<-EnsBrier(maxvalues_sim[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,9])

qsim_50<-apply(maxvalues_sim,1,quantile,0.5)
out$rmse<-sqrt(mean((qsim_50[!is.na(maxvalues_obs)]-maxvalues_obs[!is.na(maxvalues_obs)])^2))
out$Reff<-1.0-(mean((qsim_50[!is.na(maxvalues_obs)]-maxvalues_obs[!is.na(maxvalues_obs)])^2)/var(maxvalues_obs[!is.na(maxvalues_obs)]))
out$corr<-cor(maxvalues_obs[!is.na(maxvalues_obs)],qsim_50[!is.na(maxvalues_obs)])
out$crpss<-EnsCrpss(maxvalues_sim[!is.na(maxvalues_obs),], ens.ref[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)])out$brierss_qmean<-EnsBrierSs(maxvalues_sim[!is.na(maxvalues_obs),], ens.ref[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,7])
out$brierss_q5<-EnsBrierSs(maxvalues_sim[!is.na(maxvalues_obs),], ens.ref[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,8])
out$brierss_q50<-EnsBrierSs(maxvalues_sim[!is.na(maxvalues_obs),], ens.ref[!is.na(maxvalues_obs),], maxvalues_obs[!is.na(maxvalues_obs)],flood_values[ci,9])
out$csi<-csi
out$ndata<-30
out
}


