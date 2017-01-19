
plot_forecast<-function(rnr,hnr,syear,smonth,flood_values,ptype="greyshade", fpath="../data/netcdf/",
qtrans=NA,empqt=FALSE){

if (!require('RNetCDF')) {
    stop('The package RNetCDF was not installed')
  }

  if (!require('nsRFA')) {
    stop('The package nsRFA was not installed')
  }

ndays<-c(1,31,28,31,30,31,30,31,31,30,31,30)
fd<-c("01/01/","01/02/","01/03/","01/04/","01/05/","01/06/","01/07/","01/08/","01/09/","01/10/","01/11/","01/12/")
ttext<-c(paste(fd,syear,sep=""))
ci<-which(flood_values[,3]==paste(rnr,'_',hnr,sep=''))
xtt<-cumsum(ndays)[smonth:12]
xtt<-xtt-xtt[1]+1
dlabel<-ttext[smonth:12]

ldays<-c(31,28,31,30,31,30,31)
lmax<-sum(ldays[smonth:7])


nc<-RNetCDF::open.nc(paste(fpath,"seasonal_forecast_database_",rnr,'.',hnr,".nc",sep=''),main=flood_values[ci,3])
fyears<-RNetCDF::var.get.nc(nc,"Year")
fmonths<-RNetCDF::var.get.nc(nc,"Month")
scenario<-RNetCDF::var.get.nc(nc,"Scenario")
leadtime<-RNetCDF::var.get.nc(nc,"LeadTime")
qobs<-RNetCDF::var.get.nc(nc,"qobs")
qsim<-RNetCDF::var.get.nc(nc,"vf")

yi=which(syear==fyears)
mi=which(smonth==fmonths)
si=which(syear!=scenario)

scenario_si<-scenario[si]

qobs_sel<-qobs[yi,mi,]
qsim_sel<-qsim[yi,mi,si,]

if(empqt){
  empv<-read.table(file=paste("inst/empq/",rnr,'.',hnr,".txt",sep=''),sep=";",header=TRUE)  
  myfun<-approxfun(sort(empv[,1]),sort(empv[,2]))
 qsim_sel_t<-qsim_sel 
  simdim<-dim(qsim_sel_t)
  for(m in 1:simdim[1]){
    for(n in 1:simdim[2]){
        qsim_sel_t[m,n]<-myfun(qsim_sel[m,n])
    }
  }
 qsim_sel<-qsim_sel_t 
}

if(any(!is.na(qtrans))){
simdim<-dim(qsim_sel)
qsim_sel_t<-qsim_sel
tnumbers<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)),ncol=5,byrow=TRUE)
i2<-which(as.integer(tnumbers[,1])==rnr & as.integer(tnumbers[,2])==hnr)
parsel<-as.numeric(qtrans[i2,3:6])
tdist=qtrans[,7]
for(m in 1:simdim[1]){
for(n in 1:simdim[2]){
if(tdist=='gumbel'){
FF1<-nsRFA::F.gumb(qsim_sel[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-nsRFA::invF.gumb (FF1,parsel[1], parsel[2])
}
if(tdist=='gamma'){
FF1<-pgamma(qsim_sel[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-qgamma(FF1,parsel[1], parsel[2])
}
}
}

qsim_sel<-qsim_sel_t
}


ymax=max(c(qobs_sel,qsim_sel),na.rm=TRUE)
ymax2=max(flood_values[ci,4:9])
ymax=max(ymax,ymax2)

maxind<-which(max(qsim_sel[,1:lmax],na.rm=TRUE)==qsim_sel,arr.ind=TRUE)[1]
qsim_max<-na.omit(qsim_sel[maxind,])
nlead<-length(qsim_max)

windows(12,8)
plot(qobs_sel,col=1,ylim=c(0,ymax),xlim=c(1,nlead),xlab="Date",type='l',ylab="Q (m3/s)",main=flood_values[ci,2],xaxt='n')
axis(1,at=xtt,labels=dlabel)

if(ptype=="spaghetti"){
for(i in 1 : dim(qsim_sel)[1]){
lines(qsim_sel_t[i,],col="grey")
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

if(is.na(qtrans)&!empqt){
#print(qtrans)
#print(empqt)
abline(h=flood_values[ci,7],col="yellow",lty=2)
abline(h=flood_values[ci,8],col="orange",lty=2)
abline(h=flood_values[ci,9],col="red",lty=2)
}


nover<-c(1:3)
for(j in 1 : 3){
if(any(is.na(qtrans))&!empqt){
ttest<-qsim_sel>flood_values[ci,6+j]
}
else{
ttest<-qsim_sel>flood_values[ci,3+j]
}

nover_t<-rowSums(ttest,na.rm=TRUE)
nover_t[nover_t>=1]=1
nover[j]=mean(nover_t)

if(is.na(qtrans)&!empqt){
text(10,flood_values[ci,6+j],round(nover[j],2))
}
else{
text(10,flood_values[ci,3+j],round(nover[j],2))
}


}

if(ptype=="greyshade"){
if(is.na(qtrans)&!empqt){
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"50% interval","95% interval","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)",
"Mean flood (hbv)","5 years flood (hbv)", "50 years flood (hbv)"),
lty=c(1,1,1,1,1,1,1,2,2,2),col=c("blue","cyan",0,0,"yellow","orange","red","yellow","orange","red"),
fill=c(0,0,"grey30","grey80",0,0,0,0,0,0),border=c(0,0,1,1,0,0,0,0,0,0))
}else{
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"50% interval","95% interval","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)"),
lty=c(1,1,1,1,1,1,1),col=c("blue","cyan",0,0,"yellow","orange","red"),
fill=c(0,0,"grey30","grey80",0,0,0),border=c(0,0,1,1,0,0,0))

}

}

if(ptype=="spaghetti"){
if(is.na(qtrans)&!empqt){
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"Scenarios","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)",
"Mean flood (hbv)","5 years flood (hbv)", "50 years flood (hbv)"),
lty=c(1,1,1,1,1,1,2,2,2),col=c("blue","cyan","grey","grey","yellow","orange","red","yellow","orange","red"))

}
}else{
legend('topright',legend=c("Observed streamflow",paste("Maximum scenario:",scenario_si[maxind]),"50% interval","95% interval","Mean flood (obs)","5 years flood (obs)", "50 years flood (obs)"),
lty=c(1,1,1,1,1,1,1),col=c("blue","cyan",0,0,"yellow","orange","red"),
fill=c(0,0,"grey30","grey80",0,0,0),border=c(0,0,1,1,0,0,0))

}


}





analyze_forecast<-function(rnr,hnr,smonth,fpath="../data/netcdf/",flood_values,qtrans=NA,mplot=TRUE,empqt=FALSE){
if (!require('RNetCDF')) {
    stop('The package RNetCDF was not installed')
  }

  if (!require('nsRFA')) {
    stop('The package nsRFA was not installed')
  }

   if (!require('verification')) {
    stop('The package verification was not installed')
  }

  if (!require('SpecsVerification')) {
    stop('The package SpecsVerification was not installed')
  }  
  
  
ci<-which(flood_values[,3]==paste(rnr,'_',hnr,sep=''))
ndays<-c(31,28,31,30,31,30,31)
lmax<-sum(ndays[smonth:7])

nc<-RNetCDF::open.nc(paste(fpath,"seasonal_forecast_database_",rnr,'.',hnr,".nc",sep=''),main=flood_values[ci,3])
fyears<-RNetCDF::var.get.nc(nc,"Year")
fmonths<-RNetCDF::var.get.nc(nc,"Month")
scenario<-RNetCDF::var.get.nc(nc,"Scenario")
leadtime<-RNetCDF::var.get.nc(nc,"LeadTime")
qobs<-RNetCDF::var.get.nc(nc,"qobs")
qsim<-RNetCDF::var.get.nc(nc,"vf")

#yi=which(syear==fyears)
mi=which(smonth==fmonths)
mi2=which(1==fmonths)
#si=which(syear!=scenario)
maxvalues_sim<-matrix(ncol=length(scenario)-1,nrow=length(fyears))
maxvalues_obs<-rep(NA,length(fyears))
maxvalues_obs2<-rep(NA,length(fyears))
maxindex_obs<-rep(NA,length(fyears))
volumes_sim<-maxvalues_sim
volumes_obs<-maxvalues_obs
ens.ref<-matrix(ncol=(length(scenario)-1),nrow=length(fyears))
pp_qmean<-rep(NA,length(fyears))
pp_q5<-rep(NA,length(fyears))
pp_q50<-rep(NA,length(fyears))
bb_qmean<-rep(NA,length(fyears))
bb_q5<-rep(NA,length(fyears))
bb_q50<-rep(NA,length(fyears))
bb_qsimmean<-maxvalues_sim
bb_qsimq5<-maxvalues_sim
bb_qsimq50<-maxvalues_sim

bb_v_q50<- maxvalues_obs
bb_v_q20<- maxvalues_obs

bb_v_qsimq50<- maxvalues_sim
bb_v_qsimq20<- maxvalues_sim


#Pthresholds<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
#False_alarm_rate<-matrix(ncol=length(Pthresholds),nrow=3)
#Hit_rate<-False_alarm_rate
#csi<-False_alarm_rate

for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
qobs_sel<-qobs[i,mi,1:lmax]
if ( any(is.na(qobs_sel)))volumes_obs[i]<-NA
else volumes_obs[i]<-sum(qobs_sel,na.rm=TRUE)*60*60*24/1000000
}
volume_q20<-quantile(volumes_obs,0.2,na.rm=TRUE)
volume_q50<-quantile(volumes_obs,0.5,na.rm=TRUE)


for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
qobs_sel<-qobs[i,mi,1:lmax]
qobs_sel2<-qobs[i,mi2,1:lmax]
qsim_sel<-qsim[i,mi,si,1:lmax]
maxvalues_sim[i,]<-apply(qsim_sel,1,max,na.rm=TRUE)
volumes_sim[i,]<-apply(qsim_sel,1,sum,na.rm=TRUE)*60*60*24/1000000
if ( any(is.na(qobs_sel))){
maxvalues_obs[i]<-NA
maxindex_obs[i]<-NA
bb_qmean[i]<-NA
bb_q5[i]<- NA
bb_q50[i]<- NA
pp_qmean[i]<-NA
pp_q5[i]<- NA
pp_q50[i]<- NA
bb_qsimmean[i,]<-NA
bb_qsimq5[i,]<-NA
bb_qsimq50[i,]<-NA
}
else {
maxvalues_obs[i]<-max(qobs_sel,na.rm=TRUE)
maxvalues_obs2[i]<-max(qobs_sel2,na.rm=TRUE)
maxindex_obs[i]<-which.max(qobs_sel)
bb_qmean[i]<-as.integer(maxvalues_obs[i]>flood_values[ci,4])
bb_q5[i]<- as.integer(maxvalues_obs[i]>flood_values[ci,5])
bb_q50[i]<- as.integer(maxvalues_obs[i]>flood_values[ci,6])
pp_qmean[i]<-mean(maxvalues_sim[i,]>flood_values[ci,7])
pp_q5[i]<- mean(maxvalues_sim[i,]>flood_values[ci,8])
pp_q50[i]<- mean(maxvalues_sim[i,]>flood_values[ci,9])
bb_qsimmean[i,]<-maxvalues_sim[i,]>flood_values[ci,7]
bb_qsimq5[i,]<- maxvalues_sim[i,]>flood_values[ci,8]
bb_qsimq50[i,]<- maxvalues_sim[i,]>flood_values[ci,9]



bb_v_q50[i]<- as.integer(volumes_obs[i]<volume_q50)
bb_v_q20[i]<- as.integer(volumes_obs[i]<volume_q20)

bb_v_qsimq50[i,]<- as.integer(volumes_sim[i,]<volume_q50)
bb_v_qsimq20[i,]<- as.integer(volumes_sim[i,]<volume_q20)



}

}

ens.ref<-matrix(ncol=(length(na.omit(maxvalues_obs))-1),nrow=length(fyears))
ens.ref.v<-ens.ref
bb_ensmean<-ens.ref
bb_ensq5<-ens.ref
bb_ensq50<-ens.ref

bb_v_ensq50<-ens.ref
bb_v_ensq20<-ens.ref

bb_v50<-ens.ref
bb_v20<-ens.ref

#print(dim(bb_ensmean))
for(i in 1 : length(fyears)){
if(is.na(bb_qmean[i])){
ens.ref[i,]<-NA
bb_ensmean[i,]<-NA
bb_ensq5[i,]<-NA
bb_ensq50[i,]<-NA

ens.ref.v[i,]<-NA
bb_v50[i,]<-NA
bb_v20[i,]<-NA

}
else{
si=which(fyears[i]!=scenario)
ens.ref[i,]<-as.numeric(na.omit(maxvalues_obs[si]))
bb_ensmean[i,]<-ens.ref[i,]>flood_values[ci,4]
bb_ensq5[i,]<- ens.ref[i,]>flood_values[ci,5]
bb_ensq50[i,]<- ens.ref[i,]>flood_values[ci,6]

ens.ref.v[i,]<-as.numeric(na.omit(volumes_obs[si]))
bb_v_ensq50[i,]<-ens.ref.v[i,]<volume_q50
bb_v_ensq20[i,]<- ens.ref.v[i,]<volume_q20



}
}

#print(dim(bb_ensmean))

bb_qsimmean<-bb_qsimmean[!is.na(maxvalues_obs),]
bb_qsimq5<-bb_qsimq5[!is.na(maxvalues_obs),]
bb_qsimq50<-bb_qsimq50[!is.na(maxvalues_obs),]

bb_ensmean<-bb_ensmean[!is.na(maxvalues_obs),]
bb_ensq5<-bb_ensq5[!is.na(maxvalues_obs),]
bb_ensq50<-bb_ensq50[!is.na(maxvalues_obs),]


bb_ensmean<-bb_ensmean[,!is.na(bb_ensmean[1,])]
bb_ensq5<-bb_ensq5[,!is.na(bb_ensmean[1,])]
bb_ensq50<-bb_ensq50[,!is.na(bb_ensmean[1,])]




bb_qmean<-bb_qmean[!is.na(maxvalues_obs)]
bb_q5<-bb_q5[!is.na(maxvalues_obs)]
bb_q50<-bb_q50[!is.na(maxvalues_obs)]
pp_qmean<-pp_qmean[!is.na(maxvalues_obs)]
pp_q5<-pp_q5[!is.na(maxvalues_obs)]
pp_q50<-pp_q50[!is.na(maxvalues_obs)]


maxvalues_sim_o<-maxvalues_sim
if(any(!is.na(qtrans))){
simdim<-dim(maxvalues_sim)
qsim_sel_t<-maxvalues_sim
tnumbers<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)),ncol=5,byrow=TRUE)
i2<-which(as.integer(tnumbers[,1])==rnr & as.integer(tnumbers[,2])==hnr)
parsel<-as.numeric(qtrans[i2,3:6])
tdist=qtrans[i2,7]
for(m in 1:simdim[1]){
for(n in 1:simdim[2]){
if(tdist=='gumbel'){
FF1<-nsRFA::F.gumb(maxvalues_sim[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-nsRFA::invF.gumb(FF1,parsel[1], parsel[2])
}
if(tdist=='gamma'){
FF1<-pgamma(maxvalues_sim[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-qgamma(FF1,parsel[1], parsel[2])
}
}
}
maxvalues_sim<-qsim_sel_t
}

if(empqt){
#print("empqt")
 empv<-read.table(file=paste("inst/empq/",rnr,'.',hnr,".txt",sep=''),sep=";",header=TRUE)  
 myfun<-approxfun(sort(empv[,1]),sort(empv[,2]))
 simdim<-dim(maxvalues_sim)
qsim_sel_t<-maxvalues_sim

for(m in 1:simdim[1]){
for(n in 1:simdim[2]){
        qsim_sel_t[m,n]<-myfun(maxvalues_sim[m,n])
#		print(c(m,n,qsim_sel_t[m,n]))
    }
  }
 maxvalues_sim<-qsim_sel_t 
}





#maxvalues_obs<-maxvalues_obs[maxvalues_sim[,1]>0]
#ens.ref<-ens.ref[maxvalues_sim[,1]>0,]
#maxvalues_sim<-maxvalues_sim[maxvalues_sim[,1]>0,]


roc_area<-c()
roc_pvalue<-c()

roc_area_v<-c()
roc_pvalue_v<-c()
csi<-c()
csi_v<-c()

outcome<-as.integer(volumes_obs<volume_q20)
forecasted<-as.integer(apply(volumes_sim,1,quantile,0.5)<volume_q20)
forecasted<-forecasted[!is.na(outcome)]
p_volume<-rowMeans(volumes_sim<volume_q20)
p_volume<-p_volume[!is.na(outcome)]
outcome<-na.omit(outcome)
stemp=(sum(outcome) + sum(forecasted&!outcome))
csi_v[1]<- (-1.0)
if(stemp>0) csi_v[1]=sum(outcome&forecasted) / stemp
if(mplot){
windows(6,6)
roc_out<-verification::roc.plot(outcome,p_volume,thresholds=seq(0.0,1.0,by=0.1))
roc_area_v[1]=roc_out$roc.vol$Area
roc_pvalue_v[1]=roc_out$roc.vol$p.value
}
else {
roc_out<-verification::roc.area(outcome,p_volume)
roc_area_v[1]=roc_out$A
roc_pvalue_v[1]=roc_out$p.value
}

outcome<-as.integer(volumes_obs<volume_q50)
forecasted<-as.integer(apply(volumes_sim,1,quantile,0.5)<volume_q50)
forecasted<-forecasted[!is.na(outcome)]
p_volume<-rowMeans(volumes_sim<volume_q50)
p_volume<-p_volume[!is.na(outcome)]
outcome<-na.omit(outcome)
stemp=(sum(outcome) + sum(forecasted&!outcome))
csi_v[2]<- (-1.0)
if(stemp>0) csi_v[2]=sum(outcome&forecasted) / stemp
if(mplot){
windows(6,6)
roc_out<-verification::roc.plot(outcome,p_volume,thresholds=seq(0.0,1.0,by=0.1))
roc_area_v[1]=roc_out$roc.vol$Area
roc_pvalue_v[1]=roc_out$roc.vol$p.value
}
else {
roc_out<-verification::roc.area(outcome,p_volume)
roc_area_v[1]=roc_out$A
roc_pvalue_v[1]=roc_out$p.value
}





for(j in 1 : 3){
#print(dim(maxvalues_sim))
#print(length(maxvalues_obs))
#print(apply(maxvalues_sim,1,quantile,0.5))
#Use untransformed ensembles and compare to HBV model flood quantiles
outcome<-as.integer(maxvalues_obs>flood_values[ci,(3+j)])
forecasted<-as.integer(apply(maxvalues_sim_o,1,quantile,0.5)>flood_values[ci,(6+j)])
#if(sum(na.omit(outcome))==0)
#roc_area[j] <- 0
#roc_pvalue[j] <- 0
#csi[j]<- -1.0
p_flood<-rowMeans(maxvalues_sim_o>flood_values[ci,(6+j)])
p_flood<-p_flood[!is.na(outcome)]
forecasted<-forecasted[!is.na(outcome)]
outcome<-na.omit(outcome)
#print(outcome)
#print(forecasted)
stemp=(sum(outcome) + sum(forecasted&!outcome))
csi[j]<- (-1.0)
if(stemp>0) csi[j]=sum(outcome&forecasted) / stemp

if(sum(na.omit(outcome))==0){
roc_area[j] <- -0.5
roc_pvalue[j] <- 1.0
}
else{
if(mplot){
windows(6,6)
roc_out<-verification::roc.plot(outcome,p_flood,thresholds=seq(0.0,1.0,by=0.1))
roc_area[j]=roc_out$roc.vol$Area
roc_pvalue[j]=roc_out$roc.vol$p.value
}
else {
roc_out<-verification::roc.area(outcome,p_flood)
roc_area[j]=roc_out$A
roc_pvalue[j]=roc_out$p.value
}
}
}




if(mplot){
windows(12,8)
boxplot(t(maxvalues_sim),xaxt="n",ylab="Flood (m3/s)",xlab="Year")
axis(1,at=c(1:length(fyears)),labels=as.character(fyears))
points(maxvalues_obs,col="blue",pch=15)

windows(12,8)
boxplot(t(volumes_sim),xaxt="n",ylab="Volumes (M m3)",xlab="Year")
axis(1,at=c(1:length(fyears)),labels=as.character(fyears))
points(volumes_obs,col="blue",pch=15)

}




qsim_95<-apply(maxvalues_sim,1,quantile,0.95,na.rm=TRUE)
qsim_5<-apply(maxvalues_sim,1,quantile,0.05,na.rm=TRUE)
qsim_50<-apply(maxvalues_sim,1,quantile,0.5,na.rm=TRUE)
qsim_max<-apply(maxvalues_sim,1,max,na.rm=TRUE)

ssmax<-max(qsim_max,maxvalues_obs,na.rm=TRUE)
ssmin<-min(qsim_5,maxvalues_obs,na.rm=TRUE)


qsim_v_95<-apply(volumes_sim,1,quantile,0.95,na.rm=TRUE)
qsim_v_5<-apply(volumes_sim,1,quantile,0.05,na.rm=TRUE)
qsim_v_50<-apply(volumes_sim,1,quantile,0.5,na.rm=TRUE)
qsim_v_max<-apply(volumes_sim,1,max,na.rm=TRUE)

ssmax_v<-max(qsim_v_max,volumes_obs,na.rm=TRUE)
ssmin_v<-min(qsim_v_5,volumes_obs,na.rm=TRUE)



if(mplot){
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


windows(8,8)
plot(volumes_obs,qsim_v_50,xlim=c(ssmin_v,ssmax_v),ylim=c(ssmin_v,ssmax_v),xlab="Observed volume (M m3)",ylab="Scenario (M m3)")
points(volumes_obs,qsim_v_max,col="red")
points(volumes_obs,qsim_v_95,col="orange")
points(volumes_obs,qsim_v_5,col="green")
legend('bottomright',legend=c(paste("Median scenario cor:",round(cor(volumes_obs,qsim_v_50,use="pairwise.complete.obs"),2)),
paste("Maximum scenario cor:",round(cor(volumes_obs,qsim_v_max,use="pairwise.complete.obs"),2))
,paste("95% scenariocor:",round(cor(volumes_obs,qsim_v_95,use="pairwise.complete.obs"),2))
,paste("5% scenario cor:",round(cor(volumes_obs,qsim_v_5,use="pairwise.complete.obs"),2)))
,col=c(1,"red","orange","green"),pch=c(1,1,1,1))
abline(0,1)



windows(6,4)
SpecsVerification::ReliabilityDiagram(pp_qmean, bb_qmean,plot=TRUE)
windows(6,4)
SpecsVerification::ReliabilityDiagram(pp_q5, bb_q5,plot=TRUE)
windows(6,4)
SpecsVerification::ReliabilityDiagram(pp_q50, bb_q50,plot=TRUE)
}
out<-list()
out$roc_area<-roc_area
out$roc_pvalue<-roc_pvalue

out$roc_area_v<-roc_area_v
out$roc_pvalue_v<-roc_pvalue_v

out$rankhist<-SpecsVerification::Rankhist(maxvalues_sim, maxvalues_obs)
if(mplot){
windows(6,4)
SpecsVerification::PlotRankhist(out$rankhist)
}
out$sim_max<-maxvalues_sim
out$obs_max<-maxvalues_obs
out$obs_max_ind<-maxindex_obs
out$crps<-SpecsVerification::EnsCrps(maxvalues_sim, maxvalues_obs)
out$crps_v<-SpecsVerification::EnsCrps(volumes_sim, volumes_obs) 
out$brier_qmean<-SpecsVerification::EnsBrier(bb_qsimmean,bb_qmean,0.5)
out$brier_q5<-SpecsVerification::EnsBrier(bb_qsimq5,bb_q5,0.5)
out$brier_q50<-SpecsVerification::EnsBrier(bb_qsimq50,bb_q50,0.5)
qsim_50<-apply(maxvalues_sim,1,quantile,0.5,na.rm=TRUE)
qsim_v_50<-apply(volumes_sim,1,quantile,0.5,na.rm=TRUE)
out$rmse<-sqrt(mean((qsim_50-maxvalues_obs)^2))
out$Reff<-1.0-(mean((qsim_50-maxvalues_obs)^2,na.rm=TRUE)/var(maxvalues_obs,na.rm=TRUE))
out$Reff_v<-1.0-(mean((qsim_v_50-volumes_obs)^2,na.rm=TRUE)/var(volumes_obs,na.rm=TRUE))
out$corr<-cor(maxvalues_obs,qsim_50,use="pairwise.complete.obs")
out$corr_v<-cor(volumes_obs,qsim_v_50,use="pairwise.complete.obs")
out$crpss<-SpecsVerification::EnsCrpss(maxvalues_sim, ens.ref, maxvalues_obs) 
#print('test')
#print(dim(bb_ensmean))
out$brierss_qmean<-SpecsVerification::EnsBrierSs(bb_qsimmean, bb_ensmean, bb_qmean,0.5)
#print('test')
#print(dim(bb_ensq5))
out$brierss_q5<-SpecsVerification::EnsBrierSs(bb_qsimq5, bb_ensq5, bb_q5,0.5)
#print('test')
#print(dim(bb_ensq50))
out$brierss_q50<-SpecsVerification::EnsBrierSs(bb_qsimq50, bb_ensq50, bb_q50,0.5)

out$crpss_v<-SpecsVerification::EnsCrpss(volumes_sim, ens.ref.v, volumes_obs) 

out$brierss_q50<-SpecsVerification::EnsBrierSs(bb_v_qsimq50,bb_v_ensq50, bb_v_q50,0.5)

out$brierss_q20<-SpecsVerification::EnsBrierSs(bb_v_qsimq20,bb_v_ensq20, bb_v_q20,0.5)


out$csi<-csi
out$csi_v<-csi_v
out
}


analyse_all<-function(npath,flood_values,qtrans,empqt=FALSE){
	out_analyze<-list()
	temp<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)) ,ncol=5,byrow=TRUE)[,1:2]
	for(i in 1 : dim(qtrans)[1]){
		out_analyze[[i]]<-list() 
		rnr<-as.integer(temp[i,1])
		hnr<-as.integer(temp[i,2])	
 #       print(temp[i,])		
		for( j in 1 : 7){
#		print(c(i,j))
            if(empqt) out_analyze[[i]][[j]]<-analyze_forecast(rnr,hnr,j,fpath=npath,flood_values=flood_values,mplot=FALSE,empqt=TRUE)
			else out_analyze[[i]][[j]]<-analyze_forecast(rnr,hnr,j,fpath=npath,flood_values=flood_values,qtrans=qtrans,mplot=FALSE)
		}
	}
	names(out_analyze)<-paste(temp[,1],'.',temp[,2],sep='')
    out_analyze$ndata<-qtrans[,2]
	return(out_analyze)
}


get_mindex_month<-function(out,mm){
ncc<-(length(out)-1)
mindex<-rep(NA,ncc)
for(i in 1 : ncc)
mindex[i]<- median(out[[i]][[mm]]$obs_max_ind,na.rm=TRUE)
names(mindex)<-names(out)[1:ncc]
return(mindex)
}



get_mindex_all<-function(out,ndata){
ncc<-(length(out)-1)
mindex<-matrix(NA,ncol=7,nrow=ncc)
colnames(mindex)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul')
rownames(mindex)<-names(out)[1:ncc]
for(i in 1 : 7){
mindex[,i]<-get_mindex_month(out,i)
}
return(mindex)
}


get_obs_max_month<-function(out,mm){
ncc<-(length(out)-1)
obs_max<-rep(NA,ncc)
for(i in 1 : ncc)
obs_max[i]<- out[[i]][[mm]]$obs_max
names(obs_max)<-names(out)[1:ncc]
return(obs_max)
}



get_obs_max_all<-function(out,ndata){
ncc<-(length(out)-1)
obs_max<-matrix(NA,ncol=7,nrow=ncc)
colnames(obs_max)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul')
rownames(obs_max)<-names(out)[1:ncc]
for(i in 1 : 7){
obs_max[,i]<-get_obs_max_month(out,i)
}
return(obs_max)
}




get_corr_month<-function(out,mm){
  ncc<-(length(out)-1)
  corr_all<-rep(NA,ncc)
  ndata<-out[[ncc+1]]
  for(i in 1 : ncc)
    corr_all[i]<- out[[i]][[mm]]$corr
  # Beregner t-value:
  tt<-abs(qt(0.05,ndata-2))
  #beregner kritisk verdi for korrelasjon
  rr= tt/sqrt(ndata-2+tt^2)
  cout<-cbind(corr_all,rr)
  rownames(cout)<-names(out)[1:ncc]
  return(cout)
}



get_corr_all<-function(out){
  ncc<-(length(out)-1)
  corr_out<-matrix(NA,ncol=14,nrow=ncc)
  colnames(corr_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_corr','Month','R_sig','Max_corr-j_m','month1_5','max_corr_j_a','month1_4')
  rownames(corr_out)<-names(out)[1:ncc]
  for(i in 1 : 7){
    corr_out[,i]<-get_corr_month(out,i)[,1]
  }
  corr_out[,8]<-apply(corr_out[,1:7],1,max)
  corr_out[,9]<-unlist(apply(corr_out[,1:7],1,which.max))
  corr_out[,10]<-get_corr_month(out,1)[,2]
  corr_out[,11]<-apply(corr_out[,1:5],1,max)
  corr_out[,12]<-unlist(apply(corr_out[,1:5],1,which.max))
  corr_out[,13]<-apply(corr_out[,1:4],1,max)
  corr_out[,14]<-unlist(apply(corr_out[,1:4],1,which.max))
  
  return(corr_out)
}





get_corr_v_month<-function(out,mm){
ncc<-(length(out)-1)
corr_all<-rep(NA,ncc)
ndata<-out[[ncc+1]]
for(i in 1 : ncc)
corr_all[i]<- out[[i]][[mm]]$corr_v
# Beregner t-value:
tt<-abs(qt(0.05,ndata-2))
#beregner kritisk verdi for korrelasjon
rr= tt/sqrt(ndata-2+tt^2)
cout<-cbind(corr_all,rr)
rownames(cout)<-names(out)[1:ncc]
return(cout)
}



get_corr_v_all<-function(out){
  ncc<-(length(out)-1)
  corr_out<-matrix(NA,ncol=14,nrow=ncc)
  colnames(corr_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_corr','Month','R_sig','Max_corr-j_m','month1_5','max_corr_j_a','month1_4')
  rownames(corr_out)<-names(out)[1:ncc]
  for(i in 1 : 7){
    corr_out[,i]<-get_corr_v_month(out,i)[,1]
  }
  corr_out[,8]<-apply(corr_out[,1:7],1,max)
  corr_out[,9]<-unlist(apply(corr_out[,1:7],1,which.max))
  corr_out[,10]<-get_corr_month(out,1)[,2]
  corr_out[,11]<-apply(corr_out[,1:5],1,max)
  corr_out[,12]<-unlist(apply(corr_out[,1:5],1,which.max))
  corr_out[,13]<-apply(corr_out[,1:4],1,max)
  corr_out[,14]<-unlist(apply(corr_out[,1:4],1,which.max))
  return(corr_out)
}







get_Reff_month<-function(out,mm){
ncc<-(length(out)-1)
Reff_all<-rep(NA,ncc)
for(i in 1 : ncc)
Reff_all[i]<- out[[i]][[mm]]$Reff
names(Reff_all)<-names(out)[1:ncc]
return(Reff_all)
}



get_Reff_all<-function(out,ndata){
ncc<-(length(out)-1)
Reff_out<-matrix(NA,ncol=9,nrow=ncc)
colnames(Reff_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_Reff','Month')
rownames(Reff_out)<-names(out)[1:ncc]
for(i in 1 : 7){
Reff_out[,i]<-get_Reff_month(out,i)
}
Reff_out[,8]<-apply(Reff_out[,1:7],1,max)
Reff_out[,9]<-apply(Reff_out[,1:7],1,which.max)

return(Reff_out)
}

get_Reff_v_month<-function(out,mm){
  ncc<-(length(out)-1)
  Reff_all<-rep(NA,ncc)
  for(i in 1 : ncc)
    Reff_all[i]<- out[[i]][[mm]]$Reff_v
  names(Reff_all)<-names(out)[1:ncc]
  return(Reff_all)
}



get_Reff_v_all<-function(out,ndata){
  ncc<-(length(out)-1)
  Reff_out<-matrix(NA,ncol=9,nrow=ncc)
  colnames(Reff_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_Reff','Month')
  rownames(Reff_out)<-names(out)[1:ncc]
  for(i in 1 : 7){
    Reff_out[,i]<-get_Reff_v_month(out,i)
  }
  Reff_out[,8]<-apply(Reff_out[,1:7],1,max)
  Reff_out[,9]<-apply(Reff_out[,1:7],1,which.max)
  
  return(Reff_out)
}





get_csi_month<-function(out,mm,ri){
ncc<-(length(out)-1)
csi_all<-rep(NA,ncc)
for(i in 1 : ncc)
csi_all[i]<- out[[i]][[mm]]$csi[ri]
names(csi_all)<-names(out)[1:ncc]
return(csi_all)
}



get_csi_all<-function(out,ri){
ncc<-(length(out)-1)
csi_out<-matrix(NA,ncol=9,nrow=ncc)
colnames(csi_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_csi','Month')
rownames(csi_out)<-names(out)[1:ncc]
for(i in 1 : 7){
csi_out[,i]<-get_csi_month(out,i,ri)
}
csi_out[,8]<-apply(csi_out[,1:7],1,max)
csi_out[,9]<-apply(csi_out[,1:7],1,which.max)

return(csi_out)
}




get_csi_v_month<-function(out,mm,ri){
ncc<-(length(out)-1)
csi_all<-rep(NA,ncc)
for(i in 1 : ncc)
csi_all[i]<- out[[i]][[mm]]$csi_v[ri]
names(csi_all)<-names(out)[1:ncc]
return(csi_all)
}



get_csi_v_all<-function(out,ri){
ncc<-(length(out)-1)
csi_out<-matrix(NA,ncol=9,nrow=ncc)
colnames(csi_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_csi','Month')
rownames(csi_out)<-names(out)[1:ncc]
for(i in 1 : 7){
csi_out[,i]<-get_csi_v_month(out,i,ri)
}
csi_out[,8]<-apply(csi_out[,1:7],1,max)
csi_out[,9]<-apply(csi_out[,1:7],1,which.max)

return(csi_out)
}











get_crpss_month<-function(out,mm){
ncc<-(length(out)-1)
rl<-out[[ncc+1]]
crpss_all<-matrix(NA,ncol=2,nrow=ncc) 
 for(i in 1 : ncc) {
 crpss_all[i,1]<- out[[i]][[mm]]$crpss$crpss
 crpss_all[i,2]<- pt(crpss_all[i,1]/out[[i]][[mm]]$crpss$crpss.sigma,df=rl[i]-1,lower.tail=FALSE)
 
}
rownames(crpss_all)<-names(out)[1:ncc]
 return(crpss_all)
}

get_crpss_all<-function(out){
ncc<-(length(out)-1)
crpss_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(crpss_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_crpss','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','P_all')
rownames(crpss_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_crpss_month(out,mm=i)
crpss_out[,i]<-temp[,1]
crpss_out[,(9+i)]<-temp[,2]
}
crpss_out[,8]<-apply(crpss_out[,1:7],1,max)
crpss_out[,9]<-apply(crpss_out[,1:7],1,which.max)
crpss_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(crpss_out[i,crpss_out[i,9]+9])})
return(crpss_out)
}






get_crpss_v_month<-function(out,mm){
ncc<-(length(out)-1)
rl<-out[[ncc+1]]
crpss_all<-matrix(NA,ncol=2,nrow=ncc) 
 for(i in 1 : ncc) {
 crpss_all[i,1]<- out[[i]][[mm]]$crpss_v$crpss
 crpss_all[i,2]<- pt(crpss_all[i,1]/out[[i]][[mm]]$crpss_v$crpss.sigma,df=rl[i]-1,lower.tail=FALSE)
 
}
rownames(crpss_all)<-names(out)[1:ncc]
 return(crpss_all)
}

get_crpss_v_all<-function(out){
ncc<-(length(out)-1)
crpss_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(crpss_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_crpss','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','P_all')
rownames(crpss_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_crpss_v_month(out,mm=i)
crpss_out[,i]<-temp[,1]
crpss_out[,(9+i)]<-temp[,2]
}
crpss_out[,8]<-apply(crpss_out[,1:7],1,max)
crpss_out[,9]<-apply(crpss_out[,1:7],1,which.max)
crpss_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(crpss_out[i,crpss_out[i,9]+9])})
return(crpss_out)
}





get_roc_month<-function(out,mm,ri){
ncc<-(length(out)-1)
roc_out<-matrix(NA,ncol=2,nrow=ncc) 
rownames(roc_out)<-names(out)[1:ncc]
 for(i in 1 : ncc) {
 roc_out[i,1]<- out[[i]][[mm]]$roc_area[ri]
 roc_out[i,2]<- out[[i]][[mm]]$roc_pvalue[ri]
}
 return(roc_out)
}



get_roc_all<-function(out,ri){
ncc<-(length(out)-1)
roc_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(roc_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_roc','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','Sig_all')
rownames(roc_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_roc_month(out,mm=i,ri)
roc_out[,i]<-temp[,1]
roc_out[,(9+i)]<-temp[,2]
}
roc_out[,8]<-apply(roc_out[,1:7],1,max)
roc_out[,9]<-apply(roc_out[,1:7],1,which.max)
roc_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(roc_out[i,roc_out[i,9]+9])})
return(roc_out)
}





get_roc_v_month<-function(out,mm,ri){
ncc<-(length(out)-1)
roc_out<-matrix(NA,ncol=2,nrow=ncc) 
rownames(roc_out)<-names(out)[1:ncc]
 for(i in 1 : ncc) {
 roc_out[i,1]<- out[[i]][[mm]]$roc_area_v[ri]
 roc_out[i,2]<- out[[i]][[mm]]$roc_pvalue_v[ri]
}
 return(roc_out)
}



get_roc_v_all<-function(out,ri){
ncc<-(length(out)-1)
roc_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(roc_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_roc','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','Sig_all')
rownames(roc_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_roc_v_month(out,mm=i,ri)
roc_out[,i]<-temp[,1]
roc_out[,(9+i)]<-temp[,2]
}
roc_out[,8]<-apply(roc_out[,1:7],1,max)
roc_out[,9]<-apply(roc_out[,1:7],1,which.max)
roc_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(roc_out[i,roc_out[i,9]+9])})
return(roc_out)
}



get_bss_month<-function(out,mm){
ncc<-(length(out)-1)
rl<-out[[ncc+1]]
bss_out<-matrix(NA,ncol=2,nrow=ncc) 
rownames(bss_out)<-names(out)[1:ncc]
 for(i in 1 : ncc) {
 bss_out[i,1]<- out[[i]][[mm]]$brierss_qmean$bss
 bss_out[i,2]<- pt(bss_out[i,1]/out[[i]][[mm]]$brierss_qmean$bss.sigma,df=rl[i]-1,lower.tail=FALSE)
}
 return(bss_out)
}



get_bss_all<-function(out){
ncc<-(length(out)-1)
bss_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(bss_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_bss','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','P_all')
rownames(bss_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_bss_month(out,mm=i)
bss_out[,i]<-temp[,1]
bss_out[,(9+i)]<-temp[,2]
}
bss_out[,8]<-apply(bss_out[,1:7],1,max,na.rm=TRUE)
bss_out[,9]<-apply(bss_out[,1:7],1,which.max)
bss_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(bss_out[i,bss_out[i,9]+9])})
return(bss_out)
}






get_bss_v_month<-function(out,mm){
ncc<-(length(out)-1)
rl<-out[[ncc+1]]
bss_out<-matrix(NA,ncol=2,nrow=ncc) 
rownames(bss_out)<-names(out)[1:ncc]
 for(i in 1 : ncc) {
 bss_out[i,1]<- out[[i]][[mm]]$brierss_q20$bss
 bss_out[i,2]<- pt(bss_out[i,1]/out[[i]][[mm]]$brierss_q20$bss.sigma,df=rl[i]-1,lower.tail=FALSE)
}
 return(bss_out)
}



get_bss_v_all<-function(out){
ncc<-(length(out)-1)
bss_out<-matrix(NA,ncol=17,nrow=ncc)
colnames(bss_out)<-c('Jan','Feb','Mar','Apr','Mai','Jun','Jul','Max_bss','Month','P_Jan','P_Feb','P_Mar','P_Apr','P_Mai','P_Jun','P_Jul','P_all')
rownames(bss_out)<-names(out)[1:ncc]
for(i in 1 : 7){
temp<-get_bss_v_month(out,mm=i)
bss_out[,i]<-temp[,1]
bss_out[,(9+i)]<-temp[,2]
}
bss_out[,8]<-apply(bss_out[,1:7],1,max,na.rm=TRUE)
bss_out[,9]<-apply(bss_out[,1:7],1,which.max)
bss_out[,17]<-sapply(c(1:ncc),FUN=function(i){return(bss_out[i,bss_out[i,9]+9])})
return(bss_out)
}
















































analyse_all_swe<-function(npath,flood_values,qtrans,empqt=FALSE){
	out_analyze<-list()
	temp<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)) ,ncol=5,byrow=TRUE)[,1:2]
	for(i in 1 : dim(qtrans)[1]){
		out_analyze[[i]]<-list() 
		rnr<-as.integer(temp[i,1])
		hnr<-as.integer(temp[i,2])	
 #       print(temp[i,])		
		for( j in 1 : 7){
#		print(c(i,j))
            out_analyze[[i]][[j]]<-analyze_forecast_swe(rnr,hnr,j,fpath=npath,flood_values=flood_values,qtrans=NA,mplot=FALSE,empqt=FALSE)
		}
	}
	names(out_analyze)<-paste(temp[,1],'.',temp[,2],sep='')
    out_analyze$ndata<-qtrans[,2]
	return(out_analyze)
}











analyze_forecast_swe<-function(rnr,hnr,smonth,fpath="../data/netcdf/",flood_values,qtrans=NA,mplot=TRUE,empqt=FALSE){

if (!require('RNetCDF')) {
    stop('The package RNetCDF was not installed')
  }

  if (!require('nsRFA')) {
    stop('The package nsRFA was not installed')
  }

   if (!require('verification')) {
    stop('The package verification was not installed')
  }

  if (!require('SpecsVerification')) {
    stop('The package SpecsVerification was not installed')
  }  
  
  
ci<-which(flood_values[,3]==paste(rnr,'_',hnr,sep=''))
ndays<-c(31,28,31,30,31,30,31)
lmax<-sum(ndays[smonth:7])

nc<-RNetCDF::open.nc(paste(fpath,"seasonal_forecast_database_",rnr,'.',hnr,".nc",sep=''),main=flood_values[ci,3])
fyears<-RNetCDF::var.get.nc(nc,"Year")
fmonths<-RNetCDF::var.get.nc(nc,"Month")
scenario<-RNetCDF::var.get.nc(nc,"Scenario")
leadtime<-RNetCDF::var.get.nc(nc,"LeadTime")
qobs<-RNetCDF::var.get.nc(nc,"qobs")
qsim<-RNetCDF::var.get.nc(nc,"vf")
swe<-var.get.nc(nc,"swe")
precs<-var.get.nc(nc,"p")
temps<-var.get.nc(nc,"t")

#yi=which(syear==fyears)
mi=which(smonth==fmonths)
swe_start<-rep(NA,length(fyears))

maxvalues_obs<-rep(NA,length(fyears))
volumes_obs<-maxvalues_obs


for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
qobs_sel<-qobs[i,mi,1:lmax]
if ( any(is.na(qobs_sel)))volumes_obs[i]<-NA
else volumes_obs[i]<-sum(qobs_sel,na.rm=TRUE)*60*60*24/1000000

swe_start[i]<-swe[i,mi,temps[i,mi,,1]<0 & precs[i,mi,,1]<=0.0001,1][1]
qobs_sel<-qobs[i,mi,1:lmax]
if ( any(is.na(qobs_sel))){
maxvalues_obs[i]<-NA
}
else {
maxvalues_obs[i]<-max(qobs_sel,na.rm=TRUE)
}




}

#  xx <- seq(min(x,na.rm=T),max(x,na.rm=T), length=length(unique_years))
#  lines(xx, predict(fit3, data.frame(x=xx)), col="black")
#predict(out$fit, data.frame(swe_start=50))

windows(6,6)
plot(swe_start,maxvalues_obs,type="p",xlab="SWE (mm)",ylab="Flood peak",main=flood_values[ci,3])
myfit1<-lm(maxvalues_obs~swe_start)
xx <- seq(min(swe_start,na.rm=T),max(swe_start,na.rm=T), length=length(swe_start))
lines(xx, predict(myfit1, data.frame(swe_start=xx)), col="chartreuse4")

windows(6,6)
plot(swe_start,volumes_obs,type="p",xlab="SWE (mm)",ylab="Flood volume",main=flood_values[ci,3])
myfit2<-lm(volumes_obs~swe_start)
xx <- seq(min(swe_start,na.rm=T),max(swe_start,na.rm=T), length=length(swe_start))
lines(xx, predict(myfit2, data.frame(swe_start=xx)), col="chartreuse4")


ssmin=min(maxvalues_obs,myfit1$fitted.values,na.rm=T)
ssmax=max(maxvalues_obs,myfit1$fitted.values,na.rm=T)
windows(8,8)
plot(na.omit(maxvalues_obs),myfit1$fitted.values,xlim=c(ssmin,ssmax),ylim=c(ssmin,ssmax),xlab="Observed floods (m3/s)",ylab="Regression model (m3/s)")
legend('bottomright',legend=c(paste("Regression model cor:",round(cor(na.omit(maxvalues_obs),myfit1$fitted.values,use="pairwise.complete.obs"),2))),col=c(1),pch=c(1))
abline(0,1)

ssmin=min(volumes_obs,myfit2$fitted.values,na.rm=T)
ssmax=max(volumes_obs,myfit2$fitted.values,na.rm=T)
windows(8,8)
plot(na.omit(volumes_obs),myfit2$fitted.values,xlim=c(ssmin,ssmax),ylim=c(ssmin,ssmax),xlab="Observed volume (M m3)",ylab="Regression model (m3/s)")
legend('bottomright',legend=c(paste("Regression model cor:",round(cor(na.omit(volumes_obs),myfit2$fitted.values,use="pairwise.complete.obs"),2))),col=c(1),pch=c(1))
abline(0,1)
out<-list()
out$corr<-cor(na.omit(maxvalues_obs),myfit1$fitted.values,use="pairwise.complete.obs")
out$corr_v<-cor(na.omit(volumes_obs),myfit2$fitted.values,use="pairwise.complete.obs")
out$Reff<-1.0-(mean((myfit1$fitted.values-na.omit(maxvalues_obs))^2,na.rm=TRUE)/var(maxvalues_obs,na.rm=TRUE))
out$Reff_v<-1.0-(mean((myfit2$fitted.values-na.omit(volumes_obs))^2,na.rm=TRUE)/var(volumes_obs,na.rm=TRUE))
out$fit<-myfit1
return(out)

}







