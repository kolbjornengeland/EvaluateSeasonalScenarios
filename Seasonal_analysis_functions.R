

plot_forecast<-function(rnr,hnr,syear,smonth,flood_values,ptype="greyshade", fpath="../data/netcdf/",
qtrans=NA){

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

if(any(!is.na(qtrans))){
simdim<-dim(qsim_sel)
qsim_sel_t<-qsim_sel
tnumbers<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)),ncol=5,byrow=TRUE)
i2<-which(as.integer(tnumbers[,1])==rnr & as.integer(tnumbers[,2])==hnr)
parsel<-as.numeric(qtrans[i2,3:6])

for(m in 1:simdim[1]){
for(n in 1:simdim[2]){
FF1<-nsRFA::F.gumb(qsim_sel[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-nsRFA::invF.gumb (FF1,parsel[1], parsel[2])
}
}

qsim_sel<-qsim_sel_t
}



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

if(any(is.na(qtrans))){
abline(h=flood_values[ci,7],col="yellow",lty=2)
abline(h=flood_values[ci,8],col="orange",lty=2)
abline(h=flood_values[ci,9],col="red",lty=2)
}


nover<-c(1:3)
for(j in 1 : 3){
if(any(is.na(qtrans))){
ttest<-qsim_sel>flood_values[ci,6+j]
}
else{
ttest<-qsim_sel>flood_values[ci,3+j]
}

nover_t<-rowSums(ttest,na.rm=TRUE)
nover_t[nover_t>=1]=1
nover[j]=mean(nover_t)

if(any(is.na(qtrans))){
text(10,flood_values[ci,6+j],round(nover[j],2))
}
else{
text(10,flood_values[ci,3+j],round(nover[j],2))
}


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





analyze_forecast<-function(rnr,hnr,smonth,fpath="../data/netcdf/",flood_values,qtrans=NA,mplot=TRUE){
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

#Pthresholds<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
#False_alarm_rate<-matrix(ncol=length(Pthresholds),nrow=3)
#Hit_rate<-False_alarm_rate
#csi<-False_alarm_rate

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



if(any(!is.na(qtrans))){
simdim<-dim(maxvalues_sim)
qsim_sel_t<-maxvalues_sim
tnumbers<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)),ncol=5,byrow=TRUE)
i2<-which(as.integer(tnumbers[,1])==rnr & as.integer(tnumbers[,2])==hnr)
parsel<-as.numeric(qtrans[i2,3:6])

for(m in 1:simdim[1]){
for(n in 1:simdim[2]){
FF1<-nsRFA::F.gumb(maxvalues_sim[m,n],parsel[3],parsel[4])
qsim_sel_t[m,n]<-nsRFA::invF.gumb (FF1,parsel[1], parsel[2])
}
}

maxvalues_sim<-qsim_sel_t
}

for(i in 1 : length(fyears)){
si=which(fyears[i]!=scenario)
ens.ref[i,]<-maxvalues_obs[si]
}



#maxvalues_obs<-maxvalues_obs[maxvalues_sim[,1]>0]
#ens.ref<-ens.ref[maxvalues_sim[,1]>0,]
#maxvalues_sim<-maxvalues_sim[maxvalues_sim[,1]>0,]


roc_area<-c()
roc_pvalue<-c()
for(j in 1 : 3){

outcome<-as.integer(maxvalues_obs>flood_values[ci,(3+j)])
if(sum(na.omit(outcome))==0){
roc_area[j]=0
roc_pvalue[j]=0
}
else{
p_flood<-rowMeans(maxvalues_sim>flood_values[ci,(6+j)])
p_flood<-p_flood[!is.na(outcome)]
outcome<-na.omit(outcome)
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
}

qsim_95<-apply(maxvalues_sim,1,quantile,0.95)
qsim_5<-apply(maxvalues_sim,1,quantile,0.05)
qsim_50<-apply(maxvalues_sim,1,quantile,0.5)
qsim_max<-apply(maxvalues_sim,1,max,na.rm=TRUE)

ssmax<-max(qsim_max,maxvalues_obs)
ssmin<-min(qsim_5,maxvalues_obs)

if(mplot){
windows(8,8)
plot(maxvalues_obs,qsim_50,xlim=c(ssmin,ssmax),ylim=c(ssmin,ssmax),xlab="Observed floods (m3/s)",ylab="Scenario (m3/s)")
points(maxvalues_obs,qsim_max,col="red")
points(maxvalues_obs,qsim_95,col="orange")
points(maxvalues_obs,qsim_5,col="green")
legend('bottomright',legend=c(paste("Median scenario cor:",round(cor(maxvalues_obs,qsim_50),2)),
paste("Maximum scenario cor:",round(cor(maxvalues_obs,qsim_max),2))
,paste("95% scenariocor:",round(cor(maxvalues_obs,qsim_95),2))
,paste("5% scenario cor:",round(cor(maxvalues_obs,qsim_5),2)))
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

out$rankhist<-SpecsVerification::Rankhist(maxvalues_sim, maxvalues_obs)
if(mplot){
windows(6,4)
SpecsVerification::PlotRankhist(out$rankhist)
}
out$sim_max<-maxvalues_sim
out$obs_max<-maxvalues_obs
out$crps<-SpecsVerification::EnsCrps(maxvalues_sim, maxvalues_obs)
out$brier_qmean<-SpecsVerification::EnsBrier(maxvalues_sim, maxvalues_obs,flood_values[ci,7])
out$brier_q5<-SpecsVerification::EnsBrier(maxvalues_sim, maxvalues_obs,flood_values[ci,8])
out$brier_q50<-SpecsVerification::EnsBrier(maxvalues_sim, maxvalues_obs,flood_values[ci,9])
qsim_50<-apply(maxvalues_sim,1,quantile,0.5)
out$rmse<-sqrt(mean((qsim_50-maxvalues_obs)^2))
out$Reff<-1.0-(mean((qsim_50-maxvalues_obs)^2)/var(maxvalues_obs))
out$corr<-cor(maxvalues_obs,qsim_50,use="pairwise.complete.obs")
out$crpss<-SpecsVerification::EnsCrpss(maxvalues_sim, ens.ref, maxvalues_obs) 
out$brierss_qmean<-SpecsVerification::EnsBrierSs(maxvalues_sim, ens.ref, maxvalues_obs,flood_values[ci,7])
out$brierss_q5<-SpecsVerification::EnsBrierSs(maxvalues_sim, ens.ref, maxvalues_obs,flood_values[ci,8])
out$brierss_q50<-SpecsVerification::EnsBrierSs(maxvalues_sim, ens.ref, maxvalues_obs,flood_values[ci,9])

out
}


analyse_all<-function(npath,flood_values,qtrans){
	out_analyze<-list()
	for(i in 1 : dim(qtrans)[1]){
		out_analyze[[i]]<-list() 
		for( j in 1 : 7){
		print(c(i,j))
			temp<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)) ,ncol=5,byrow=TRUE)
			rnr<-as.integer(temp[i])
			hnr<-as.integer(temp[i,2])
			out_analyze[[i]][[j]]<-analyze_forecast(rnr,hnr,j,fpath=npath,flood_values=flood_values,qtrans=qtrans,mplot=FALSE)
		}
	}
	return(out_analyze)
}

#get_corr_month(seasonal_evaluation,qtransform_sel[,2],3)


get_corr_month<-function(out,ndata,mm){
corr_all<-rep(NA,length(out))
for(i in 1 : length(out))
corr_all[i]<- out[[i]][[mm]]$corr
# Beregner t-value:
tt<-abs(qt(0.05,ndata-2))
#beregner kritisk verdi for korrelasjon
rr= tt/sqrt(ndata-2+tt^2)
return(cbind(corr_all,rr))
}




get_crpss_month<-function(out,mm){

crpss_all<-rep(NA,dim(qtransform_sel)[1]) 

 for(i in 1 : dim(qtransform_sel)[1]) 
 
 crpss_all[i]<- out_analyze[[i]][[mm]][[1]]$crpss$crpss 

crpss_all 

}

get_corr_month<-function(out,mm){
corr_all<-rep(NA,dim(qtransform_sel)[1])
for(i in 1 : dim(qtransform_sel)[1])
corr_all[i]<- out_analyze[[i]][[mm]][[1]]$corr
corr_all
}

get_rmse_month<-function(out,mm){
rmse_all<-rep(NA,dim(qtransform_sel)[1])
for(i in 1 : dim(qtransform_sel)[1])
rmse_all[i]<- out_analyze[[i]][[mm]][[1]]$rmse
rmse_all
}

get_brierss_qmean_month<-function(out,mm){
brierss_all<-rep(NA,dim(qtransform_sel)[1])
for(i in 1 : dim(qtransform_sel)[1])
brierss_all[i]<- out_analyze[[i]][[mm]][[1]]$brierss_qmean$bss
brierss_all
}

get_csi<-function(out,mm){
csi_all<-rep(NA,dim(qtransform_sel)[1])
for(i in 1 : dim(qtransform_sel)[1])
csi_all[i]<- out_analyze[[i]][[mm]][[1]]$csi
csi_all
}

get_ROC_mean<-function(out,mm){
ROC_mean_all<-rep(NA,dim(qtransform_sel)[1])
for(i in 1 : dim(qtransform_sel)[1])
ROC_mean_all[i]<- out_analyze[[i]][[mm]][[1]]$roc_area[1]
ROC_mean_all
}
