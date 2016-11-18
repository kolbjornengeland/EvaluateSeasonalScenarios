

create_transformations<-function(HBVsim_path='inst/Hydra2_6004',Qobs_path='inst/Disc',outfile='inst/Transformations.txt',tdist='gumbel') {
  if (!require('nsRFA')) {
    stop('The package nsRFA was not installed')
  }
  
    if (!require('fitdistrplus')) {
    stop('The package fitdistrplus was not installed')
  }
  
  
#filepath<-('E:/SveinTore/Hydra2_6004')
stations<-list.files(HBVsim_path)

dim.stations<-length(stations)

shn<-strsplit(stations,'.')

stemp<-unlist(strsplit(stations,".",fixed=TRUE))
rnr<-stemp[seq(from=1, to=dim.stations*5-4,by=5)]
hnr<-stemp[seq(from=2, to=dim.stations*5-3,by=5)]

nobs<-rep(NA,dim.stations)
param_qsim<-matrix(ncol=2,nrow=dim.stations)
param_qobs<-matrix(ncol=2,nrow=dim.stations)

for(i in 1 : dim.stations){
print(i)
#print(stations[i])
qsim<-read.table(paste(HBVsim_path,'/',stations[i],sep=''),header=FALSE,skip=8)
#print(dim(qsim))
years<-substr(qsim[,1],1,4)
qfile<-paste(Qobs_path,'/q.',rnr[i],'~1.',hnr[i],sep='')
qobs<-read.table(qfile,header=FALSE)
# Remove 31 january 2015 since it is not in the simulated dataset
qobs<-qobs[qobs[,1]!='20151231/1200',]

# Use only positive values
qobs<-qobs[qobs[,2]>=0.0,]
qtest<-qobs[1:2000,]
si<-match(qobs[,1],qsim[,1])
qobs<-qobs[!is.na(si),]
qsim_emp<-qsim[as.numeric(na.omit(si)),]

#qsim_emp<-qsim[qobs[,2]>0.0,]
#qobs<-qobs[qobs[,2]>0.0,]
years_obs<-substr(qobs[,1],1,4)
years_ll<- as.numeric(by(qobs[,2],years_obs,length))
years_u<-unique(years_obs)
years_sel<-years_u[years_ll>364]
qobs<-qobs[years_obs%in%years_sel,]
qsim_emp<-qsim_emp[years_obs%in%years_sel,]
years_obs<-substr(qobs[,1],1,4)
nobs[i]<-length(years_sel)

if(tdist=='gumbel'){
ams_qsim<- by(qsim[,2],years,max)
ams_qsim_lmom <- Lmoments(ams_qsim)
param_qsim[i,] <- as.numeric(par.gumb(ams_qsim_lmom[1],ams_qsim_lmom[2]))
ams_qobs<- by(qobs[,2],years_obs,max)
ams_qobs_lmom <- Lmoments(ams_qobs)
param_qobs[i,] <- as.numeric(par.gumb(ams_qobs_lmom[1],ams_qobs_lmom[2]))
}

if(tdist=='gamma'){
gfit<-fitdist(qsim[,2], "gamma",method='qme',probs=c(0.7,0.9999))
param_qsim[i,] <- as.numeric(gfit$estimate)
gfit<-fitdist(qobs[,2], "gamma")
param_qobs[i,] <- as.numeric(gfit$estimate,method='qme',probs=c(0.7,0.9999))
}

if(tdist=='empq'){
dd<-dim(qsim_emp)[1]-dim(qobs)[1]
if(dd!=0)print(dd)
#if(rnr[i]==2 & hnr[i]==28){
#print(dim(qsim_emp))
#print(dim(qobs))
#}
xtab<-cbind(sort(qsim_emp[,2]),sort(qobs[,2]))
ofile<-paste("inst/empq/",rnr[i],'.',hnr[i],".txt",sep='')
#print(ofile)

  write.table(xtab,file=ofile,sep=";",col.names=c("qsim","qobs"),row.names=FALSE)

}
}

write.table(cbind(stations,nobs,param_qobs,param_qsim,tdist),file=outfile)

}
