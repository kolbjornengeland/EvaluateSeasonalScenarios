

create_transformations<-function(HBVsim_path='inst/Hydra2_6004',Qobs_path='inst/Disc',outfile='inst/Transformations.txt') {
  if (!require('nsRFA')) {
    stop('The package nsRFA was not installed')
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

qsim<-read.table(paste(HBVsim_path,'/',stations[i],sep=''),header=FALSE,skip=8)
years<-substr(qsim[,1],1,4)
ams_qsim<- by(qsim[,2],years,max)
ams_qsim_lmom <- Lmoments(ams_qsim)
param_qsim[i,] <- as.numeric(par.gumb(ams_qsim_lmom[1],ams_qsim_lmom[2]))

qfile<-paste(Qobs_path,'/q.',rnr[i],'~1.',hnr[i],sep='')
qobs<-read.table(qfile,header=FALSE)
qobs<-qobs[qobs[,2]>0.0,]
years_obs<-substr(qobs[,1],1,4)
ams_qobs<- by(qobs[,2],years_obs,max)
ams_ll<- as.numeric(by(qobs[,2],years_obs,length))
ams_qobs<-ams_qobs[ams_ll>364]
nobs[i]<-length(ams_qobs)

ams_qobs_lmom <- Lmoments(ams_qobs)
param_qobs[i,] <- as.numeric(par.gumb(ams_qobs_lmom[1],ams_qobs_lmom[2]))

}

write.table(cbind(stations,nobs,param_qobs,param_qsim),file=outfile)

}
