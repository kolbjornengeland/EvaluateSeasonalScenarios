



# Read metadata from file:
ccp<-read.table("Metadata.txt",sep=";",header=TRUE)

# order the catchemnt characteristics to be in the sam order as the output fromthe analyses
catchment_properties<-ccp[match(names(seasonal_evaluation)[1:(length(seasonal_evaluation)-1)],paste(ccp[,1],'.',ccp[,2],sep='')),]



# Exctracting climate characteristics
Mprec<-read.table('inst/metdata/SaveP.txt',header=TRUE,row.names=1,check.names = FALSE)
Mrain<-read.table('inst/metdata/SaveR.txt',header=TRUE,row.names=1,check.names = FALSE)
Mtemp<-read.table('inst/metdata/SaveT.txt',header=TRUE,row.names=1,check.names = FALSE)
NSubZero<-read.table("inst/metdata/SubZero.txt",check.names = FALSE) 


# Annual values
Aprec<-colMeans(Mprec)
Arain<-colMeans(Mrain)
Atemp<-colMeans(Mtemp)
ASZ<-colSums(NSubZero)

# Ratio rain precipitation
ARRatio<-Arain/Aprec

# Winter values (November - April)
Wprec<-colMeans(Mprec[c(1,2,3,4,11,12),])
Wrain<-colMeans(Mrain[c(1,2,3,4,11,12),])
Wtemp<-colMeans(Mtemp[c(1,2,3,4,11,12),])
WSZ<-colSums(NSubZero[c(1,2,3,4,11,12),])/(30+31+31+28.25+31+30)

# Ratio rain precipitation
WRRatio<-Wrain/Wprec






##########################################################################################

Scripts below are for preparing metadata

##########################################################################################


# Reading daily values
Dtemp<-read.table('inst/metdata/DailyValues/aveT.txt',header=TRUE,row.names=1,check.names = FALSE)
Mdates<-rownames(Dtemp)
Myears<-substr(Mdates,1,4)
Mmonth<-substr(Mdates,6,7)

# Number of days with sub-zero average temperature
Dtemp_binary<-(Dtemp<=0.0)
Nmonthly<-apply(Dtemp_binary,2,FUN=by,Mmonth,sum)/length(unique(Myears))
write.table(Nmonthly,"inst/SubZero.txt") 







#Extract metadata for stations:

# installering av pakker
install.packages("devtools")
install.packages("curl")
library(devtools)
install_github("NVE/NVEDATA")

install_github("NVE/NVEHYDROTOOLS")

setwd('C:/Users/koe/Documents/EvaluateSeasonalScenarios/')
source('databasefunctions.R')
source('Seasonal_analysis_functions.R')
source('creating_transformations.R')

# The following code returns metadata for all available stations:

library(NVEDATA)
metadata <- get_metadata()



# Get the metadata for the flood forecasting stations
qtrans<-read.table("inst/Transformations_gumbel.txt")
temp<-matrix(unlist(strsplit(as.character(qtrans[,1]),'.',fixed=TRUE)) ,ncol=5,byrow=TRUE)[,1:2]

mms<-metadata[,1]*100000+metadata[,2]
ffs<-as.integer(temp[,1])*100000+as.integer(temp[,2])

fi<-match(ffs,as.integer(unlist(mms)))
 
metadata_flood_stations<-metadata[fi,]

# List the available metadata
names(metadata_flood_stations)

# List the regulated volum (reservoir volume / annual ruoff ) 
reguleringsgrad<-metadata_flood_stations$regulation_part_reservoirs 


library(NVEHYDROTOOLS)
example_stations<-paste(temp[,1],".",temp[,2],".0",sep="")
shapef <- '//nve/fil/h/HM/Interne Prosjekter/Flomkart/Data/GISData/Hydrologi_TotalNedborfeltMalestasjon.shp'
slayer<- 'Hydrologi_TotalNedborfeltMalestasjon'
snr_t="C:/Users/koe/Documents/Flomkart/NVEHYDROTOOLS/inst/Complete_data/CatchmentCharacteristics/Feltnr_flomkart_til_feltnr_GIS.txt"
outfile<-'inst/CID.txt'
grid_id_example_catchments<-gridcell_list(c_ids=example_stations,c_shape=shapef,snr_translation=snr_t,c_layer=slayer,outfile=outfile)



grid_id<-'inst/CID.txt'
outf<-'inst/metdata/'
metinf<-get_metdataforfloods(gridid=grid_id,first_day=as.Date("1962/1/1"),last_day=as.Date("2015/12/31"),
station_file=NA,metfolder="U:/metdata/met_obs_v2.1/",snowfolder="U:/snowsim/snowsim_v2.0.1/",hbvfolder="Z:/gwbsim/",outfolder=outf)
