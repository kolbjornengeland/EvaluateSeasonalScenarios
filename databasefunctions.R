
create_netcdffiles<-function(filepath, qobspath, outpath,fyears=c(1962:2015),fvariables=c("vf","t","swe","sm","p","gv","e")){
if (!require('RNetCDF')) {
    stop('The package RNetCDF was not installed')
  }

#Read the content of the data folder
stations<-list.files(filepath)

#Defines the dimensions
dim.stations<-length(stations)-2
dim.years<-length(fyears)
fmonths<-seq(1:7)
dim.months<-length(fmonths)

dim.variables<-length(fvariables)

alldates<-seq(as.Date("1962/1/1"), as.Date("2015/12/31"), "days")
dim.dates<-length(alldates)
dim.days<-366
dim.characters <- 64


shn<-strsplit(stations,'.')
stemp<-unlist(strsplit(stations,".",fixed=TRUE))
rnr<-stemp[seq(from=1, to=dim.stations*2-1,by=2)]
hnr<-stemp[seq(from=2, to=dim.stations*2,by=2)]



for(i in 1 : dim.stations){
nc <- create.nc(paste(outpath,"/seasonal_forecast_database_",stations[i],".nc",sep=''))  # CHECK DIR
create_database(nc)


qfile<-paste(qobspath,'/q.',rnr[i],'~1.',hnr[i],sep='')
qtemp<-read.table(qfile,header=FALSE,sep="")
qyear<-as.integer(substr(qtemp[,1],1,4))
qmonth<-as.integer(substr(qtemp[,1],5,6))


for(j in 1 : dim.variables){
for (k in 1 : dim.years){
for (l in 1 : dim.months){

if(j==1){
qsel<-qtemp[qyear==fyears[k]&qmonth>=fmonths[l],2]
        var.put.nc(nc,"qobs", data = qsel, start=c(k, l, 1), 
                   count=c(1, 1, length(qsel)))
}

# Read the forecast
filename<-paste(filepath,'/',stations[i],'/res.',fvariables[j],'.',stations[i],'.',fyears[k],'0',fmonths[l],'01',sep='')
data.temp<-t(read.table(filename,header=TRUE,row.names=1))

        var.put.nc(nc, fvariables[j], data = data.temp, start=c(k, l, 1, 1), 
                   count=c(1, 1, dim(data.temp)))
				  
sync.nc(nc)
}
}
}
close.nc(nc)
}
}




create_database<-function(nc){
# Function to create the netcdf-filename
# It has four dimensjons for storing the forecasts and one extra dimensjon for storing text strings.
# Each forecast is defined by its issue date that for the scenario is the first of the months january to july
# The dimensjon "year" describes the year of the issue date, and the dimension "month" describes the month of the issue date.
# Each forecast has 54 scenarious. The scenario is the third dimensjon
# Each forecast has lead time that is up to 366 days. Lead time is the fourth dimensjon

att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", "Seasonal forecasts")
att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR", paste("Created on", base::date()))

#dim.def.nc(nc, "station", dim.stations)
dim.def.nc(nc, "year", dim.years)  
dim.def.nc(nc, "month", dim.months)
#dim.def.nc(nc, "dates", dim.dates)
dim.def.nc(nc, "max_string_length", dim.characters)
dim.def.nc(nc,"day",dim.days)
dim.def.nc(nc,"scenario",dim.years)


# Creating variables that describes the dimensions in the database
var.def.nc(nc, varname = "Month", vartype = "NC_FLOAT", dimensions = c("month"))
att.put.nc(nc, "Month", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "Month", "short_name", "NC_CHAR", "Month")
att.put.nc(nc, "Month", "long_name", "NC_CHAR", "Month for the issue date")
var.put.nc(nc, "Month", fmonths)
sync.nc(nc)


# Creating variables that describes the dimensions in the database
var.def.nc(nc, varname = "Year", vartype = "NC_FLOAT", dimensions = c("year"))
att.put.nc(nc, "Year", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "Year", "short_name", "NC_CHAR", "Month")
att.put.nc(nc, "Year", "long_name", "NC_CHAR", "Year for the issue date")
var.put.nc(nc, "Year", fyears)
sync.nc(nc)

# Creating variables that describes the dimensions in the database
var.def.nc(nc, varname = "Scenario", vartype = "NC_FLOAT", dimensions = c("scenario"))
att.put.nc(nc, "Scenario", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "Scenario", "short_name", "NC_CHAR", "Scenario")
att.put.nc(nc, "Scenario", "long_name", "NC_CHAR", "Year used to create the scenario")
var.put.nc(nc, "Scenario", fyears)
sync.nc(nc)

# Creating variables that describes the dimensions in the database
var.def.nc(nc, varname = "LeadTime", vartype = "NC_FLOAT", dimensions = c("day"))
att.put.nc(nc, "LeadTime", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "LeadTime", "short_name", "NC_CHAR", "Lead time")
att.put.nc(nc, "LeadTime", "long_name", "NC_CHAR", "Lead time. Forecast is always issued the first each month")
var.put.nc(nc, "LeadTime", fyears)
sync.nc(nc)


temp.dat<-array(NA,dim=c(dim.years,dim.months,dim.years,dim.days))
temp2.dat<-array(NA,dim=c(dim.years,dim.months,dim.days))

var.def.nc(nc, varname = "vf", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "vf", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "vf", "short_name", "NC_CHAR", "Streamflow scenario")
att.put.nc(nc, "vf", "long_name", "NC_CHAR", "Streamflow scenario for seasonal forecasting (m3/s)")
var.put.nc(nc, "vf", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "qobs", vartype = "NC_FLOAT", dimensions = c("year","month","day"))
att.put.nc(nc, "qobs", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "qobs", "short_name", "NC_CHAR", "Observed streamflow")
att.put.nc(nc, "qobs", "long_name", "NC_CHAR", "Observed streamflow (m3/s)")
var.put.nc(nc, "qobs", temp2.dat)
sync.nc(nc)

var.def.nc(nc, varname = "e", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "e", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "e", "short_name", "NC_CHAR", "Evapotranspiration scenario")
att.put.nc(nc, "e", "long_name", "NC_CHAR", "Evapotranspiration scenario for seasonal forecasting (mm/day)")
var.put.nc(nc, "e", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "swe", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "swe", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "swe", "short_name", "NC_CHAR", "SWE scenario")
att.put.nc(nc, "swe", "long_name", "NC_CHAR", "SWE scenario for seasonal forecasting (mm)")
var.put.nc(nc, "swe", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "sm", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "sm", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "sm", "short_name", "NC_CHAR", "Soil moisture scenario")
att.put.nc(nc, "sm", "long_name", "NC_CHAR", "Soil moisture  scenario for seasonal forecasting (mm)")
var.put.nc(nc, "sm", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "gv", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "gv", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "gv", "short_name", "NC_CHAR", "Groundwater scenario")
att.put.nc(nc, "gv", "long_name", "NC_CHAR", "Groundwater scenario for seasonal forecasting (mm)")
var.put.nc(nc, "gv", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "p", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "p", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "p", "short_name", "NC_CHAR", "Precipiation scenario")
att.put.nc(nc, "p", "long_name", "NC_CHAR", "Precipiation scenario for seasonal forecasting (mm/day)")
var.put.nc(nc, "p", temp.dat)
sync.nc(nc)

var.def.nc(nc, varname = "t", vartype = "NC_FLOAT", dimensions = c("year","month","scenario","day"))
att.put.nc(nc, "t", "missing_value", "NC_FLOAT", -9999)
att.put.nc(nc, "t", "short_name", "NC_CHAR", "Temperature scenario")
att.put.nc(nc, "t", "long_name", "NC_CHAR", "Temperature scenario for seasonal forecasting (mm)")
var.put.nc(nc, "t", temp.dat)
sync.nc(nc)


}

