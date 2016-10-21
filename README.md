# EvaluateSeasonalScenarios

R package for evaluation seasonal scenarios givn by the HBV model at NVE. This includes to 

## Installation

Download the *.R files into your own directory and 'source' them:

```R
setwd('C:/Users/koe/Documents/EvaluateSeasonalScenarios/')
source('databasefunctions.R')
source('Seasonal_analysis_functions.R')
```

You will need the RNetCDF package:

```R
install('RNetCDF')
install('nsRFA')
install('SpecsVerification')
install('verification')
```R

The data used here is available upon request

## Preprosessing of files: Creating NCDF files from original text files
The first step is to create a smalle number of NetCDF files form the original ascii files. 
filepath is here the path where the original files are stored
qobspath is path to where the observed streamflows are stored
outpath is here the path where the NetCDF files should be written
```R
filepath<-('X:/Hbv/ScenarioerFou/vf')
qobspath<-c('inst/Disc')
outpath<-('M:/Dokumenter/Sesongvarsler')
create_netcdffiles(filepath,outpath)
```

## Preprosessing of files: Creating Transformations for flood values
Fits gumbel distributions to annual maxima of observed and simulated streamflow. Write the parameters to file 
HBVsim_path is here the path to folder with HBV simulations of streamflow
Qobs_path is path to where the observed streamflows are stored
outfile is here the name of the output file

```R
HBV='inst/Hydra2_6004'
Qobs='inst/Disc'
outf='inst/Transformations.txt'
create_transformations(HBV,Qobs,outf)
```


## Plotting scenarios
Plots the seasonal scenarios for a selected catchment, for a spcific issue year and month. Can choos between a spaghetti plot and grey-shade plot 

flood_values is a file with estimated mean, 5-years and 50 years floods for observations and for HBV simulations
qtransform is a file with parameterw for using a gumbel distributio to transform simulated floods to become mor similar to observed floods
NetCDFfolder is the path to where the forecast NetCDF files are stored
Rnr is Regine area
Hnr is main number
Syear is year


```R
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")
NetCDFfolder="M:/Dokumenter/Sesongvarsler/"
Rnr=2
Hnr=11
Syear=1974
Smonth=4
plot_forecast(Rnr,Hnr,Syear,Smonth,flood_values,ptype="greyshade",fpath=NetCDFfolder)
plot_forecast(Rnr,Hnr,Syear,Smonth,flood_values,ptype="greyshade",fpath=NetCDFfolder,qtrans=qtransform)
plot_forecast(Rnr,Hnr,Syear,Smonth,flood_values,ptype="spaghetti",fpath=NetCDFfolder,qtrans=qtransform)
```

## Analysing scenarios for one specific catchment issued at one specific month
Evaluate the performance of the scenariou using standard verification tools and verification metricsRnr

Returns: 
1. crps, crpss
2. For mean 5-years and 50 years floods the Brier score, the Brier skill score and the, ROC-area 
3. For the median ensemble: rmse, correlation and Reff
		
flood_values is a file with estimated mean, 5-years and 50 years floods for observations and for HBV simulations		
NetCDFfolder is the path to where the forecast NetCDF files are stored
Rnr is Regine area
Hnr is main number
Syear is year
Qobs_path is path to where the observed streamflows are stored
outfile is here the name of the output file

```R
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")
NetCDFfolder="M:/Dokumenter/Sesongvarsler/"
Rnr=2
Hnr=11
Smonth=4
analyze_forecast(Rnr,Hnr,Smonth,fpath=NetCDFfolder,flood_values,qtrans=NA,mplot=TRUE)
```



NetCDFfolder="M:/Dokumenter/Sesongvarsler/"
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")

analyse_all(NetCDFfolder,flood_values,qtransform)




