# EvaluateSeasonalScenarios

R package for evaluation seasonal scenarios givn by the HBV model at NVE. This includes to 

## Installation

Download the *.R files into your own directory and 'source' them:

```R
setwd'C:/Users/koe/Documents/EvaluateSeasonalScenarios/')
source('databasefunctions.R')
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
qobspath<-c('X:/Hbv/ScenarioerFou/data/Disc')
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
Plots the seasonal scenarios for a selected catchment, year and month
H
Qobs_path is path to where the observed streamflows are stored
outfile is here the name of the output file

```R
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")

Rnr=2
Hnr=11
Syear
Smonth=4
plot_forecast(Rnr,Hnr,Syear,Smonth,flood_values,ptype="greyshade",fpath="../data/netcdf/")

```