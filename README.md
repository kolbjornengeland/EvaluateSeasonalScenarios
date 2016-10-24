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
filepath<-'X:/Hbv/ScenarioerFou/vf'
qobspath<-'inst/Disc'
outpath<-'M:/Dokumenter/Sesongvarsler'
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
Evaluate the performance of the scenarios using standard verification tools and verification metricsRnr

Returns: 
1. crps, crpss
2. For mean 5-years and 50 years floods the Brier score, the Brier skill score and the, ROC-area 
3. For the median ensemble: rmse, correlation and Reff


Rnr is Regine area
Hnr is main number
Smonth is the issue month		
flood_values is a matrix with estimated mean, 5-years and 50 years floods for observations and for HBV simulations		
qtrans is a matrix with the parameters for transforming the scenarios. Created by 'create_transformations'
fpath is the path to where the forecast NetCDF files are stored
mplot is TRUE for showing evaluation plots.
 
```R
Rnr=2
Hnr=11
Smonth=4
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")
NetCDFfolder="M:/Dokumenter/Sesongvarsler/"

out<-analyze_forecast(Rnr,Hnr,Smonth,fpath=NetCDFfolder,flood_values,qtrans=NA,mplot=TRUE)
```


## Analysing scenarios for all catchments and all months
Evaluate the performance of the scenarios using standard verification tools and verification metrics

Returns: 
1. crps, crpss
2. For mean 5-years and 50 years floods the Brier score, the Brier skill score and the, ROC-area 
3. For the median ensemble: rmse, correlation and Reff



Rnr is Regine area
Hnr is main number
Smonth is the issue month		
flood_values is a matrix with estimated mean, 5-years and 50 years floods for observations and for HBV simulations		
qtrans is a matrix with the parameters for transforming the scenarios. Created by 'create_transformations'
fpath is the path to where the forecast NetCDF files are stored
mplot is TRUE for showing evaluation plots.
 


```R
NetCDFfolder="M:/Dokumenter/Sesongvarsler/"
flood_values<-read.table('inst/flomtabell_ny.txt',sep="")
qtransform<-read.table("inst/Transformations.txt")
qtransform_sel<-qtransform[qtransform[,2]>30,]
# Hack. need to exclude 22.16, 41.1 and 55.5 for the moment. Problems with missing data for scenarios,
qtransform_sel<-qtransform_sel[-c(74,94,100),]

seasonal_evaluation<-analyse_all(NetCDFfolder,flood_values,qtransform_sel)
```



## Extracting correlations for all catchments for one month:

Returns: matrix with correlation in first colomn, and minimum significant correlation in the second.

out is the output from 'analyse_all'
ndata is the number of data used for evaluation. Necessary for calculating significance of correlation
mm is the month

```R
corr_march<-get_corr_month(out=seasonal_evaluation,ndata=qtransform_sel[,2],mm=3)
```



## Extracting correlations for all catchments and all months:

Returns: matrix with correlation in first colomn, and minimum significant correlation in the second.
out is the output from 'analyse_all'
ndata is the number of data used for evaluation. Necessary for calculating significance of correlation

```R
corr_all<-get_corr_all(out=seasonal_evaluation,ndata=qtransform_sel[,2])
```


## Extracting crpss for all catchments for one month:

Returns: matrix with crpss in first colomn, and s.e. in the second.

out is the output from 'analyse_all'
ndata is the number of data used for evaluation. Necessary for calculating significance of correlation
mm is the month

```R
crpss_march<-get_crpss_month(out=seasonal_evaluation,mm=3)
```


## Extracting crpss for all catchments for all month:

Returns: matrix with crpss for months Jan-July, max crpss for all months, month for max crpss and sigma-s for crpss

```R
crpss_all<-get_crpss_all(out=seasonal_evaluation)
```




## Extracting roc_area for all catchments for one month:

Returns: matrix with roc area in first colomn, and p-value in the second.

out is the output from 'analyse_all'
mm is the month
ri is the index for the thershold used for evaluating forecasts. The thresholds are given in 'flood_values'

```R
roc_march<-get_roc_month(out=seasonal_evaluation,mm=3,ri=1)
```




## Extracting roc_area for all catchments for all months:

Returns: matrix with roc area in first colomn, and p-value in the second.

out is the output from 'analyse_all'
ri is the index for the thershold used for evaluating forecasts. The thresholds are given in 'flood_values'

```R
roc_all<-get_roc_all(out=seasonal_evaluation,ri=1)
```





## Extracting Brier skill score for all catchments for one month:

Returns: matrix with Brier skill score  in first colomn, and p-value in the second.

out is the output from 'analyse_all'
mm is the month

```R
bss_march<-get_bss_month(out=seasonal_evaluation,mm=3)
```




## Extracting Brier skill score  for all catchments for all months:

Returns: matrix with Brier skill score  in first colomn, and p-value in the second.

out is the output from 'analyse_all'

```R
bss_all<-get_bss_all(out=seasonal_evaluation)
```

