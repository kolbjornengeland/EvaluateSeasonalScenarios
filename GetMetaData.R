# installering av pakker
install.packages("devtools")
install.packages("curl")
library(devtools)
install_github("NVE/NVEDATA")

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



