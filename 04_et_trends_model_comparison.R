# *------------------------------------------------------------------
# | PROGRAM NAME: et_trends_model_comparison
# | FILE NAME: 04_et_trends_model_comparison.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This script calculates trends in reference evapotranspiration
# |				using different ET0 equations.
# |				Reproduces Figure 4 from "Observed Drought Indices ..."
# |
# *------------------------------------------------------------------
# | COMMENTS:               
# |
# |  1:  
# |  2:  
# |  3: 
# |*------------------------------------------------------------------
# | DATA USED:               
# | 
# |
# |*------------------------------------------------------------------
# | CONTENTS:               
# |
# |  PART 1:  
# |  PART 2: 
# |  PART 3: 
# *-----------------------------------------------------------------
# | UPDATES:               
# |
# |
# *------------------------------------------------------------------
 


### This function clears the memory in R, in case there are old variables from previous work.
rm(list=ls())


## packages
require(ncdf4)
require(ggplot2)
require(zoo)
#require(raster)
require(Rmpi)
require(snow)


###########################################################################
## Introduce Functions
###########################################################################
source("NCDFFunctions.R")


accum.func <- function(x, period, process="mean") {
  time.scale <- period
  x <- as.numeric(x)
  nn <- length(x)
  time.index <- 1:nn
  
    if(time.scale>1){
    ## compute running mean
    ffrom <- sapply((1:nn) - time.scale + 1, function(x) max(x,1))
    tto <- 1:nn
    elements <- apply(cbind(ffrom,tto), 1,function(x) seq(x[1], x[2]) )
    if(is.matrix(elements))
      elements <- as.data.frame(elements)
	if (process=="mean") {
        funct <- function(which,what) mean(what[which],na.rm=FALSE)  
    } else if (process=="sum") {
		funct <- function(which,what) sum(what[which],na.rm=FALSE) 
	}	
	x <- sapply(elements, funct, what=x)
    x[1:(time.scale-1)] <- NA
  }
  return(x)
}

season.func <- function(x, ref.test, method="mean") {
nn <- length(x)
### Create daily index
 mday <- rep(1:365,length.out=nn)

 seas.ts <- rep(NA,365)
 
 for(mm in 1:365){
	xm <- x[mday==mm] ## select month
    xm.fit <- x[mday==mm&ref.test] ## select data in the reference time interval

	if (method=="mean"){
	 seas.ts[mm] <- mean(xm.fit)
	} else if (method=="median") {
	  seas.ts[mm] <- median(xm.fit)
	}
	
  x[mday==mm] <- xm - seas.ts[mm]
}
results.list <- list(ts=x, seas.ts=seas.ts)
return(results.list)
}



season.func2 <- function(x, ref.test, method="mean") {
nn <- length(x)
### Create daily index
 mday <- rep(1:365,length.out=nn)

 seas.ts <- rep(NA,365)
 
 for(mm in 1:365){
	xm <- x[mday==mm] ## select month
    xm.fit <- x[mday==mm&ref.test] ## select data in the reference time interval

	if (method=="mean"){
	 seas.ts[mm] <- mean(xm.fit)
	} else if (method=="median") {
	  seas.ts[mm] <- median(xm.fit)
	} else if (method=="perc20") {
	  seas.ts[mm] <- quantile(xm.fit,0.2)[[1]]
	} else if (method=="perc80") {
	  seas.ts[mm] <- quantile(xm.fit,0.8)[[1]]
	}
}
return(seas.ts)
}

read.func <- function(pe.folder, pe.file, land.list, date.list) {

PEPath <- file.path(mainPath,paste("PET_Calc/",pe.folder,sep=""))

i <- 1
leap.year <- 0
PEFile <- file.path(PEPath,paste(pe.file,date.list[i],".nc", sep=""))

## Read in files
PE.nc <- ReadNCDF(PEFile, print=FALSE)
PE.mat <- PE.nc$var$PET$vals
PE.land <- PE.nc$dim$land$vals

PE.mat.Europe <- matrix(PE.mat[PE.land %in% land.list],length(land.list),dim(PE.mat)[2])

for(i in 2:length(date.list)) {
#for(i in 2:5) {

year <- substring(date.list[i],nchar(date.list[i])-8,nchar(date.list[i])-2)
month <- substring(date.list[i],nchar(date.list[i])-1,nchar(date.list[i])-0)

PEFile <- file.path(PEPath,paste(pe.file,date.list[i],".nc", sep=""))

## Read in files
PE.nc <- ReadNCDF(PEFile, print=FALSE)
new.PE <- PE.nc$var$PET$vals
PE.land <- PE.nc$dim$land$vals

if(dim(new.PE)[2]==29) {
	leap.year <- 1
}

if (month==12){
	if(leap.year==1) {
		new.PE <- new.PE[,1:30]
	} else {
		new.PE <- new.PE  }
	leap.year <- 0
}
  
new.PE.Europe <- matrix(new.PE[PE.land %in% land.list],length(land.list),dim(new.PE)[2])

 
#PE.mat <- cbind(PE.mat, new.PE)
PE.mat.Europe <- cbind(PE.mat.Europe, new.PE.Europe)

}
return(PE.mat.Europe)

}

###########################################################################
## Set the Paths
###########################################################################
### path for Land Mask
landPath <- "WFD-land-lat-long-z-corrected.nc"

systemPath <- "~/nobackup/WFD"
mainPath <- file.path(systemPath,"WATCH_Forcing_Data")


###########################################################################
## Get and prepare land mask
###########################################################################
land.nc <- nc_open(landPath)
print(land.nc)

## Extract the land (latitude longitude) from the file
myLand <- data.frame(
                     land=ncvar_get(land.nc,"land"),
                     Longitude=ncvar_get(land.nc,"Longitude"),
                     Latitude=ncvar_get(land.nc,"Latitude"),
                     Z=ncvar_get(land.nc,"Z")
                     )
nc_close(land.nc);rm(land.nc)
## Based on spatial extend of EWA data set
## LatLim <- sort(ceiling(c(37.066,70.073)))
## LonLim <- sort(ceiling(c(-9.535, 30.767)))

### Based on visual approximations
LatLim <- c(33,72)
LonLim <- c(-28,48)

LandInEurope <- with(myLand,{
    Latitude >= LatLim[1] & Latitude <= LatLim[2] &
    Longitude >= LonLim[1] & Longitude <= LonLim[2]
})
 plot(Latitude~Longitude,data=myLand[LandInEurope,],pch=19)

 
IcelandMask <- with(myLand,{
    Latitude >= 62.5 & Latitude <= LatLim[2] &
    Longitude >= LonLim[1] & Longitude <= -5
})
  plot(Latitude~Longitude,data=myLand[IcelandMask,],pch=19)

  
AzoresMask <- with(myLand,{
    Latitude >= LatLim[1] & Latitude <= 40 &
    Longitude >= LonLim[1] & Longitude <= -15
})
  plot(Latitude~Longitude,data=myLand[AzoresMask,],pch=19)

  
LandEdit <- LandInEurope & !IcelandMask & !AzoresMask
   plot(Latitude~Longitude,data=myLand[LandEdit,],pch=19)

 
#Lon.seq <- seq(min(myLand$Longitude[LandInEurope]), max(myLand$Longitude[LandInEurope]), 0.5)
#Lat.seq <- seq(min(myLand$Latitude[LandInEurope]), max(myLand$Latitude[LandInEurope]), 0.5)
#Lat.seq <- seq(LatLim[1]+.25,LatLim[2]-.25,0.5)
LandInEurope <- myLand$land[LandInEurope]
LandEdit <- myLand$land[LandEdit]

Land.frame <- data.frame(Land = myLand$land[myLand$land %in% LandInEurope], Lon = myLand$Longitude[myLand$land %in% LandInEurope], Lat = myLand$Latitude[myLand$land %in% LandInEurope])

LandEdit <- Land.frame$Land %in% LandEdit
LandEdit.frame <- Land.frame[LandEdit,]

#ggplot(LandEdit.frame, aes(x=Lon, y=Lat, fill=Land)) + geom_tile() + theme_bw()


###########################################################################
##  Prepare time sequence without leap days
###########################################################################
WFD.dates <- seq(as.Date("1958/1/1"), as.Date("2001/12/31"), "days")
WFDEI.dates <- seq(as.Date("1979/1/1"), as.Date("2012/12/31"), "days")

WFD.month.day <- substring(WFD.dates,6,10)
WFDEI.month.day <- substring(WFDEI.dates,6,10)

### Create time index with no leap years
WFD.noleap <- WFD.dates[WFD.month.day != "02-29"]
WFDEI.noleap <- WFDEI.dates[WFDEI.month.day != "02-29"]


################################################################
###  Cut a reference period
################################################################
ref.start <- 1970
ref.end <- 1999

x.zoo <- zoo(seq(1,length(WFD.noleap)),WFD.noleap)
htime <- time(x.zoo)
htime <- as.POSIXlt(htime)$year + 1900
hi <- 1:length(htime) 

 if(missing(ref.start)){
     cat("Missing Ref.start")
    ref.start <- min(hi)
  } else {
    if(ref.start>=min(htime)&ref.start<=max(htime)){
	ref.start <- min(hi[htime==ref.start])
    } else {
      stop("'ref.start' outside of observed time period")
    }
  }
if(missing(ref.end)){
    ref.end <- max(hi)
  } else {
    if(ref.end<=max(htime)&ref.end>=min(htime)){
      ref.end <- max(hi[htime==ref.end])
    } else {
      stop("'ref.start' outside of observed time period")
    }
  }
  
  
 if(ref.start>=ref.end)
  stop("'ref.start' should be < 'ref.end'")     
 
  nn <- length(WFD.noleap)
  time.index <- 1:nn
  fit.time <- time.index>=ref.start&time.index<=ref.end
  

###########################################################################
##  Read in WFD and WFDEI data
###########################################################################

year.list <- seq(1958,2001,1)
date.list <- as.numeric(paste(rep(year.list,each=12),sprintf("%02d", rep(1:12,length(year.list))), sep=""))

PE.list <- list()
PE.models <- c("Thornthwaite","Hargreaves","P-M No Rad","Priestly Taylor","P-M Rad")

## Thornthwaite
PE.list[[1]] <- read.func(pe.folder="PET_Thornth", pe.file="PE_Thornth_", land.list=LandEdit.frame$Land, date.list=date.list)

## Hargreaves
PE.list[[2]] <- read.func(pe.folder="PET_Hargreaves", pe.file="PE_Hargreaves_", land.list=LandEdit.frame$Land, date.list=date.list)

## P-M No Rad
PE.list[[3]] <- read.func(pe.folder="PET_Penman_NoRad", pe.file="PE_Penman_norad_", land.list=LandEdit.frame$Land, date.list=date.list)

## Priestly Taylor
PE.list[[4]] <- read.func(pe.folder="PET_Priest_Taylor", pe.file="PE_Priest_", land.list=LandEdit.frame$Land, date.list=date.list)

## P-M Rad
PE.list[[5]] <- read.func(pe.folder="PET_Penman_Rad", pe.file="PE_Penman_rad_", land.list=LandEdit.frame$Land, date.list=date.list)


## Test plots
test.df <- data.frame(LandEdit.frame, values=PE.list[[3]][,200])
ggplot(test.df, aes(x=Lon, y=Lat, fill=values)) + geom_tile() + theme_bw()+scale_fill_gradient(low="white", high="red")

plot(WFD.noleap,PE.list[[3]][3000,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))




###########################################################################
## Apply a 6 month accumulation
###########################################################################
#Run in parallel

## Start the MPI worker processes:
## Version for Abel
numWorkers <- as.numeric(Sys.getenv("SLURM_NTASKS")) - 1
cluster <- makeCluster(numWorkers, type = "MPI")

Sys.time()
for(k in 1:5){
PE.list[[k]] <- parApply(cluster, PE.list[[k]],1,accum.func,period=183, process="sum")
PE.list[[k]] <- t(PE.list[[k]])
}
Sys.time()


plot(WFD.noleap,PE.list[[3]][50,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))

###########################################################################
## Plot the Spatial Mean over Europe
###########################################################################
WFD.mean.ts <- list()

for(k in 1:5){
WFD.mean.ts[[k]] <- apply(PE.list[[k]], 2, mean, na.rm=TRUE)
if(k==1) {
	Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.ts[[k]], Data=PE.models[[k]])
} else {
	Temp.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.ts[[k]], Data=PE.models[[k]])
	Plot.df <- rbind(Plot.df, Temp.df)
}
}

write.csv(Plot.df, "PET_All_TimeSeriesPlot.csv",row.names=FALSE)


###########################################################################
## Option 1: Remove seasonality from European mean
###########################################################################

for(k in 1:5){
WFD.noseason <- season.func(WFD.mean.ts[[k]], fit.time, method="mean")
WFD.noseason <- WFD.noseason$ts
if(k==1) {
	Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.noseason, Data=PE.models[[k]])
} else {
	Temp.df <- data.frame(Date=WFD.noleap, Value=WFD.noseason, Data=PE.models[[k]])
	Plot.df <- rbind(Plot.df, Temp.df)
}
}

write.csv(Plot.df, "PET_All_Option1_df.csv",row.names=FALSE)


###########################################################################
## Option 2: Remove seasonality from each grid cell
###########################################################################

for(k in 1:5){
Ref.seas.mean <- t(apply(PE.list[[k]],1,season.func2, ref.test=fit.time, method="mean"))

WFD.seas.mean <- t(apply(Ref.seas.mean, 1, rep, length.out=dim(PE.list[[k]])[2]))

WFD.noseason.grid <- PE.list[[k]] - WFD.seas.mean
WFD.mean.grid <- apply(WFD.noseason.grid, 2, mean, na.rm=TRUE)

if(k==1) {
	Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.grid, Data=PE.models[[k]])
} else {
	Temp.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.grid, Data=PE.models[[k]])
	Plot.df <- rbind(Plot.df, Temp.df)
}
}


write.csv(Plot.df, "PET_All_Option2_df.csv",row.names=FALSE)


###########################################################################
## Option 4: Area above 20th percentile
###########################################################################


for(k in 1:5){

Ref.seas.perc <- t(apply(PE.list[[k]],1,season.func2, ref.test=fit.time, method="perc80"))

WFD.seas.perc <- t(apply(Ref.seas.perc, 1, rep, length.out=dim(PE.list[[k]])[2]))

WFD.perc.grid <- WFD.accum > WFD.seas.perc
WFD.perc.area.grid <- apply(WFD.perc.grid, 2, mean, na.rm=TRUE)

if(k==1) {
	Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.perc.area.grid, Data=PE.models[[k]])
} else {
	Temp.df <- data.frame(Date=WFD.noleap, Value=WFD.perc.area.grid, Data=PE.models[[k]])
	Plot.df <- rbind(Plot.df, Temp.df)
}
}


write.csv(Plot.df, "PET_All_Option4_df.csv",row.names=FALSE)


