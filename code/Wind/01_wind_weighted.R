### This function clears the memory in R, in case there are old variables from previous work.
rm(list=ls())


## packages
require(ncdf4)
#require(ggplot2)
require(zoo)
#require(raster)
require(Rmpi)
require(snow)
require(ggplot2)

require(SDMTools)
require(RcppRoll)
require(Hmisc)	

###########################################################################
## Introduce Functions
###########################################################################
source("NCDFFunctions.R")

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


###########################################################################
## Set the Paths
###########################################################################	

### path for Land Mask
WFDPath <- "~/nobackup/WFD/WATCH_Forcing_Data"
WFDEIPath <- "~/nobackup/WFDEI/WATCH_Forcing_Data"

landPath <- file.path(WFDPath,"WFD-land-lat-long-z-corrected.nc")
borderPath <- file.path(WFDPath,"Borders")

### path for WFD and WFDEI
WFD.name <- file.path(WFDPath,"Daily_Climate/WFD_wind_daily_2d_wide.nc")
WFDEI.wind.name <- file.path(WFDEIPath,"Daily_Climate/New/WFDEI_wind_daily_2d_wide.nc")

resultsPath <- file.path(WFDPath,"Clustering/SPI_6_Jan")

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

grid.data <- grid.info(lats=Land.frame$Lat, 0.5)
Land.frame <- data.frame(Land.frame, grid.data)

LandEdit <- Land.frame$Land %in% LandEdit
LandEdit.frame <- Land.frame[LandEdit,]

#ggplot(LandEdit.frame, aes(x=Lon, y=Lat, fill=Land)) + geom_tile() + theme_bw()
ggplot(LandEdit.frame, aes(x=Lon, y=Lat, fill=area)) + geom_tile() + theme_bw()

###########################################################################
##  Prepare time sequence without leap days
###########################################################################
WFD.dates <- seq(as.Date("1958/1/1"), as.Date("2001/12/31"), "days")
WFDEI.dates <- seq(as.Date("1979/1/1"), as.Date("2014/12/31"), "days")

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
## Read in WFD and WFDEI data
###########################################################################
# Read in WFD
WFD_Wind.nc <- nc_open(WFD.name)
WFD_Wind.nc
WFD_Wind.mat <- ncvar_get(WFD_Wind.nc, "Wind")
WFD.landid <- ncvar_get(WFD_Wind.nc, "land")
nc_close(WFD_Wind.nc)
rm(WFD_Wind.nc)

### Extract WFD to Wide size
WFD_Wind.mat.Europe <- matrix(WFD_Wind.mat[WFD.landid %in% LandInEurope],length(LandInEurope),dim(WFD_Wind.mat)[2])
#### Extract the bad parts
WFD_Wind.mat.Edit <- WFD_Wind.mat.Europe[LandEdit,]

rm(WFD_Wind.mat)
rm(WFD_Wind.mat.Europe)

# Read in WFDEI Wind
WFDEI_Wind.nc <- nc_open(WFDEI.wind.name)
WFDEI_Wind.nc
WFDEI_Wind.mat <- ncvar_get(WFDEI_Wind.nc, "Wind")
WFDEI.landid <- ncvar_get(WFDEI_Wind.nc, "landid")
nc_close(WFDEI_Wind.nc)
rm(WFDEI_Wind.nc)
 
 ### Extract WFDEIEI to Wide size
WFDEI_Wind.mat.Europe <- matrix(WFDEI_Wind.mat[WFDEI.landid %in% LandInEurope],length(LandInEurope),dim(WFDEI_Wind.mat)[2])
#### Extract the bad parts
WFDEI_Wind.mat.Edit <- WFDEI_Wind.mat.Europe[LandEdit,]
rm(WFDEI_Wind.mat)
rm(WFDEI_Wind.mat.Europe)


## Test plot
plot(WFD.noleap,WFD_Wind.mat.Edit[3000,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))
lines(WFDEI.noleap,WFDEI_Wind.mat.Edit[3000,], col="red")


###########################################################################
## Run everything for Wind
###########################################################################

## Start the MPI worker processes:
## Version for Abel
#numWorkers <- as.numeric(Sys.getenv("SLURM_NTASKS")) - 1
#cluster <- makeCluster(numWorkers, type = "MPI")

Sys.time()

WFD.accum <- t(apply(WFD_Wind.mat.Edit,1,roll_mean, n=183, fill=NA, align="right"))
WFDEI.accum <- t(apply(WFDEI_Wind.mat.Edit,1,roll_mean, n=183, fill=NA, align="right"))

#WFD.accum <- t(parApply(cluster, WFD_Wind.mat.Edit,1,roll_mean, n=183, fill=NA, align="right"))
#WFDEI.accum <- t(parApply(cluster, WFDEI_Wind.mat.Edit,1,roll_mean, n=183, fill=NA, align="right"))
Sys.time()

plot(WFD.noleap,WFD.accum[50,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))
lines(WFDEI.noleap,WFDEI.accum[50,], col="red")



###########################################################################
## Plot the Spatial Mean over Europe
###########################################################################
WFD.mean.ts <- apply(WFD.accum, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)
WFDEI.mean.ts <- apply(WFDEI.accum, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)

Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.ts, Data="WFD")
Plot.df <- rbind(Plot.df, data.frame(Date=WFDEI.noleap, Value=WFDEI.mean.ts, Data="WFDEI"))

#write.csv(Plot.df, "Wind_TimeSeriesPlot_Europe.csv",row.names=FALSE)


###########################################################################
## Option 1: Remove seasonality from European mean
###########################################################################
 WFD.noseason <- season.func(WFD.mean.ts, fit.time, method="mean")
 WFDEI.noseason <-  WFDEI.mean.ts - rep(WFD.noseason$seas.ts, length.out=length(WFDEI.mean.ts))
 WFD.noseason <- WFD.noseason$ts
 
Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.noseason, Data="WFD")
Plot.df <- rbind(Plot.df, data.frame(Date=WFDEI.noleap, Value=WFDEI.noseason, Data="WFDEI"))

write.csv(Plot.df, "Wind_Option1_Europe.csv",row.names=FALSE)


###########################################################################
## Option 2: Remove seasonality from each grid cell
###########################################################################

Ref.seas.mean <- t(apply(WFD.accum,1,season.func2, ref.test=fit.time, method="mean"))

WFD.seas.mean <- t(apply(Ref.seas.mean, 1, rep, length.out=dim(WFD.accum)[2]))
WFDEI.seas.mean <- t(apply(Ref.seas.mean, 1, rep, length.out=dim(WFDEI.accum)[2]))

WFD.noseason.grid <- WFD.accum - WFD.seas.mean
WFDEI.noseason.grid <- WFDEI.accum - WFDEI.seas.mean

# Test plot
plot(WFD.noleap, WFD.noseason.grid[50,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))
lines(WFDEI.noleap, WFDEI.noseason.grid[50,], col="red")

WFD.mean.grid <- apply(WFD.noseason.grid, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)
WFDEI.mean.grid <- apply(WFDEI.noseason.grid, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)

Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.mean.grid, Data="WFD")
Plot.df <- rbind(Plot.df, data.frame(Date=WFDEI.noleap, Value=WFDEI.mean.grid, Data="WFDEI"))

write.csv(Plot.df, "Wind_Option2_Europe.csv",row.names=FALSE)


###########################################################################
## Option 3: Use above, but plot percentiles
###########################################################################
prob.list <- c(0.01,0.05, 0.1,0.25,0.5,0.75,0.9,0.95,0.99)

WFD.full.test <- apply(WFD.noseason.grid, 2, function(x) sum(is.na(x)) != length(x))
WFD.quant.grid <- matrix(NA, length(prob.list), dim(WFD.noseason.grid)[2])
WFD.quant.grid[,WFD.full.test] <- apply(WFD.noseason.grid[,WFD.full.test], 2, wtd.quantile, weights=LandEdit.frame$area, probs=prob.list, normwt=TRUE, na.rm=TRUE) 
WFD.quant.grid <- t(WFD.quant.grid)
	
WFDEI.full.test <- apply(WFDEI.noseason.grid, 2, function(x) sum(is.na(x)) != length(x))
WFDEI.quant.grid <- matrix(NA, length(prob.list), dim(WFDEI.noseason.grid)[2])
WFDEI.quant.grid[,WFDEI.full.test] <- apply(WFDEI.noseason.grid[,WFDEI.full.test], 2, wtd.quantile, weights=LandEdit.frame$area, probs=prob.list, normwt=TRUE, na.rm=TRUE) 
WFDEI.quant.grid <- t(WFDEI.quant.grid)
	
percentile.name <- paste(prob.list*100, "%", sep="")
colnames(WFD.quant.grid) <- percentile.name
colnames(WFDEI.quant.grid) <- percentile.name

rm(Quant.df)
for (k in c(1,2,4,5,6,8,9)) {

column.num <- k
percentile.name <- colnames(WFD.quant.grid)

temp.wfd <- data.frame(Date=WFD.noleap, Value=WFD.quant.grid[,column.num], Data="WFD", Percentile=percentile.name[column.num])
temp.wfdei <- data.frame(Date=WFDEI.noleap, Value=WFDEI.quant.grid[,column.num], Data="WFDEI", Percentile=percentile.name[column.num])

if ( exists("Quant.df") == TRUE) {
	Quant.df <- rbind(Quant.df, temp.wfd, temp.wfdei)
} else {
	Quant.df <- rbind(temp.wfd, temp.wfdei)
}
}

write.csv(Quant.df, "Wind_Option3_Europe.csv",row.names=FALSE)


###########################################################################
## Option 4: Area above 20th percentile
###########################################################################
Ref.seas.perc <- t(apply(WFD.accum,1,season.func2, ref.test=fit.time, method="perc80"))

WFD.seas.perc <- t(apply(Ref.seas.perc, 1, rep, length.out=dim(WFD.accum)[2]))
WFDEI.seas.perc <- t(apply(Ref.seas.perc, 1, rep, length.out=dim(WFDEI.accum)[2]))

WFD.perc.grid <- WFD.accum > WFD.seas.perc
WFDEI.perc.grid <- WFDEI.accum > WFDEI.seas.perc

# Test plot
plot(WFD.noleap, WFD.perc.grid[50,], type="l", xlim=c(as.Date("1955-01-01"), as.Date("2015-01-01")))
lines(WFDEI.noleap, WFDEI.perc.grid[50,], col="red")

WFD.perc.area.grid <- apply(WFD.perc.grid, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)
WFDEI.perc.area.grid <- apply(WFDEI.perc.grid, 2, weighted.mean, LandEdit.frame$area, na.rm=TRUE)

Plot.df <- data.frame(Date=WFD.noleap, Value=WFD.perc.area.grid, Data="WFD")
Plot.df <- rbind(Plot.df, data.frame(Date=WFDEI.noleap, Value=WFDEI.perc.area.grid, Data="WFDEI"))


write.csv(Plot.df, "Wind_Option4_Europe.csv",row.names=FALSE)








