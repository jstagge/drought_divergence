
TwoDtoThreeD <- function(landID, TwoD.mat, Lat.seq, Lon.seq, orientation){

### orientation = 1 for time along y axis, 2 for time along x axis

### Create the index array for Watch Forcing Data
Lat.world <- seq(89.75,-89.75,-0.5)
Lon.world <- seq(-179.75,179.75,0.5)

totalcells <- length(Lon.world)*length(Lat.world)
indexarray <- matrix(seq(1,totalcells,1), nrow=length(Lat.world), ncol=length(Lon.world), byrow=TRUE)

### Subset the index array to the region described in Lat.seq and Lon.seq 
indexarray <- indexarray[which(Lat.world==max(Lat.seq)):which(Lat.world==min(Lat.seq)),which(Lon.world==min(Lon.seq)):which(Lon.world==max(Lon.seq))]

### Calculate time steps for matrix.  Catches error if only one time step.
if (is.vector(TwoD.mat)==TRUE) {
	t.length <- 1
} else {
	t.length <- dim(TwoD.mat)[orientation]
}

### Create empty 3-d matrix
ThreeD.mat <- array(NA,c(length(Lat.seq),length(Lon.seq),t.length))

### Loop through time steps, filling 3-D metrix
for(t in seq(1,t.length,1)) {
results.vec <- rep(NA,length(Lat.seq)*length(Lon.seq))

if (orientation==1) {
	data.for.vec <- TwoD.mat[t,]
	length(data.for.vec) <- length(landID)
	results.vec[indexarray %in% landID] <- data.for.vec
} else if (orientation==2) {
	data.for.vec <- TwoD.mat[,t]
	length(data.for.vec) <- length(landID)
	results.vec[indexarray %in% landID] <- data.for.vec
}

ThreeD.mat[,,t] <- matrix(results.vec, nrow=length(Lat.seq), ncol=length(Lon.seq))
}

return(ThreeD.mat)
}


ThreeDtoTwoD <- function(ThreeD.mat, orientation){

### orientation = 1 for time along y axis, 2 for time along x axis

TwoD.vector <- as.numeric(ThreeD.mat)

land.length <- dim(ThreeD.mat)[1]*dim(ThreeD.mat)[2]

mat_dim <- dim(ThreeD.mat)
if (length(mat_dim) == 2) {
time.length <- 1
} else {
time.length <- dim(ThreeD.mat)[3]
}

if (orientation==1) {
	TwoD.mat <- matrix(TwoD.vector, nrow=land.length, ncol=time.length)
} else if (orientation==2) {
	TwoD.mat <- matrix(TwoD.vector, nrow=time.length, ncol=land.length)
}

return(TwoD.mat)
}


ReadNCDF <- function(ncdf, print=TRUE){
### This function reads in ncdf files.  
### I am still in early stages of developing it, so it may fail for extremely large data sets.  
### It should be safe for monthly data

ncdf.file <- nc_open(ncdf)
if (print==TRUE) { print(ncdf.file) }

var.list <- list()
dim.list <- list()
att.list <- list()

for (i in 1:ncdf.file$nvars) {
var.name <- ncdf.file$var[[i]]$name
var.list[[var.name]] <- list(name=ncdf.file$var[[i]]$name, longname=ncdf.file$var[[i]]$longname, units=ncdf.file$var[[i]]$units, vals=ncvar_get(ncdf.file, var.name))
read_time <- format(Sys.time(), "%a %b %d %H:%M:%S %Y")
if (print==TRUE) { cat(read_time,"Variable ",var.name, "read in succesfully \n")}
}

for (i in 1:ncdf.file$ndims) {
dim.name <- ncdf.file$dim[[i]]$name
dim.list[[dim.name]] <- list(name=ncdf.file$dim[[i]]$name, units=ncdf.file$dim[[i]]$units, vals=ncdf.file$dim[[i]]$vals)
read_time <- format(Sys.time(), "%a %b %d %H:%M:%S %Y")
if (print==TRUE) { cat(read_time,"Dimension ",dim.name, "read in succesfully \n") }
}

atts <- ncatt_get(ncdf.file,0)
for (i in 1:ncdf.file$natts) {
read_time <- format(Sys.time(), "%a %b %d %H:%M:%S %Y")
if (print==TRUE) { cat(read_time,"Attribute ",names(atts)[i], "read in succesfully \n") }
}

full.nc <- list(var=var.list, dim=dim.list, att=atts)
nc_close(ncdf.file);rm(ncdf.file)

return(full.nc)
}
