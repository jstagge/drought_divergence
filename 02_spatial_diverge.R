# *------------------------------------------------------------------
# | PROGRAM NAME: spatial_diverge
# | FILE NAME: 02_spatial_diverge.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This script calculates trends in drought occurence
# |				using a linear model at the grid scale.
# |				Reproduces Figure 2 from "Observed Drought Indices ..."
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
 
require(ncdf4)
require(ggplot2)
require(scales)
require(grid)
require(gridExtra)
require(zoo)
require(RColorBrewer)
require(reshape2)
require(Rmpi)
require(snow)
require(snowfall)

#source("CustomThemes.R")

dev.new()

###########################################################################
## Set the Paths
###########################################################################
	
### path for Land Mask
WFDPath <- "~/nobackup/WFD/WATCH_Forcing_Data/"
WFDEIPath <- "~/nobackup/WFDEI/WATCH_Forcing_Data/"

landPath <- file.path(WFDPath,"WFD-land-lat-long-z-corrected.nc")
borderPath <- file.path(WFDPath,"Borders")

resultsPath <- paste(WFDPath,"Clustering/SPI_6_Jan", sep="")

combinedPath <- file.path(resultsPath,"Combined6")
SPIPath <- file.path(resultsPath,"SPI6")
SPEIPath <- file.path(resultsPath,"SPEI6")

###########################################################################
## Get national borders for plotting
###########################################################################
border <- read.table(file.path(borderPath,"Borders_MWDB3.cno"),sep=",",na.strings=9999,fill=T)
around <- which(abs(diff(border[,1]))>180)+1
border[around,] <- NA
colnames(border) <- c("x","y")


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

ggplot(LandEdit.frame, aes(x=Lon, y=Lat, fill=Land)) + geom_tile() + theme_bw()

SPI.LandID.list <- c(paste("SPI_",LandEdit.frame$Land, sep=""))
SPEI.LandID.list <- c(paste("SPEI_",LandEdit.frame$Land, sep=""))
LandID.list <- c(SPI.LandID.list, SPEI.LandID.list)

SPI.land <- data.frame(LandEdit.frame,Subset="SPI")
SPEI.land <- data.frame(LandEdit.frame,Subset="SPEI")


###########################################################################
## Read in SPI
###########################################################################

trendcoefdec.df <- read.csv("trend_coef_decade.csv")
percperdec.df <- read.csv("perc_per_decade.csv")
trendsig.df <- read.csv("trend_sig.csv")
pval.df <- read.csv("trend_pval.csv")

trendcoefdec.df <- data.frame(LandEdit.frame, trendcoefdec.df)
percperdec.df <- data.frame(LandEdit.frame, percperdec.df)
trendsig.df <- data.frame(LandEdit.frame, trendsig.df)
pval.df <- data.frame(LandEdit.frame, pval.df)

############################
#### Plot Results Map
############################
plot.border <- 2
LatLim <- c(33,71)
LonLim <- c(-13,48)

require(directlabels)


p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=SPEI))
p <- p + scale_fill_gradient2(low="blue", high="red")
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#p <- p + geom_path(data=border, aes(x=x, y=y), colour="black", alpha=0.25, size=.3)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + coord_cartesian(xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))
p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2.5,2.5), seq(2.5, 20, 2.5)), colour="grey50")
p <- p + geom_contour(aes(z=SPEI), breaks=c(-20,-10,10,20), colour="black")
p <- p + theme_bw()
p

direct.label(p)

#p <- p  + coord_map(projection = "albers", mean(LandEdit.frame$Lon), mean(LandEdit.frame$Lat), xlim=c(min(LandEdit.frame$Lon)-plot.border,max(LandEdit.frame$Lon)+plot.border), ylim=c(min(LandEdit.frame$Lat)-plot.border,max(LandEdit.frame$Lat)+plot.border))

#dev.new(width=3.5, height=3) 
plot(p)
#dev.off()




p <- ggplot(pval.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=SPEI))
p <- p + scale_fill_gradient(low="red", high="white", limits=c(0,0.1))
#p <- p + geom_path(data=border, aes(x=x, y=y), colour="black", alpha=0.25, size=.3)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + coord_cartesian(xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))
#p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2.5,2.5), seq(2.5, 20, 2.5)), colour="grey50")
#p <- p + geom_contour(aes(z=SPEI), breaks=c(-20,-10,10,20), colour="black")
p <- p + theme_bw()
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p







p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=SPEI))
p <- p + scale_fill_gradient2(low="blue", high="red")
#p <- p + geom_path(data=border, aes(x=x, y=y), colour="black", alpha=0.25, size=.3)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + coord_cartesian(xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))
p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2.5,2.5), seq(2.5, 20, 2.5)), colour="grey50")
p <- p + geom_contour(aes(z=SPEI), breaks=c(-20,-10,10,20), colour="black")
p <- p + geom_point(data=subset(pval.df, SPEI < 0.05), aes(x=Lon, y=Lat), size=0.1)
p <- p + theme_bw()
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p



p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=Diff))
p <- p + scale_fill_gradient2(low="blue", high="red", limits=c(-4,4))
#p <- p + geom_path(data=border, aes(x=x, y=y), colour="black", alpha=0.25, size=.3)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + coord_cartesian(xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))
p <- p + geom_contour(aes(z=Diff), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey50")
#p <- p + geom_contour(aes(z=Diff), breaks=c(-20,-10,10,20), colour="black")
#p <- p + geom_point(data=subset(pval.df, Diff < 0.1), aes(x=Lon, y=Lat), size=0.1)
p <- p + theme_bw()
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p



percperdec.df$TrueDiff <- percperdec.df$SPEI - percperdec.df$SPI
limit.list <- c(-abs(max(percperdec.df$TrueDiff, na.rm=TRUE)), abs(max(percperdec.df$TrueDiff, na.rm=TRUE)))
break.list <- pretty(limit.list, n = 9, min.n = 7)
 
 x.lim.list <- c(LonLim[1]-plot.border,LonLim[2]+plot.border)
 y.lim.list <- c(LatLim[1]-plot.border,LatLim[2]+plot.border)
p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=TrueDiff))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
#p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_gradientn(name="diff blah", colours = col.scheme, limits=limit.list, breaks=break.list, na.value = "grey70", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_gradient2(name="Difference in Median SPI3", low = "red", mid = "white", high = "blue", na.value = "grey80", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5), breaks=break.list)
#p <- p + scale_fill_gradient2(name="Median SPI3 Difference", low = "#B2182B", mid = "white", high = "#2166AC", na.value = "grey80", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5), breaks=break.list)
#p <- p + scale_fill_gradientn(colours=rainbow(5), limits = c(0,0.2))
p <- p + geom_contour(aes(z=TrueDiff), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p



p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
#p <- p + coord_cartesian(xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))

#p <- p + geom_contour(aes(z=Diff), breaks=c(-20,-10,10,20), colour="black")
#p <- p + geom_point(data=subset(pval.df, Diff < 0.1), aes(x=Lon, y=Lat), size=0.1)
p <- p + theme_bw()
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p



limit.list <- c(-abs(max(percperdec.df$SPI, na.rm=TRUE)), abs(max(percperdec.df$SPI, na.rm=TRUE)))
break.list <- pretty(limit.list, n = 9, min.n = 7)
 
 
p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=SPI))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="black", alpha=0.25, size=.3)
#p <- p  + coord_map(projection = "albers", 10, 45, xlim=c(-25,50), ylim=c(28,72))
p <- p + coord_fixed(ratio = 1, xlim=c(LonLim[1]-plot.border,LonLim[2]+plot.border), ylim=c(LatLim[1]-plot.border,LatLim[2]+plot.border))
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_gradientn(name="diff blah", colours = col.scheme, limits=limit.list, breaks=break.list, na.value = "grey70", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_gradient2(name="Difference in Median SPI3", low = "red", mid = "white", high = "blue", na.value = "grey80", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5), breaks=break.list)
#p <- p + scale_fill_gradient2(name="Median SPI3 Difference", low = "#B2182B", mid = "white", high = "#2166AC", na.value = "grey80", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5), breaks=break.list)
#p <- p + scale_fill_gradientn(colours=rainbow(5), limits = c(0,0.2))
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p







require(scales)
require(gridExtra)
source("CustomThemes.R")

##########################################
#### Plotting Scales
##########################################
RdBuPalette = colorRampPalette(brewer.pal(9, "RdBu"))
SpecPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplotPalette <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

col.scheme <- rev(RdBuPalette(100))



################################
### Plot SPI
##############################
limit.list <- c(-abs(max(percperdec.df$SPI, na.rm=TRUE)), abs(max(percperdec.df$SPI, na.rm=TRUE)))
limit.list <- c(-9,9)
break.list <- pretty(limit.list, n =9, min.n = 7)

break.list <- c(-12, break.list, 12)
 limit.list <- c(-abs(max(percperdec.df$SPI, na.rm=TRUE)), abs(max(percperdec.df$SPI, na.rm=TRUE)))

pval.plot <- pval.df 
pval.plot$SPIsig <- ""
pval.plot$SPIsig[pval.plot$SPI < 0.1] <-  "< 10%"
pval.plot$SPIsig[pval.plot$SPI < 0.05] <-  "< 5%"
pval.plot$SPIsig <- factor(pval.plot$SPIsig)
 
 
 x.lim.list <- c(LonLim[1]-plot.border,LonLim[2]+plot.border)
 y.lim.list <- c(LatLim[1]-plot.border,LatLim[2]+plot.border)
p <- ggplot(percperdec.df, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=SPI))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
#p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_gradientn(name="SPI6 Drought Likelihood Trend   (%/Decade)", colours = col.scheme, limits=limit.list, breaks=break.list, na.value = "grey70", guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + geom_contour(aes(z=SPI), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
#p <- p + geom_point(data=subset(pval.df, SPI < 0.05), aes(x=Lon, y=Lat, ), size=0.2, shape=47)
p <- p + geom_point(data=pval.plot, aes(x=Lon, y=Lat, shape=SPIsig))
p <- p + scale_shape_manual(name="Trend Significance", values=c(NA, 92, 4), guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p

ggsave("SPI_spatial_signif.png", p, width=5.5, height=5, dpi=600)
ggsave("SPI_spatial_signif.pdf", p, width=5.5, height=5)


p <- ggplot(percperdec.df, aes(x=SPI))
p <- p + geom_bar(fill="grey30", aes(y = ..density..), binwidth=1)
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("SPI6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", labels=percent)
p <- p + coord_cartesian(xlim=c(-10,13), ylim=c(0,0.19))
p
ggsave("SPI_spatial_hist.png", p, width=5.5, height=3, dpi=600)
ggsave("SPI_spatial_hist.pdf", p, width=5.5, height=3)


p <- ggplot(percperdec.df, aes(x=SPEI))
p <- p + geom_bar(fill="grey30", aes(y = ..density..), binwidth=1)
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("SPEI6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", labels=percent)
p <- p + coord_cartesian(xlim=c(-10,13), ylim=c(0,0.19))
p
ggsave("SPEI_spatial_hist.png", p, width=5.5, height=3, dpi=600)
ggsave("SPEI_spatial_hist.pdf", p, width=5.5, height=3)




p <- ggplot(percperdec.df, aes(x=SPI))
p <- p + geom_density(fill="grey30")
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("SPI6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", labels=percent)
p

ggsave("SPI_spatial_densiy.png", p, width=5.5, height=3, dpi=600)
ggsave("SPI_spatial_density.pdf", p, width=5.5, height=3)


###################################
###  Plot SPI with breaks
##################################
r <- percperdec.df
break.list <- seq(-12,12,3)
r$value.bin <- cut(r$SPI, break.list)
RdBuPalette = colorRampPalette(brewer.pal(length(levels(r$value.bin)), "RdBu"))
 
p <- ggplot(r, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=value.bin))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
#p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_manual(name="SPI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(length(levels(r$value.bin)))), drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(14))[c(1:6,9:14)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(10))[c(1:4,7:10)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + geom_contour(aes(z=SPI), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
p <- p + geom_point(data=subset(pval.df, SPI < 0.05), aes(x=Lon, y=Lat), size=0.1, shape=3)
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p
ggsave("SPI_spatial_v1.png", p, width=5.5, height=5, dpi=600)
ggsave("SPI_spatial_v1.pdf", p, width=5.5, height=5)




###################################
###  Plot SPEI with breaks
##################################
r <- percperdec.df
break.list <- seq(-12,12,3)
r$value.bin <- cut(r$SPEI, break.list)
RdBuPalette = colorRampPalette(brewer.pal(length(levels(r$value.bin)), "RdBu"))
 
p <- ggplot(r, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=value.bin))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
#p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(length(levels(r$value.bin)))), drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(14))[c(1:6,9:14)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(10))[c(1:4,7:10)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
p <- p + geom_point(data=subset(pval.df, SPEI < 0.05), aes(x=Lon, y=Lat), size=0.1, shape=3)
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p
ggsave("SPEI_spatial_v1.png", p, width=5.5, height=5, dpi=600)
ggsave("SPEI_spatial_v1.pdf", p, width=5.5, height=5)






#####################################
###  Run for Diff
#######################################
## Calculate difference as difference in trends
percperdec.df$TrueDiff <- percperdec.df$SPEI - percperdec.df$SPI


p <- ggplot(percperdec.df, aes(x=TrueDiff))
p <- p + geom_bar(fill="grey30", aes(y = ..density..), binwidth=1)
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("Diff6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", breaks=seq(0,1,0.05), labels=percent)
p <- p + coord_cartesian(xlim=c(-10,13), ylim=c(0,0.34))
p
ggsave("Diff_spatial_hist.png", p, width=5.5, height=3, dpi=600)
ggsave("Diff_spatial_hist.pdf", p, width=5.5, height=3)



p <- ggplot(percperdec.df, aes(x=SPI))
p <- p + geom_bar(fill="grey30", aes(y = ..density..), binwidth=1)
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("SPI6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", breaks=seq(0,1,0.05), labels=percent)
p <- p + coord_cartesian(xlim=c(-10,13), ylim=c(0,0.34))
p
ggsave("SPI_spatial_hist_v2.png", p, width=5.5, height=3, dpi=600)
ggsave("SPI_spatial_hist.pdf", p, width=5.5, height=3)


p <- ggplot(percperdec.df, aes(x=SPEI))
p <- p + geom_bar(fill="grey30", aes(y = ..density..), binwidth=1)
p <- p + geom_vline(y=0, colour="red", size=0.7)
p <- p + theme_pub_majgrid(9)
p <- p  +  scale_x_continuous("SPEI6 Drought Likelihood Trend   (%/Decade)") + scale_y_continuous("Proportion", breaks=seq(0,1,0.05), labels=percent)
p <- p + coord_cartesian(xlim=c(-10,13), ylim=c(0,0.34))
p
ggsave("SPEI_spatial_hist.png", p, width=5.5, height=3, dpi=600)
ggsave("SPEI_spatial_hist.pdf", p, width=5.5, height=3)










###################################
###  Plot Diff with breaks
##################################
r <- percperdec.df
break.list <- seq(-12,12,3)
r$value.bin <- cut(r$TrueDiff, break.list)
RdBuPalette = colorRampPalette(brewer.pal(length(levels(r$value.bin)), "RdBu"))
 
p <- ggplot(r, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=value.bin))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
#p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_manual(name="Diff6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(length(levels(r$value.bin)))), drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(14))[c(1:6,9:14)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(10))[c(1:4,7:10)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
#p <- p + geom_point(data=subset(pval.df, Diff < 0.05), aes(x=Lon, y=Lat), size=0.1, shape=3)
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p
ggsave("Diff_spatial_v1.png", p, width=5.5, height=5, dpi=600)
ggsave("Diff_spatial_v1.pdf", p, width=5.5, height=5)


r <- percperdec.df
break.list <- c(-12,seq(-6,6,2),12)
r$value.bin <- cut(r$TrueDiff, break.list)
RdBuPalette = colorRampPalette(brewer.pal(length(levels(r$value.bin)), "RdBu"))
 
p <- ggplot(r, aes(x=Lon,y=Lat))
p <- p + geom_tile(aes(fill=value.bin))
p <- p + geom_path(data=border, aes(x=x, y=y), colour="grey50", alpha=0.25, size=.3)
p <- p  + coord_map(projection = "albers",  mean( x.lim.list),  mean( y.lim.list), xlim=x.lim.list, ylim=y.lim.list)
#p <- p + coord_fixed(ratio = 1, xlim= x.lim.list, ylim= y.lim.list)
p <- p  +  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p <- p + theme_pub_map_grey()
p <- p + scale_fill_manual(name="Diff6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(length(levels(r$value.bin)))), drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(14))[c(1:6,9:14)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + scale_fill_manual(name="SPEI6 Drought Likelihood Trend   (%/Decade)", values = rev(RdBuPalette(10))[c(1:4,7:10)], drop=FALSE, guide=guide_legend(label.position = "bottom", title.position="top", label.hjust = 0.5))
#p <- p + geom_contour(aes(z=SPEI), breaks=c(seq(-20,-2,2), seq(2, 20, 2)), colour="grey10")
#p <- p + geom_point(data=subset(pval.df, Diff < 0.05), aes(x=Lon, y=Lat), size=0.1, shape=3)
p <- p + theme(legend.position="bottom", legend.direction="horizontal")
p <- p + theme(legend.key.size = unit(0.9, "lines"))
p <- p + theme(legend.title = element_text(face = "bold", hjust = 0))
p <- p + theme(legend.margin = unit(0, "line"))
p <- p + theme(legend.key.width = unit(1.7, "lines"))
p <- p + theme(plot.margin = unit(c(0.5,1,0.5,0.2), "lines"))
p <- p + theme(panel.background =  element_rect(fill = "white", colour=NA))
p
ggsave("Diff_spatial_v2.png", p, width=5.5, height=5, dpi=600)
ggsave("Diff_spatial_v2.pdf", p, width=5.5, height=5)



