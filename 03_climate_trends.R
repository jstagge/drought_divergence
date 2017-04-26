# *------------------------------------------------------------------
# | PROGRAM NAME: climate_trends
# | FILE NAME: 04_climate_trends.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This code calculates trends for each of the consitutent
# |				climate variables.
# |				Reproduces Figure 3 from "Observed Drought Indices ..."
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
 
### Clear any existing data or functions.
rm(list=ls())

###########################################################################
## Run code to process evapotranspiration
###########################################################################
source(file.path("code/ET_O", "01_et0_weighted.R"))
source(file.path("code/ET_O", "02_et0_trend.R"))


###########################################################################
## Run code to process Temperature
###########################################################################
source(file.path("code/Temp", "01_temp_weighted.R"))
source(file.path("code/Temp", "02_tair_trend.R"))
source(file.path("code/Temp", "03_tdiff_trend.R"))


###########################################################################
## Run code to process Wind Speed
###########################################################################
source(file.path("code/Wind", "01_wind_weighted.R"))
source(file.path("code/Wind", "02_wind_trend.R"))


