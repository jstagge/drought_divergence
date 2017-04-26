# *------------------------------------------------------------------
# | PROGRAM NAME: prepare_file_system
# | FILE NAME: 00_prepare_file_system.R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  This script will install any packages needed for the subsequent
# |				code and will create the output folder structure to hold results.     
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
###  Check for necessary packages and install if needed
###########################################################################
### Set a list of packages
list_of_packages <- c("ggplot2", "zoo", "mgcv")


### Determine which packages are missing
package_list <- installed.packages()[,"Package"]
installed_test <- (list_of_packages %in% package_list)
packages_needed <- list_of_packages[!installed_test]

### If packages are missing, install them
if(length(packages_needed)) install.packages(packages_needed)



###########################################################################
## Set the Paths
###########################################################################
### Path for Data and Output	
data_path <- "./data"
output_path <- "./output"
function_path <- "./functions"

###########################################################################
###  Generate folder structure
###########################################################################
### Create output folder
dir.create(output_path)

### Create subfolders
dir.create(file.path(output_path, "figures"))
dir.create(file.path(output_path, "continental_diverge"))
dir.create(file.path(output_path, "spatial_diverge"))
dir.create(file.path(output_path, "climate_trends"))
dir.create(file.path(output_path, "et_comparison"))



