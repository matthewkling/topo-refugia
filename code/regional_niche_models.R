# This script fits climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA occurrence data


library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)


# focal species list
spp <- c("Acer macrophyllum", "Pseudotsuga menziesii", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

# climate data
clim <- list.files("big_data/climate", full.names=T) %>%
      raster()

# clip climate data to study area
WNAx <- extent(clim)
WNAx@xmin <- -150
WNAx@xmax <- -90
WNAx@ymin <- 20
WNAx@ymax <- 80
clim <- crop(clim,WNAx)
plot(clim)

# sample background pseudoabsence data in propotion to spatial climate frequency

# load FIA data for focal species

# fit a binomial GAM describing species occurrence as a function of climate

# find the mean and mode of the response surface


