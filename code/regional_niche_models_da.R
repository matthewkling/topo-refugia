
# This script fits multivariate climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA and CCH occurrence data, and projects them to
# various future scenarios.

rm(list=ls())
library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)
library(rgeos)
library(mgcv)
library(corrplot)
library(dismo)
library(rJava)

select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")


# load fia data for focal species
# note: leaving the data in granular tree form, to model abundance
fia <- read_csv("data/fia.csv")
dim(fia)

# load harbarium species occurrences for non-fia focal species
cch_spp <- c("Adenostoma fasciculatum")
cch <- readRDS("/Users/david/Documents/Projects/CalDiversity/CCH_2016_master/California_Species_clean_All_epsg_3310.Rdata") %>%
      rename(gs = current_name_binomial, lon = longitude, lat = latitude) %>%
      select(gs, lon, lat) %>%
      filter(gs %in% cch_spp)
head(cch)

## climate data used to fit models
climate <- stack("big_data/climate/historic.gri")
climate$ppt <- log10(climate$ppt)
names(climate)

ext <- extent(climate)
ext[2] <- -115
clim <- crop(climate, ext)
cor(values(clim), use="pairwise.complete.obs") %>% "^"(2) %>%
      corrplot::corrplot(order="AOE")

var_sets <- list(a=c("cwd", "aet", "tminmin"))
clim <- climate[[unique(unlist(var_sets))]]

## climate datsets for which to project suitability
#scen_files <- c("f:/chelsa/toporefugia/historic/derived/historic.gri",
#                list.files("f:/chelsa/toporefugia/future/derived", full.names=T, pattern="\\.gri"))
scen_files <- c(list.files("data/climate/future", full.names=T, pattern="\\.Rdata"))
scen_names <- substr(scen_files,21,27)

scenarios <- list()
for (i in 1:length(scen_files)) scenarios[[i]] <- readRDS(scen_files[i])
names(scenarios) <- scen_names

head(data.frame(values(scenarios[[1]])))

algo <- 'gam' # choose modeling algorithm

# loop through species
sp <- 'Aesculus californica' # test value to debug inside loop
res <- c()
for(sp in sort(unique(c(fia$gs, as.character(cch$gs))))){
      #message(sp)
      
      # species occurrences
      fs <- switch(as.character(sp %in% cch_spp), "TRUE"=cch, "FALSE"=fia) %>%
            filter(gs==sp, lon < -95)
      
      # sample background pseudoabsence data in propotion to spatial climate frequency,
      # from buffer region around occurrence points
      hull <- fs %>%
            dplyr::select(lon, lat) %>%
            chull() %>% # fitting convex hull first to avoid buffering huge point set
            slice(fs, .) %>%
            dplyr::select(lon, lat) %>%
            as.matrix()
      domain <- SpatialPolygons(list(Polygons(list(Polygon(hull)), ID=1))) %>%
            SpatialPolygonsDataFrame(data=data.frame(ID=1)) %>%
            gBuffer(width=5) # in degrees
      sp_clim <- clim %>%
            mask(domain) %>% 
            crop(domain)
      bg <- sp_clim %>% 
            values() %>% 
            na.omit() %>%
            as.data.frame() %>%
            sample_n(10000)
      
      # presences
      fs <- sample_n(fs, min(nrow(fs), 10000))
      coordinates(fs) <- c("lon", "lat")
      cs <- extract(sp_clim, fs) %>%
            na.omit()
      
      md <- rbind(cs %>% cbind(pres=1),
                  bg %>% cbind(pres=0))
      
      vars <- "a" # test value for debugging
      for(vars in names(var_sets)){
            
            # fit a binomial GAM describing species occurrence as a function of climate
            if (algo=='gam') {
              formula <- as.formula(paste0("pres ~ ", paste0("s(", var_sets[[vars]], ")", collapse=" + ")))
              outfile <- paste0("data/pwd_distributions/models/",
                                sp, "_", vars, ".rds")
              if(file.exists(outfile)) {
                fit <- readRDS(outfile) 
              } else {
                fit <- gam(formula, data=md, family=binomial(logit))
                saveRDS(fit, paste0("data/pwd_distributions/models/",
                                    sp, "_", vars, ".rds"))
              }
              # what is peak probability value, and how do values compare
              # for presence vs. absence sites. Only for manual examination
              if (FALSE) {
                boxplot(md$pres,fit$fitted.values)
                plot(md$cwd,fit$fitted.values)
                plot(md$aet,fit$fitted.values)
                plot(md$tminmin,fit$fitted.values)
              }
              maxval <- max(fit$fitted.values[md$pres==1],na.rm=T)
              message(paste(sp,round(summary(fit)$dev.expl,4),round(maxval,3)))
              res <- rbind(res,c(sp,round(summary(fit)$dev.expl,4),round(maxval,3)))
            } else {
              fit <- maxent(md[,var_sets[[vars]]], md$pres,
                            args=c("-a", "-z", "outputformat=raw", 
                                   "nothreshold", "nohinge"))
            }
            # model predictions (the slow step)
            scen <- "climAAA" # test value for debugging
            for(scen in names(scenarios)){
                  outfile <- paste0("data/pwd_distributions/rasters/",
                                    sp, "__", vars, "__", scen, ".tif")
                  if(file.exists(outfile)) next()
                  message(paste("      ", scen))
                  
                  pred <- predict(fit, data.frame(values(scenarios[[scen]])), type="response") %>% as.vector()
                  template <- scenarios[[1]]$cwd
                  template[] <- pred
                  names(template) <- 'pred'
                  writeRaster(template, 
                              outfile,
                              overwrite=T)
            }
      }
}

res
write.csv(data.frame(sp=res[,1],dev.expl=as.numeric(res[,2]),max.fit.value=as.numeric(res[,3])),'data/gam_varseta_modelstats.csv')

# check results are readable
mean(values(raster(outfile)))

#var_set 'a'
# $a
# [1] "cwd"     "aet"     "tminmin"

# sp - dev.expl - max fitted.value
# Aesculus californica 0.4985
# Pseudotsuga menziesii 0.4081
# Notholithocarpus densiflorus 0.7168
# Sequoia sempervirens 0.8542
# Arbutus menziesii 0.6613
# Umbellularia californica 0.6495
# Acer macrophyllum 0.6684
# Quercus garryana 0.6008
# Quercus kelloggii 0.6843
# Quercus agrifolia 0.735
# Quercus douglasii 0.7037
# Quercus lobata 0.4923
# Adenostoma fasciculatum 0.5073