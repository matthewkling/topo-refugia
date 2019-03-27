
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
for(sp in unique(c(fia$gs, as.character(cch$gs)))){
      message(sp)
      
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
              message(paste("    ",round(summary(fit)$dev.expl,4)))
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

# check results are readable
mean(values(raster(outfile)))

## extract results
dcwds <- c(0,40,80,120)
daets <- c(-25,0,25,50)
dtminmins <- c(0,1,2,3)

infiles <- list.files('data/pwd_distributions/rasters')
length(infiles)
length(sp)
gam_pwd_suit <- data.frame(fname=infiles,sp=NA,var_set=NA,scen=NA,dcwd=NA,daet=NA,dtmintmin=NA,suit=NA)

hypx <- readRDS('data/HYPshapefile-geo/hyp.boundary.Rdata')

i=1
for (i in 1:length(infiles)){
  infile <- infiles[i]
  ras <- raster(paste0('data/pwd_distributions/rasters/',infile))
  chs <- strsplit(infile,'__')
  gam_pwd_suit$sp[i] <- chs[[1]][1]
  gam_pwd_suit$var_set[i] <- chs[[1]][2]
  scen <- substr(chs[[1]][3],1,7)
  gam_pwd_suit$scen[i] <- scen
  gam_pwd_suit$dcwd[i] <- dcwds[match(substr(scen,5,5),c('A','B','C','D'))]
  gam_pwd_suit$daet[i] <- daets[match(substr(scen,6,6),c('A','B','C','D'))]
  gam_pwd_suit$dtmintmin[i] <- dtminmins[match(substr(scen,7,7),c('A','B','C','D'))]
  gam_pwd_suit$suit[i] <- mean(extract(ras,hypx)[[1]])
}
head(gam_pwd_suit)
tail(gam_pwd_suit)

# add column for historical suitability (ABA) to calculate changes in suitability
hsuit <- gam_pwd_suit[gam_pwd_suit$scen=='climABA',]
hsuit
gam_pwd_suit$hsuit <- hsuit$suit[match(gam_pwd_suit$sp,hsuit$sp)]
head(gam_pwd_suit)
tail(gam_pwd_suit)
gam_pwd_suit$dsuit <- gam_pwd_suit$suit - gam_pwd_suit$hsuit
hist(gam_pwd_suit$dsuit)
summary(gam_pwd_suit$dsuit)

species <- unique(gam_pwd_suit$sp)
climcf <- data.frame(sp=species,cwd.cf=NA,aet.cf=NA,tminmin.cf=NA)

i=2
for (i in 1:length(species)){
  message(species[i])
  xx <- gam_pwd_suit[gam_pwd_suit$sp==species[i],]
  fit <- glm(suit~dcwd+daet+dtmintmin,data=xx)
  summary(fit)$coeff
  climcf$sp[i] <- species[i]
  climcf$cwd.cf[i] <- summary(fit)$coeff[2,1] * 120
  climcf$aet.cf[i] <- summary(fit)$coeff[3,1] * 75
  climcf$tminmin.cf[i] <- summary(fit)$coeff[4,1] * 4
  plot(xx$suit[order(xx$daet,xx$dtmintmin,xx$dcwd)])
}
climcf
barplot(t(as.matrix(climcf[,-1])),beside=T)

topo <- read.csv('data/pwd_niche_means.csv')
head(topo)
d <- merge(climcf,topo,by.x='sp',by.y='sci.names')
d

plot(d$cwd.mean,d$cwd.cf)
summary(lm(d$cwd.cf~d$cwd.mean))
abline(h=0)
sort(d$cwd.cf)

plot(d$south.mean,d$aet.cf)
summary(lm(d$aet.cf~d$cwd.mean))
abline(h=0)

plot(d$south.mean,d$tminmin.cf)
abline(h=0)

plot(d$topoid.mean,d$tminmin.cf)
abline(h=0)

pairs(d[,c('cwd.mean','cwd.cf','aet.cf','tminmin.cf')])

##### After GAM models have been fit, they can be reloaded here without going through the entire loop below - just to
if (FALSE) {
  for(sp in unique(c(fia$gs, as.character(cch$gs)))){
    outfile <- paste0("data/pwd_distributions/models/",
                      sp, "_", vars, ".rds")
    fit <- readRDS(outfile) 
    message(paste(sp,round(summary(fit)$dev.expl,4)))
  }
}
#GAM dev.expl values with var_set 'a'
# $a
# [1] "cwd"     "aet"     "tminmin"

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