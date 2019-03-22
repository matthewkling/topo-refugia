
# This script fits multivariate climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA and CCH occurrence data, and projects them to
# various future scenarios.


library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)
library(rgeos)
library(mgcv)

select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")


# load fia data for focal species
# note: leaving the data in granular tree form, to model abundance
fia <- read_csv("data/fia.csv")

# load harbarium species occurrences for non-fia focal species
cch_spp <- c("Adenostoma fasciculatum")
cch <- read_csv("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv") %>%
      rename(gs = current_name_binomial, lon = longitude, lat = latitude) %>%
      select(gs, lon, lat) %>%
      filter(gs %in% cch_spp)


## climate data used to fit models
climate <- stack("f:/chelsa/toporefugia/historic/derived/historic.gri")
climate$ppt <- log10(climate$ppt)

ext <- extent(climate)
ext[2] <- -115
clim <- crop(climate, ext)
cor(values(clim), use="pairwise.complete.obs") %>% "^"(2) %>%
      corrplot::corrplot(order="AOE")

var_sets <- list(a=c("cwd", "aet", "tminmin"))
clim <- climate[[unique(unlist(var_sets))]]


## climate datsets for which to project suitability
scen_files <- c("f:/chelsa/toporefugia/historic/derived/historic.gri",
                list.files("f:/chelsa/toporefugia/future/derived", full.names=T, pattern="\\.gri"))
scenarios <- scen_files %>%
      lapply(stack) %>%
      lapply(subset, subset=unique(unlist(var_sets))) %>%
      lapply(function(x){
            if(! "ppt" %in% names(x)) return(x)
            x$ppt <- log10(x$ppt)
            return(x)
      }) %>%
      lapply(values) %>%
      lapply(as.data.frame)
names(scenarios) <- basename(scen_files) %>% sub("\\.gri", "", .)



# loop through species
for(sp in unique(c(fia$gs, cch$gs))){
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
      
      for(vars in names(var_sets)){
            
            # fit a binomial GAM describing species occurrence as a function of climate
            formula <- as.formula(paste0("pres ~ ", paste0("s(", var_sets[[vars]], ")", collapse=" + ")))
            fit <- gam(formula, data=md, family=binomial(logit))
            saveRDS(fit, paste0("data/regional_distributions/models/",
                                sp, "_", vars, ".rds"))
            
            # model predictions (the slow step)
            for(scen in names(scenarios)){
                  outfile <- paste0("data/regional_distributions/rasters/",
                                    sp, "__", vars, "__", scen, ".tif")
                  if(file.exists(outfile)) next()
                  message(paste("      ", scen))
                  
                  pred <- predict(fit, scenarios[[scen]], type="response") %>% as.vector()
                  template <- climate[[1]] 
                  template[] <- pred
                  writeRaster(template, 
                              outfile,
                              overwrite=T)
            }
      }
}



