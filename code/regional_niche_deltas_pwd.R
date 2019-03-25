
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
library(randomForest)
library(dismo)

select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")


# load fia data for focal species
fia <- read_csv("data/fia.csv") %>%
      group_by(plt_cn, lat, lon, gs) %>%
      summarize(n=n()) %>%
      ungroup()
      

# load harbarium species occurrences for non-fia focal species
cch_spp <- c("Adenostoma fasciculatum")
cch <- read_csv("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv") %>%
      rename(gs = current_name_binomial, lon = longitude, lat = latitude) %>%
      select(gs, lon, lat) %>%
      filter(gs %in% cch_spp)




## historic climate data used to fit models
climate <- stack("f:/chelsa/toporefugia/historic/derived/historic.gri")
climate$ppt <- log10(climate$ppt)

ext <- extent(climate)
ext[2] <- -115
clim <- crop(climate, ext)
cor(values(clim), use="pairwise.complete.obs") %>% "^"(2) %>%
      corrplot::corrplot(order="AOE")

var_sets <- list(a=c("cwd", "aet", "tminmin"),
                 b=c("cwd", "ppt", "djf", "jja"))
for(i in 1:length(var_sets)) names(var_sets)[i] <- paste(var_sets[[i]], collapse="_")
clim <- climate[[unique(unlist(var_sets))]]




## pwd climate data to project suitability
hyp <- readRDS("data/HYPshapefile-teale-albers/hyp.boundary.Rdata") %>%
      spTransform(crs(clim))
scen_files <- c("f:/chelsa/toporefugia/historic/derived/historic.gri",
                list.files("f:/chelsa/toporefugia/future/derived", full.names=T, pattern="\\.gri"))
scenarios <- scen_files %>%
      lapply(stack) %>%
      lapply(subset, subset=unique(unlist(var_sets))) %>%
      lapply(crop, y=hyp) %>%
      lapply(mask, mask=hyp) %>%
      lapply(function(x){
            if(! "ppt" %in% names(x)) return(x)
            x$ppt <- log10(x$ppt)
            return(x)
      }) %>%
      lapply(values) %>%
      lapply(as.data.frame) %>%
      lapply(na.omit) %>%
      lapply(function(x) mutate(x, cell=1:nrow(x)))
names(scenarios) <- basename(scen_files) %>% sub("\\.gri", "", .)
for(i in 1:length(scenarios)) scenarios[[i]]$scenario <- names(scenarios)[i]
scenarios <- do.call("rbind", scenarios)


# loop through species
results <- data.frame()
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
      #fs <- sample_n(fs, min(nrow(fs), 10000))
      fs <- sample_n(fs, 10000, replace=T)
      coordinates(fs) <- c("lon", "lat")
      cs <- extract(sp_clim, fs) %>%
            na.omit()
      
      md <- rbind(cs %>% cbind(pres=1),
                  bg %>% cbind(pres=0))
      
      
      for(vars in names(var_sets)){
            
            for(algo in c("gam", "maxent", "randomforest")){
                  message(paste("      ", vars, "     ", algo))
                  
                  if(algo=="gam"){
                        formula <- as.formula(paste0("pres ~ ", paste0("s(", var_sets[[vars]], ")", collapse=" + ")))
                        fit <- gam(formula, md, family=binomial(logit))
                  }
                  if(algo=="randomforest"){
                        formula <- as.formula(paste0("pres ~ ", paste0(var_sets[[vars]], collapse=" + ")))
                        fit <- randomForest(formula=formula, data=md %>% mutate(pres=as.logical(pres)))
                  }
                  if(algo=="maxent") fit <- maxent(md[,var_sets[[vars]]], md$pres,
                                                   args=c("-a", "-z", "outputformat=raw", 
                                                          "nothreshold", "nohinge"))
                  
                  # model predictions
                  results <- scenarios %>%
                        mutate(algorithm=algo,
                               var_set=vars,
                               species=sp,
                               pred=predict(fit, scenarios, type="response")) %>%
                        rbind(results, .)
            }
      }
}


write_csv(results, "data/pwd_niche_values.csv")

deltas <- results %>%
      select(cell:pred) %>%
      spread(scenario, pred) %>%
      gather(scenario, pred, -cell, -algorithm, 
             -var_set, -species, -historic) %>%
      mutate(delta = pred - historic)
      
write_csv(deltas, "data/pwd_niche_deltas.csv")

deltas %>%
      select(cell:species, scenario, delta) %>%
      spread(algorithm, delta) %>%
      ggplot(aes(randomforest, maxent, color=species)) + 
      geom_point() +
      facet_grid(scenario~var_set, scales="free")
