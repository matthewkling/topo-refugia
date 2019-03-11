
library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)


# load focal species list
spp <- c("Acer macrophyllum", "Pseudotsuga menziesii", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

# load climate data
clim_files <- list.files("E:/phycon/data/bcm_normals", pattern="filled3x", full.names=T) 
clim_vars <- substr(basename(clim_files), 1, 3)
clim <- lapply(clim_files, readRDS) %>% do.call("stack", .)
names(clim) <- clim_vars
clim <- clim[[sort(names(clim))]]
clim_fut <- readRDS("e:/phycon/data/bcm_normals/BCM_FUT_810m_filled/CCSM4_rcp85_2040_2069.rds")

# load pepperwood boundary
pwd <- readOGR("e:/fia/PPshapefile-teale-albers", "Pepperwood")
crs(pwd) <- crs(clim)

# load occurrence and background data
occ <- fread("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv") %>%
      dplyr::select(clade, current_name_binomial, longitude, latitude) %>%
      filter(current_name_binomial %in% spp)
bg <- readRDS("E:/phycon/data/occurrences/10000_CA_plants_bg_810m_occbias.rdata")




# cluster setup
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)


results <- foreach(pres = split(occ, occ$current_name_binomial), .combine="c") %dopar% {
      
      options(java.parameters = "-Xmx1g" )
      require(dismo)
      require(raster) 
      
      # format presence data
      pres <- as.data.frame(pres)
      coordinates(pres) <- c("longitude", "latitude")
      projection(pres) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
      pres <- spTransform(pres, crs(clim))
      
      # drop occurrences that fall in the water
      valid <- !is.na(raster::extract(clim[[1]], pres))
      pres <- pres[valid,]
      
      # fit and project maxent model
      mx <- try(maxent(clim, pres, bg,
                       args=c("-a", "-z", "outputformat=raw", "maximumbackground=10000", 
                              "nothreshold", "nohinge")))
      if(class(mx)=="try-error") return("failure")
      suit <- predict(mx, clim)
      suit_fut <- predict(mx, clim_fut)
      
      # fit and project distance model
      sigma <- 50 # km
      dst <- distanceFromPoints(clim[[1]], pres)
      dst <- mask(dst, clim[[1]])
      gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
      prox <- calc(dst, function(x) gauss(x, sigma*1000))
      
      # combine and save
      pred <- suit
      pred_fut <- suit_fut
      saveRDS(pred, paste0("E:/fia/pepperwood/ranges/present/", pres$current_name_binomial[1], ".rds"))
      saveRDS(pred_fut, paste0("E:/fia/pepperwood/ranges/future/", pres$current_name_binomial[1], ".rds"))
      saveRDS(prox, paste0("E:/fia/pepperwood/ranges/distance_kernel/", pres$current_name_binomial[1], ".rds"))
      saveRDS(dst, paste0("E:/fia/pepperwood/ranges/distance/", pres$current_name_binomial[1], ".rds"))
      return("success")
}



# values at pepperwood

p <- foreach(x=list.files("E:/fia/pepperwood/ranges/distance/", full.names=T), 
             .combine="rbind", .packages=c("raster", "tidyverse")) %dopar% {
                   
                   spp <- basename(x)
                   paste0(c("E:/fia/pepperwood/ranges/present/",
                            "E:/fia/pepperwood/ranges/future/",
                            "E:/fia/pepperwood/ranges/distance/",
                            "E:/fia/pepperwood/ranges/distance_kernel/"), 
                          basename(x)) %>%
                         lapply(readRDS) %>%
                         do.call("stack", .) %>%
                         crop(pwd) %>%
                         mask(pwd) %>%
                         rasterToPoints() %>%
                         as.data.frame() %>%
                         rename(suit_hist=layer.1, suit_fut=layer.2,
                                distance_m=layer.3, distance_kernel=layer.4) %>%
                         mutate(species=sub("\\.rds", "", spp)) %>%
                         dplyr::select(-x, -y) %>%
                         group_by(species) %>%
                         summarize_all(funs(min, median, max)) %>%
                         gather(var, value, -species) %>%
                         separate(var, c("v1", "v2", "stat")) %>%
                         unite(variable, v1, v2) %>%
                         spread(variable, value) %>%
                         mutate(suit_ratio = suit_fut / suit_hist,
                                suit_diff = suit_fut - suit_hist)
             }

write.csv(p, "e:/fia/pepperwood/pepperwood_results.csv", row.names=F)

stopCluster(cl)
