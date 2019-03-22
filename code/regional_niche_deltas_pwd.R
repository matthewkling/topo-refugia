

# regional niche deltas at pepperwood

r <- list.files("data/regional_distributions/rasters/",
                full.names=T)

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica",
         "Adenostoma fasciculatum")

# load pepperwood boundary and get its climate
pwd <- readOGR("data/PPshapefile-teale-albers", "Pepperwood")
crs(pwd) <- crs("+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0 ")
pwd <- pwd %>% spTransform(crs(raster(r[1])))


results <- data.frame()
for(s in spp){
      message(s)
      
      h <- raster(r[grepl(s, r) & grepl("historic", r)]) %>% crop(pwd)
      rf <- r[grepl(s, r) & grepl("rcp", r)]
      f <- stack(rf) %>% crop(pwd)
      d <- f - h
      v <- extract(d, pwd)[[1]] %>%
            as.data.frame()
      names(v) <- rf %>% basename() %>% 
            str_split("__", simplify=T) %>% 
            "["(,3) %>% 
            sub("\\.tif", "", .)
      v <- v %>%
            mutate(species=s, cell = 1:nrow(.)) %>%
            gather(scenario, delta, -cell, -species)
      
      results <- rbind(results, v)
}

write_csv(results, "data/pwd_niche_deltas.csv")
