# This script fits climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA occurrence data


library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)
library(rgeos)
library(mgcv)


# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

# climate data
clim <- list.files("big_data/climate", full.names=T) %>%
      raster()

# clip climate data to conterminous US, to match FIA
usa <- raster::getData("GADM", country="USA", level=1)
usa <- usa[! usa$NAME_1 %in% c("Alaska", "Hawaii"),]
clim <- crop(clim, usa)
clim <- mask(clim, usa)
writeRaster(clim, "big_data/climate/cwd_masked.tif", overwrite=T)

# load FIA data for focal species
# note: leaving the data in granular tree form, to model abundance
fia <- read_csv("data/fia.csv")

# loop through species
results <- data.frame()
for(sp in unique(fia$gs)){
      message(sp)
      
      # species occurrences
      fs <- filter(fia, gs==sp,
                   lon < -95)
      
      # sample background pseudoabsence data in propotion to spatial climate frequency,
      # from buffer region around occurrence points
      hull <- fs %>%
            dplyr::select(lon, lat) %>%
            chull() %>%
            slice(fs, .) %>%
            dplyr::select(lon, lat) %>%
            as.matrix()
      domain <- SpatialPolygons(list(Polygons(list(Polygon(hull)), ID=1))) %>%
            SpatialPolygonsDataFrame(data=data.frame(ID=1)) %>%
            gBuffer(width=5) # in degrees
      bg <- clim %>%
            mask(domain) %>% 
            crop(domain) %>% 
            values() %>% 
            na.omit() %>%
            sample(10000)
      
      # presences
      fs <- sample_n(fs, min(nrow(fs), 10000))
      coordinates(fs) <- c("lon", "lat")
      cs <- extract(clim, fs)
      
      # fit a binomial GAM describing species occurrence as a function of climate
      md <- data.frame(clim=cs, pres=1) %>%
            rbind(data.frame(clim=bg, pres=0))
      fit <- gam(pres ~ s(clim), data=md, family=binomial(logit))
      
      # grid approximation of response surface, summary stats
      results <- data.frame(clim=seq(0, round(max(values(clim), na.rm=T)), 1)) %>%
            mutate(fit=predict(fit, ., type="response"),
                   fit=fit/sum(fit, na.rm=T),
                   species=sp) %>%
            rbind(results)
}


# summarize and plot results
r <- results %>%
      group_by(species) %>%
      summarize(mean=weighted.mean(clim, fit),
                median=clim[cumsum(fit)>.5][1],
                max=clim[fit==max(fit)][1]) %>%
      arrange(mean) %>%
      mutate(species = factor(species, levels=species)) %>%
      gather(stat, clim, mean:max)

results <- mutate(results, species=factor(species, levels=levels(r$species)))

p <- ggplot() +
      geom_area(data=results, aes(clim, fit), fill="gray") +
      geom_vline(data=r, aes(xintercept=clim, color=stat)) +
      geom_text(data=r, aes(x=1500, y=max(results$fit)*.5, label=species),
                hjust=1) +
      facet_grid(species~.) +
      xlim(0, 1500) +
      labs(x="CWD (mm)",
           y="relative frequency") +
      theme_minimal() +
      theme(legend.position="bottom", panel.grid=element_blank(),
            strip.text=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank())
ggsave("figures/regional_response_curves.png", p, width=8, height=8, units="in")




