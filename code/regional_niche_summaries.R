

# This script fits univariate climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA and CCH occurrence data, and uses them to
# calculate niche summary statistics for each species for each climate variable.


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


# climate data
clim <- stack("f:/chelsa/toporefugia/historic/derived/historic.gri")
clim$ppt <- log10(clim$ppt)


# load fia data for focal species
# note: leaving the data in granular tree form, to model abundance
fia <- read_csv("data/fia.csv")

# load harbarium species occurrences for non-fia focal species
cch_spp <- c("Adenostoma fasciculatum")
cch <- read_csv("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv") %>%
      rename(gs = current_name_binomial, lon = longitude, lat = latitude) %>%
      select(gs, lon, lat) %>%
      filter(gs %in% cch_spp)



# loop through species
results <- data.frame()
#dev <- c()
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
      
      # fit a separate model for each variable
      for(var in names(clim)){
            
            # fit a binomial GAM describing species occurrence as a function of climate
            md <- data.frame(value=cs[,var], pres=1) %>%
                  rbind(data.frame(value=bg[,var], pres=0))
            fit <- gam(pres ~ s(value), data=md, family=binomial(logit))
            #dev <- c(dev,summary(fit)$dev.expl)
            
            # grid approximation of response surface
            results <- data.frame(var=var, 
                                  value=seq(min(md$value), max(md$value), length.out=1000)) %>%
                  mutate(fit=predict(fit, ., type="response"),
                         fit=fit/sum(fit, na.rm=T),
                         species=sp) %>%
                  rbind(results)
      }
}


# summarize and plot results
r <- results %>%
      group_by(species, var) %>%
      summarize(mean=weighted.mean(value, fit),
                median=value[cumsum(fit)>.5][1],
                max=value[fit==max(fit)][1]) %>%
      ungroup()

r

write.csv(r, "data/regional_niche_stats.csv", row.names = F)


for(v in unique(r$var)){
      rv <- r %>%
            filter(var==v) %>%
            arrange(mean) %>%
            mutate(species = factor(species, levels=species)) %>%
            gather(stat, value, mean:max)
      
      resultsv <- results %>%
            filter(var==v) %>%
            mutate(species=factor(species, levels=levels(rv$species)))
      
      p <- ggplot() +
            geom_area(data=resultsv, aes(value, fit), fill="gray") +
            geom_vline(data=rv, aes(xintercept=value, color=stat)) +
            geom_text(data=rv, aes(x=max(resultsv$value), y=max(resultsv$fit)*.5, label=species),
                      hjust=1) +
            facet_grid(species~.) +
            labs(x=v,
                 y="relative frequency") +
            theme_minimal() +
            theme(legend.position="bottom", panel.grid=element_blank(),
                  strip.text=element_blank(), axis.ticks.y=element_blank(),
                  axis.text.y=element_blank(), axis.title.y=element_blank())
      ggsave(paste0("figures/regional_response_curves/", v, ".png"), 
             p, width=8, height=6, units="in")
}



