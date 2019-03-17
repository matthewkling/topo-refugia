# This script fits climatic niche models for Pepperwood tree species
# using CHELSA climate data and FIA and CCH occurrence data


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
if(F){
      clim <- list.files("big_data/climate/raw", full.names=T) %>%
            stack()
      names(clim) <- c("cwd", "ppt", "bio5", "bio6")
      usa <- raster::getData("GADM", country="USA", level=1)
      usa <- usa[! usa$NAME_1 %in% c("Alaska", "Hawaii"),]
      clim <- crop(clim, usa)
      clim <- mask(clim, usa)
      clim <- crop(clim, extent(-124.7668, -95, 24.51653, 49.38319 ))
      clim$ppt <- log10(clim$ppt)
      writeRaster(clim, "big_data/climate/climate.tif", overwrite=T)
}
clim <- stack("big_data/climate/climate.tif")
names(clim) <- c("cwd", "ppt", "bio5", "bio6")
clim <- clim$cwd


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
      
      # grid approximation of response surface
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

write.csv(r, "data/regional_niche_stats.csv", row.names = F)




#################### multivariate ###############

stop("incomplete code below")


# loop through species
models <- list()
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
            crop(domain) %>% 
            mask(domain) %>% 
            values() %>% 
            na.omit() %>%
            as.data.frame() %>%
            sample_n(10000) %>%
            mutate(pres=0) %>%
            sample_n(min(nrow(fs), 10000))
      
      # presences
      fs <- sample_n(fs, min(nrow(fs), 10000))
      coordinates(fs) <- c("lon", "lat")
      cs <- extract(clim, fs) %>%
            as.data.frame() %>%
            mutate(pres=1)
      
      md <- rbind(bg, cs)
      
      # fit a binomial GAM describing species occurrence as a function of climate
      formula <- as.formula(paste0("pres ~ ", paste0("s(", names(clim), ")", collapse=" + ")))
      fit <- gam(formula, data=md, family=binomial(logit), select=T)
      #vis.gam(fit, c("cwd", "ppt"), theta=50, phi=50, type="response")
      
      # grid approximation of response surface, summary stats
      n <- 31
      g <- data.frame(row=1:n)
      for(i in names(select(bg, -pres))){
            g[,i] <- seq(min(bg[,i]), max(bg[,i]), length.out=n)
      }
      g <- g %>% select(-row) %>% 
            expand.grid() %>%
            mutate(fit = predict(fit, ., type="response"))
      
      results <- g %>%
            mutate(fit=fit/sum(fit, na.rm=T),
                   species=sp) %>%
            rbind(results)
      
      next()
      
      g %>%
            filter(bio6 == bio6[fit==max(fit)],
                   bio5 == bio5[fit==max(fit)]) %>%
            ggplot(aes(cwd, ppt, fill=fit)) + 
            geom_raster()
      
      g %>%
            group_by(cwd, ppt) %>%
            summarize(fit = mean(fit)) %>%
            ggplot(aes(cwd, ppt, fill=fit)) + 
            geom_raster()
}


# summarize and plot results
mean <- results %>%
      group_by(species) %>%
      summarize_at(vars(cwd, ppt, bio5, bio6), 
                   list(mean=weighted.mean), w=.$fit) %>%
      gather(stat, value, -species)

median <- results %>%
      group_by(species) %>%
      summarize_at(vars(cwd, ppt, bio5, bio6), 
                   list(median=clim[cumsum(fit)>.5][1])) %>%
      gather(stat, value, -species)

mean <- results %>%
      group_by(species) %>%
      summarize_at(vars(cwd, ppt, bio5, bio6), 
                   list(mean=weighted.mean), w=.$fit) %>%
      gather(stat, value, -species)


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

write.csv(r, "data/regional_niche_stats.csv", row.names = F)
