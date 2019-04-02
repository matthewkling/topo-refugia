
library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)


spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Lithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")


# load Little range maps
spf <- list.files("F:/little_trees/raw_data", 
                  recursive=T, full.names=T, pattern="shp")
spf <- lapply(spp, function(x) spf[grepl(x, spf)])

d <- spf %>% 
      lapply(readOGR) %>%
      lapply(gSimplify, tol=.01) %>%
      lapply(broom::tidy)

for(i in 1:length(d)) d[[i]]$gs <- spp[i]
d <- do.call("rbind", d)
d$gs <- sub("Lithocarpus", "Notholithocarpus", d$gs)
md <- map_data("usa")

w2d <- c("Pseudotsuga menziesii", "Acer macrophyllum", "Notholithocarpus densiflorus", 
         "Sequoia sempervirens",  "Arbutus menziesii", "Umbellularia californica", 
         "Quercus garryana", "Quercus kelloggii", "Quercus agrifolia",
         "Quercus lobata", "Aesculus californica","Quercus douglasii")

p <- ggplot(d %>% mutate(gs = factor(gs, levels=w2d)), 
            aes(long, lat, group=paste(gs, group))) +
      geom_polygon(data=md, aes(long, lat, group=group), 
                   color=NA, fill="gray90") +
      geom_polygon() +
      coord_cartesian(xlim=c(-125, -118), ylim=c(34, 45)) +
      annotate(geom="point", y=38.576906, x=-122.703292, color="red", size=3) +
      theme_void() +
      facet_wrap(~gs) +
      theme(legend.position="none")
ggsave("figures/map_ranges_faceted.png", 
       p, width=8, height=8, units="in")

ggplot(d, aes(long, lat, group=paste(gs, group), color=gs)) +
      geom_polygon(fill=NA) +
      coord_cartesian(xlim=c(-125, -118), ylim=c(34, 45)) +
      annotate(geom="point", y=38.576906, x=-122.703292, color="black", shape=10, size=4) +
      theme_void()

ggplot(d, aes(long, lat, group=paste(gs, group))) +
      geom_polygon(fill="black", color=NA, alpha=.1) +
      coord_cartesian(xlim=c(-125, -118), ylim=c(33, 45)) +
      annotate(geom="point", y=38.576906, x=-122.703292, 
               color="orangered", shape=21, size=4) +
      theme_void()


# cwd map
cwd <- stack("big_data/climate/climate.tif")[[1]] %>%
      rasterToPoints() %>% 
      as.data.frame()

p <- ggplot() +
      geom_raster(data=cwd %>% filter(x > -125, x < -118, y > 33, y < 45),
                  aes(x, y, fill=climate.1)) +
      scale_fill_gradientn(colours=c("darkblue", "forestgreen", "yellow3", 
                                     "darkgoldenrod2", "chocolate")) +
      annotate(geom="point", y=38.576906, x=-122.703292, 
               color="red", shape=21, size=4) +
      coord_cartesian(xlim=c(-125, -118), ylim=c(33, 45),
                      expand=c(0,0)) +
      theme_void() +
      theme(legend.position=c(.3, .1),
            legend.direction="horizontal") +
      labs(fill="CWD (mm)")
ggsave("figures/map_cwd.png", 
       p, width=6, height=8, units="in")


r <- spf %>% 
      lapply(readOGR) %>%
      lapply(gSimplify, tol=.01) %>%
      lapply(rasterize, y=raster("big_data/climate/climate.tif")) %>%
      stack() %>%
      reclassify(c(0, Inf, 1)) %>%
      sum(na.rm=T)

rd <- r %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      filter(layer > 0)

md <- map_data("usa")

p <- ggplot() +
      geom_polygon(data=md, aes(long, lat, group=group), 
                   color=NA, fill="gray95") +
      geom_raster(data=rd %>% filter(x > -125, x < -118, y > 33, y < 45),
                  aes(x, y, fill=layer)) +
      scale_fill_gradientn(colours=c("lightyellow", "orange", "red", "darkred"),
                           breaks=seq(1, 12, 3), limits=c(0, 12)) +
      annotate(geom="point", y=38.576906, x=-122.703292, 
               color="cyan", shape=1, size=4) +
      annotate(geom="point", y=38.576906, x=-122.703292, 
               color="cyan", shape=3, size=6) +
      coord_cartesian(xlim=c(-125, -118), ylim=c(33, 45),
                      expand=c(0,0)) +
      theme_void() +
      theme(legend.position=c(.3, .1),
            legend.direction="horizontal") +
      labs(fill="focal species\nrichness")
ggsave("figures/map_species_richness.png", 
       p, width=6, height=8, units="in")


