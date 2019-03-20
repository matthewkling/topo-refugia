
# Attenuata-type scatterplots for all pwd species. 
# Occupancy heat maps.

library(tidyverse)
library(raster)

# climate data
clim <- stack("f:/chelsa/toporefugia/historic/derived/historic.gri")
clim$ppt <- log10(clim$ppt)

# load pepperwood boundary and get its climate
pwd <- readOGR("data/PPshapefile-teale-albers", "Pepperwood")
crs(pwd) <- crs("+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0 ")
pwd <- pwd %>% spTransform(crs(clim)) %>%
      extract(clim, .)
pwd <- pwd[[1]]
write_csv(as.data.frame(pwd), "data/pwd_climate_1km.csv")

# load FIA and add climate data to each record
fia <- read_csv("data/fia.csv") %>%
      filter(!is.na(lon))
coordinates(fia) <- c("lon", "lat")
clim <- clim$cwd
fia$cwd <- extract(clim, fia)
fia <- as.data.frame(fia)


# calculate northness
d <- fia %>%
      mutate(aspect_deg = 360-aspect,
             aspect_rad = aspect_deg/360*2*pi,
             slope_rad = atan(slope/100),
             slope_deg = slope_rad * 180 / pi,
             northness = cos(aspect_rad) * sin(slope_rad),
             southness = -northness) %>%
      filter(is.finite(northness))

# southness~cwd regression line plots
p <- d %>% 
      group_by(species) %>% 
      #filter(ecdf(cwd)(cwd) > .025,
      #       ecdf(cwd)(cwd) < .975) %>%
      ggplot(aes(cwd, southness, color=gs)) +
      geom_vline(xintercept=pwd[[1]], alpha=.5) +
      geom_smooth(method=lm, se=F) +
      theme_minimal() +
      labs(color=NULL, 
           x="CWD (mm); vertical lines are Pepperwood 1km pixels")
ggsave("figures/regional_southness_cwd_lm.png", p, width=8, height=8, units="in")


# southness~cwd scatterplots
w2d <- c("Pseudotsuga menziesii", "Acer macrophyllum", "Notholithocarpus densiflorus", 
         "Sequoia sempervirens",  "Arbutus menziesii", "Umbellularia californica", 
         "Quercus garryana", "Quercus kelloggii", "Quercus agrifolia",
         "Quercus lobata", "Aesculus californica","Quercus douglasii")
sd <- d %>% 
      mutate(gs=factor(gs, levels=w2d)) %>% # arrange wet to dry
      group_by(plot, subp, gs) %>%
      summarize(n=n(),
                cwd=mean(cwd, na.rm=T),
                southness=mean(southness))
p <- sd %>%
      ggplot(aes(cwd, southness, size=n)) +
      geom_vline(xintercept=pwd[[1]], alpha=.5, color="dodgerblue") +
      geom_point(alpha=.25) +
      geom_smooth(method=lm, se=F, color="darkred") +
      facet_wrap(~gs, nrow=6) +
      theme_minimal() +
      theme(legend.position="none") +
      labs(color=NULL, 
           x="CWD (mm); vertical lines are Pepperwood 1km pixels",
           y="subplot southness")
ggsave("figures/regional_southness_cwd_scatter.png", p, width=8, height=8, units="in")





# next todo: occurrence heatmaps, relative to climate background

# do this with 2d GAM, or bins, or what? GAM would be most elegant/consistent
