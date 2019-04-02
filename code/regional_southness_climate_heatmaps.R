
library(tidyverse)
library(raster)
select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

# climate data, clipped to conterminous US
# (created in regional niche model script)
#clim <- stack("big_data/climate/climate.tif")
clim <- stack("f:/chelsa/toporefugia/historic/derived/historic.gri")

# load FULL fia dataset, summarize to subplot-species level
d <- readRDS("e:/fia/topoclimate/data_munged.rds") %>%
      select(-c(bio1:bio19)) %>%
      mutate(southness = -northness,
             gs = paste(genus, species),
             gs = sub("Lithocarpus densiflorus", 
                      "Notholithocarpus densiflorus", gs)) %>%
      filter(lon < -95) %>%
      group_by(gs, plt_cn, plot, subp, lon, lat) %>%
      summarize(southness = mean(southness),
                n=n())

# climate at each plot
coordinates(d) <- c("lon", "lat")
d$cwd <- extract(clim$cwd, d)
d$tmax <- extract(clim$tmaxmax, d)
d <- as.data.frame(d)



# binned frequencies
d$site <- paste(d$plt_cn, d$plot, d$subp)

sites <- d %>% select(site, cwd, tmax, southness) %>% distinct()

pres <- d %>% select(gs, site) %>% mutate(occ = 1)

b <- d %>%
      expand(gs, site) %>%
      left_join(sites) %>%
      left_join(pres) %>%
      mutate(occ = !is.na(occ)) %>%
      mutate(cwd_bin = plyr::round_any(cwd, 50),
             tmax_bin = plyr::round_any(tmax, 1),
             southness_bin = plyr::round_any(southness, .1)) %>%
      group_by(gs, southness_bin, cwd_bin) %>%
      #group_by(gs, southness_bin, tmax_bin) %>%
      summarize(occ = mean(occ)) %>%
      filter(gs %in% spp)


w2d <- c("Pseudotsuga menziesii", "Acer macrophyllum", "Notholithocarpus densiflorus", 
         "Sequoia sempervirens",  "Arbutus menziesii", "Umbellularia californica", 
         "Quercus garryana", "Quercus kelloggii", "Quercus agrifolia",
         "Quercus lobata", "Aesculus californica","Quercus douglasii")

b <- b %>% 
      ungroup() %>%
      filter(abs(southness_bin)<=.5) %>%
      mutate(gs=factor(gs, levels=w2d))

txt <- b %>%
      select(gs) %>%
      distinct()

p <- ggplot() +
      geom_tile(data=b %>% group_by(gs) %>% mutate(occ = occ/max(occ)), 
                aes(cwd_bin, southness_bin, fill=occ)) +
      geom_text(data=txt, aes(x=1450, y=0, label=sub(" ", "\n", gs)),
                color="black", hjust=1, lineheight=.75) +
      scale_fill_gradientn(colours=c("gray90", "forestgreen", "black"), 
                           #trans="log10", 
                           na.value="gray95") +
      facet_grid(gs~.) +
      scale_y_continuous(breaks=c(-.5, 0, .5)) +
      xlim(0, 1500) +
      labs(x="CWD (mm)",
           y="southness",
           fill="occupancy") +
      theme_minimal() +
      theme(strip.text=element_blank())
ggsave("figures/regional_southness_cwd_occupancy.png", p, width=8, height=8, units="in")




##############################################





# load FULL fia dataset, summarize to subplot-species level
d <- readRDS("e:/fia/topoclimate/data_munged.rds") %>%
      select(-c(bio1:bio19)) %>%
      mutate(southness = -northness,
             gs = paste(genus, species),
             gs = sub("Lithocarpus densiflorus", 
                      "Notholithocarpus densiflorus", gs)) %>%
      filter(lon < -95) %>%
      mutate(site = paste(plt_cn, plot, subp)) %>%
      select(gs, site, lon, lat, southness) %>%
      distinct()

# climate at each plot
sites <- d %>% 
      select(site, lon, lat, southness) %>%
      distinct()
coordinates(sites) <- c("lon", "lat")
sites$cwd <- extract(clim$cwd, sites)
sites$tmax <- extract(clim$tmaxmax, sites)
sites <- as.data.frame(sites)


slopes <- data.frame()
for(s in spp){
      sd <- d %>%
            filter(gs == s) %>%
            left_join(sites, .) %>%
            mutate(pres = !is.na(gs)) %>%
            filter(abs(southness) <= .5)
      
      for(half in c("high", "low")){
            
            if(half == "high") lim <- max(sd$cwd[sd$pres], na.rm=T)
            if(half == "low") lim <- min(sd$cwd[sd$pres], na.rm=T)
            
            if(half == "high") mid <- quantile(sd$cwd[sd$pres], .9, na.rm=T)
            if(half == "low") mid <- quantile(sd$cwd[sd$pres], .1, na.rm=T)
            
            if(half == "high") md <- filter(sd, cwd > mid)
            if(half == "low") md <- filter(sd, cwd < mid)
            
            fit <- glm(pres ~ cwd + southness, data=md,
                       family=binomial(link="logit"))
            
            pred <- expand.grid(cwd=seq(mid, lim, length.out=21),
                                southness = seq(-.5, .5, length.out=21))
            pred$pred <- predict(fit, pred, type="response")
            
            pred <- pred %>%
                  mutate(
                        gs = s,
                        half = half,
                        mid = mid,
                        lim = lim,
                        slope = - coef(fit)["cwd"] / coef(fit)["southness"])
            slopes <- rbind(slopes, pred)
      }
}


p <- ggplot() +
      geom_tile(data= b %>% 
                      group_by(gs) %>% mutate(occ = occ/max(occ)) %>% 
                      ungroup() %>% mutate(gs=factor(gs, levels=w2d)), 
                aes(cwd_bin, southness_bin, fill=occ)) +
      geom_text(data=txt %>% mutate(gs=factor(gs, levels=w2d)), 
                aes(x=1450, y=0, label=sub(" ", "\n", gs)),
                color="black", hjust=1, lineheight=.75, fontface="italic") +
      geom_rect(data=slopes %>% mutate(gs=factor(gs, levels=w2d)), 
                 aes(xmin=mid, xmax=lim, ymin=-.5, ymax=.5),
                color="black", fill=NA) +
      geom_contour(data=slopes %>% mutate(gs=factor(gs, levels=w2d)),
                   aes(cwd, southness, z=pred, group=half),
                   color="black") +
      scale_fill_gradientn(colours=c("gray95", "darkred"), 
                           na.value="white") +
      scale_x_continuous(limits=c(0, 1500), expand=c(0,0)) +
      scale_y_continuous(breaks=c(-.5, 0, .5)) +
      facet_grid(gs~.) +
      labs(x="CWD (mm)",
           y="southness",
           fill="relative\noccupancy") +
      theme_minimal() +
      theme(strip.text=element_blank(),
            panel.grid=element_blank(),
            legend.position=c(.12, .15),
            legend.background=element_rect(fill="white", color=NA))
ggsave("figures/regional_southness_cwd_occupancy_contours.png", p, width=8, height=8, units="in")

