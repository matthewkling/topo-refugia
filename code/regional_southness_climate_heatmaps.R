

library(tidyverse)
select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

# climate data, clipped to conterminous US
# (created in regional niche model script)
clim <- raster("big_data/climate/cwd_masked.tif")

# load FULL fia dataset, summarize to subplot-species level
d <- readRDS("e:/fia/topoclimate/data_munged.rds") %>%
      select(-c(bio1:bio19)) %>%
      mutate(southness = -northness,
             gs = paste(genus, species),
             gs = sub("Lithocarpus densiflorus", 
                      "Notholithocarpus densiflorus", gs))


# climate at each plot
coordinates(d) <- c("lon", "lat")
d$cwd <- extract(clim, d)
d <- as.data.frame(d)


# binned frequencies

b <- d %>%
      filter(lon < -95) %>%
      mutate(cwd_bin = plyr::round_any(cwd, 50),
             southness_bin = plyr::round_any(southness, .1)) %>%
      expand(cwd_bin, southness_bin, gs)

b <- d %>%
      filter(lon < -95) %>%
      mutate(cwd_bin = plyr::round_any(cwd, 50),
             southness_bin = plyr::round_any(southness, .1)) %>%
      left_join(b, .) %>% ungroup() %>%
      group_by(cwd_bin, southness_bin) %>%
      mutate(n_plots = length(unique(na.omit(plot)))) %>%
      filter(gs %in% spp) %>%
      group_by(cwd_bin, southness_bin, gs) %>%
      summarize(n_trees = length(na.omit(southness)),
                n_plots = mean(n_plots, na.rm=T)) %>%
      mutate(occ = n_trees/n_plots,
             occ = ifelse(!is.finite(occ), 0, occ))

w2d <- c("Pseudotsuga menziesii", "Acer macrophyllum", "Notholithocarpus densiflorus", 
         "Sequoia sempervirens",  "Arbutus menziesii", "Umbellularia californica", 
         "Quercus garryana", "Quercus kelloggii", "Quercus agrifolia",
         "Quercus lobata", "Aesculus californica","Quercus douglasii")

b <- b %>% filter(abs(southness_bin)<=.5) %>%
      mutate(gs=factor(gs, levels=w2d)) %>% # arrange wet to dry
      group_by(gs) %>%
      mutate(occ_norm = scales::rescale(log10(occ)),
             occ_norm = ifelse(!is.finite(occ_norm), 0, occ_norm))

txt <- b %>%
      select(gs) %>%
      distinct()

p <- ggplot() +
      geom_tile(data=b, 
                aes(cwd_bin, southness_bin, fill=occ_norm)) +
      geom_text(data=txt, aes(x=1450, y=0, label=sub(" ", "\n", gs)),
                color="black", hjust=1, lineheight=.75) +
      scale_fill_gradientn(colours=c("gray80", "forestgreen", "black")) +
      facet_grid(gs~.) +
      scale_y_continuous(breaks=c(-.5, 0, .5)) +
      xlim(0, 1500) +
      labs(x="CWD (mm)",
           y="southness",
           fill="log trees\nper plot") +
      theme_minimal() +
      theme(strip.text=element_blank())
ggsave("figures/regional_southness_cwd_occupancy.png", p, width=8, height=8, units="in")

