
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
            filter(abs(southness) <= .5,
                   !is.na(southness),
                   !is.na(cwd))
      
      for(half in c("high", "low")){
            for(quant in c(.5, .6, .7, .8, .85, .9, .95)){
                  
                  if(half == "high") lim <- max(sd$cwd[sd$pres], na.rm=T) + 100
                  if(half == "low") lim <- min(sd$cwd[sd$pres], na.rm=T) - 100
                  
                  if(half == "high") mid <- quantile(sd$cwd[sd$pres], quant, na.rm=T)
                  if(half == "low") mid <- quantile(sd$cwd[sd$pres], 1-quant, na.rm=T)
                  
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
                              quantile=quant,
                              edge = half,
                              mid = mid,
                              lim = lim,
                              n_pres = sum(md$pres),
                              n_subplots = nrow(md),
                              slope_southness = coef(fit)["southness"],
                              slope_cwd = coef(fit)["cwd"],
                              slope = - coef(fit)["cwd"] / coef(fit)["southness"],
                              p_southness = summary(fit)$coefficients["southness", "Pr(>|z|)"],
                              p_cwd = summary(fit)$coefficients["cwd", "Pr(>|z|)"]
                        )
                  slopes <- rbind(slopes, pred)
            }
      }
}




q <- .9

lines <- slopes %>% 
      filter(quantile==q) %>% 
      select(gs, edge, slope, mid, lim, p_southness) %>% 
      distinct() %>% 
      mutate(gs=factor(gs, levels=w2d),
             intercept=-slope*(mid+(lim-mid)/2),
             sig=p_southness < 0.05) %>%
      mutate(bezel_end_upper = case_when(edge=="high" & ((.5 - intercept)/slope) < mid ~ mid,
                                         edge=="high" & ((.5 - intercept)/slope) > mid ~ ((.5 - intercept)/slope),
                                         edge=="low" & ((.5 - intercept)/slope) > mid ~ mid,
                                         edge=="low" & ((.5 - intercept)/slope) < mid ~ ((.5 - intercept)/slope)),
             bezel_end_lower = case_when(edge=="high" & ((-.5 - intercept)/slope) > lim ~ lim,
                                         edge=="high" & ((-.5 - intercept)/slope) < lim ~ ((-.5 - intercept)/slope),
                                         edge=="low" & ((-.5 - intercept)/slope) < lim ~ lim,
                                         edge=="low" & ((-.5 - intercept)/slope) > lim ~ ((-.5 - intercept)/slope))
      )


p <- ggplot() +
      
      # heatmaps
      geom_tile(data= b %>% 
                      group_by(gs) %>% mutate(occ = occ/max(occ)) %>% 
                      ungroup() %>% mutate(gs=factor(gs, levels=w2d)), 
                aes(cwd_bin, southness_bin, fill=occ)) +
      
      # bounding boxes of analysis regions
      geom_rect(data=lines, aes(xmin=mid, xmax=lim, ymin=-.5, ymax=.5),
                color="gray80", fill=NA) +
  
      # species names
      geom_text(data=txt %>% mutate(gs=factor(gs, levels=w2d)), 
                aes(x=1450, y=0, label=sub(" ", "\n", gs)),
                color="black", hjust=1, lineheight=.75, fontface="italic") +
      
      # sloped isoclines
      geom_segment(data=lines, 
                   aes(x=mid, xend=lim, y=slope*mid+intercept, yend=slope*lim+intercept,
                       linetype = sig),
                   color="black") +
      
      # upper lines connecting isocline end to data boundary
      geom_segment(data=lines,
                   aes(x=mid, xend=bezel_end_upper, y=.5, yend=.5,
                       linetype = sig),
                   color="black") +
      
      # lower lines connecting isocline end to data boundary
      geom_segment(data=lines, 
                   aes(x=mid, xend=bezel_end_lower, y=-.5, yend=-.5,
                       linetype = sig),
                   color="black") +
      
      # styling
      scale_linetype_manual(values=c(2, 1), guide=F) +
      scale_fill_gradientn(colours=c("gray95", "orange", "red", "darkred"), 
                           na.value="white") +
      coord_cartesian(ylim=c(-.5,.5), xlim=c(0, 1550), expand=0) +
      scale_y_continuous(breaks=c(-.3, 0, .3)) +
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

# export slope data
slopes %>%
      select(-cwd, -southness, -pred) %>%
      distinct() %>%
      write_csv("data/logistic_regression_coefficients.csv")



slopes <- read_csv("data/logistic_regression_coefficients.csv") %>%
      filter(quantile==0.9) %>%
      select(gs, edge, slope_southness:p_cwd) %>%
      distinct() %>%
      mutate(ratio=-1/slope)


slopes %>%
      mutate(edge = paste(edge, "CWD")) %>%
      rename(species = gs,
             niche_edge = edge,
             southness_coef = slope_southness,
             CWD_coef = slope_cwd,
             southness_pvalue = p_southness,
             CWD_pvalue = p_cwd) %>%
      select(species, niche_edge, CWD_coef, CWD_pvalue,
             southness_coef, southness_pvalue) %>%
      mutate_at(vars(CWD_coef:southness_pvalue), signif, digits=3) %>%
      write_csv("data/table_s1.csv")





table(slopes$edge, sign(slopes$slope_cwd))
table(slopes$edge, sign(slopes$slope_southness))
table(sign(slopes$slope_southness)==sign(slopes$slope_cwd))

table(sign(slopes$slope_southness)==sign(slopes$slope_cwd), slopes$p_southness < 0.05)

table(slopes$p_cwd < 0.05)

table(slopes$p_cwd < 0.05, slopes$p_southness < 0.05)

slopes %>% filter(p_cwd > 0.05)


mean(slopes$ratio)
median(slopes$ratio)
summary(slopes$ratio[slopes$p_southness < 0.05])

hist(slopes$slope[slopes$p_southness < 0.05])
summary(slopes$slope[slopes$p_southness < 0.05])
summary(slopes$slope)


# no sig diff in slopes between wet and dry edges
t.test(slopes$ratio[slopes$edge=="high"],
       slopes$ratio[slopes$edge=="low"], paired=T)

ss <- slopes %>%
      group_by(gs) %>%
      mutate(sig = all(p_southness < 0.05)) %>%
      filter(sig)

ggplot(ss, aes(edge, slope, group=gs)) + geom_line()

t.test(ss$slope[ss$edge=="high"],
       ss$slope[ss$edge=="low"], paired=T)


################## every species #################


# climate data, clipped to conterminous US
# (created in regional niche model script)
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
      filter(abs(southness_bin)<=.5)



# climate at each plot
sites <- d %>% 
      select(site, lon, lat, southness) %>%
      distinct()
coordinates(sites) <- c("lon", "lat")
sites$cwd <- extract(clim$cwd, sites)
sites$tmax <- extract(clim$tmaxmax, sites)
sites <- as.data.frame(sites)



for(s in unique(b$gs)){
      message(s)
      
      sd <- d %>%
            filter(gs == s) %>%
            left_join(sites, .) %>%
            mutate(pres = !is.na(gs)) %>%
            filter(abs(southness) <= .5,
                   !is.na(southness),
                   !is.na(cwd))
      if(sum(sd$pres) < 100) next()
      
      
      slopes <- data.frame()
      for(half in c("high", "low")){
            quant <- .9
            
            if(half == "high") lim <- max(sd$cwd[sd$pres], na.rm=T) + 100
            if(half == "low") lim <- min(sd$cwd[sd$pres], na.rm=T) - 100
            
            if(half == "high") mid <- quantile(sd$cwd[sd$pres], quant, na.rm=T)
            if(half == "low") mid <- quantile(sd$cwd[sd$pres], 1-quant, na.rm=T)
            
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
                        quantile=quant,
                        edge = half,
                        mid = mid,
                        lim = lim,
                        n_pres = sum(md$pres),
                        n_subplots = nrow(md),
                        slope_southness = coef(fit)["southness"],
                        slope_cwd = coef(fit)["cwd"],
                        slope = - coef(fit)["cwd"] / coef(fit)["southness"],
                        p_southness = summary(fit)$coefficients["southness", "Pr(>|z|)"],
                        p_cwd = summary(fit)$coefficients["cwd", "Pr(>|z|)"]
                  )
            slopes <- rbind(slopes, pred)
      }
      
      lines <- slopes %>% 
            select(gs, edge, slope, mid, lim) %>% 
            distinct() %>% 
            mutate(intercept=-slope*(mid+(lim-mid)/2))
      
      p <- ggplot() +
            geom_tile(data= b %>% filter(gs==s) %>% mutate(occ = occ/max(occ)), 
                      aes(cwd_bin, southness_bin, fill=occ)) +
            geom_segment(data=lines, 
                         aes(x=mid, xend=lim, y=slope*mid+intercept, yend=slope*lim+intercept),
                         color="black") +
            geom_segment(data=lines, 
                         aes(x=mid, xend=(.5 - intercept)/slope, y=.5, yend=.5),
                         color="black") +
            geom_segment(data=lines, 
                         aes(x=mid, xend=(-.5 - intercept)/slope, y=-.5, yend=-.5),
                         color="black") +
            geom_segment(data=lines, 
                         aes(x=mid, xend=(-.5 - intercept)/slope, y=-.5, yend=-.5),
                         color="black") +
            geom_vline(data=slopes, 
                       aes(xintercept=mid),
                       color="black", linetype=3) +
            geom_vline(data=slopes, 
                       aes(xintercept=lim),
                       color="black", linetype=3) +
            scale_fill_gradientn(colours=c("gray95", "orange", "red", "darkred"), 
                                 na.value="white") +
            coord_cartesian(ylim=c(-.5,.5), xlim=c(0, 1500), expand=0) +
            scale_y_continuous(breaks=c(-.3, 0, .3)) +
            labs(x="CWD (mm)",
                 y="southness",
                 fill="relative\noccupancy    ",
                 title=s) +
            theme_minimal() +
            theme(panel.grid=element_blank(),
                  legend.position="top")
      
      ggsave(paste0("figures/southness_heatmaps/", s, ".png"), p, width=8, height=3, units="in")
      
}







