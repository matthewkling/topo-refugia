rm(list=ls())
library(tidyverse)
library(raster)
library(rgdal)
library(data.table)
library(doParallel)
library(rgeos)
library(mgcv)
library(corrplot)
library(dismo)
library(rJava)

select <- dplyr::select

# focal species list
spp <- c("Acer macrophyllum", "Adenostoma fasiculatum",
         "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", 
         "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")

## extract results
dcwds <- c(0,40,80,120)
daets <- c(-25,0,25,50)
dtminmins <- c(0,1,2,3)

infiles <- list.files('data/pwd_distributions/rasters')
length(infiles)/length(spp)
gam_pwd_suit <- data.frame(fname=infiles,sp=NA,var_set=NA,scen=NA,dcwd=NA,daet=NA,dtmintmin=NA,suit=NA)

hypx <- readRDS('data/HYPshapefile-geo/hyp.boundary.Rdata')

i=1
for (i in 1:length(infiles)){
  message(i)
  infile <- infiles[i]
  ras <- raster(paste0('data/pwd_distributions/rasters/',infile))
  chs <- strsplit(infile,'__')
  gam_pwd_suit$sp[i] <- chs[[1]][1]
  gam_pwd_suit$var_set[i] <- chs[[1]][2]
  scen <- substr(chs[[1]][3],1,7)
  gam_pwd_suit$scen[i] <- scen
  gam_pwd_suit$dcwd[i] <- dcwds[match(substr(scen,5,5),c('A','B','C','D'))]
  gam_pwd_suit$daet[i] <- daets[match(substr(scen,6,6),c('A','B','C','D'))]
  gam_pwd_suit$dtmintmin[i] <- dtminmins[match(substr(scen,7,7),c('A','B','C','D'))]
  gam_pwd_suit$suit[i] <- mean(extract(ras,hypx)[[1]])
}
head(gam_pwd_suit)
tail(gam_pwd_suit)

# add column for historical suitability (ABA) to calculate changes in suitability
hsuit <- gam_pwd_suit[gam_pwd_suit$scen=='climABA',]
hsuit
gam_pwd_suit$hsuit <- hsuit$suit[match(gam_pwd_suit$sp,hsuit$sp)]
head(gam_pwd_suit)
tail(gam_pwd_suit)
gam_pwd_suit$dsuit <- gam_pwd_suit$suit - gam_pwd_suit$hsuit
hist(gam_pwd_suit$dsuit)
summary(gam_pwd_suit$dsuit)

# compare historical suitability for maximum across the entire species range
res <- read.csv('data/gam_varseta_modelstats.csv',as.is=T)
res <- merge(res,hsuit)
head(res)
cbind(res$sp,round(res$suit,3),res$max.fit.value,round(res$suit/res$max.fit.value,3))
summary(res$suit)
summary(round(res$suit/res$max.fit.value,3))

# look at futures vs. historical
plot(gam_pwd_suit$suit~gam_pwd_suit$hsuit)
abline(0,1)

# over all combinations, how many are have lower suitability than historic
length(which(gam_pwd_suit$suit < gam_pwd_suit$hsuit))-13
(length(which(gam_pwd_suit$suit < gam_pwd_suit$hsuit))-13)/(832-13)

# look at historical conditions
gam_pwd_suit[gam_pwd_suit$scen=='climABA',]

species <- unique(gam_pwd_suit$sp)
climcf <- data.frame(sp=species,cwd.cf=NA,aet.cf=NA,tminmin.cf=NA,cwdtmin.cf=NA)

i=2
for (i in 1:length(species)){
  message(species[i])
  xx <- gam_pwd_suit[gam_pwd_suit$sp==species[i],]
  fit <- glm(suit~dcwd+daet+dtmintmin,data=xx)
  summary(fit)$coeff
  climcf$sp[i] <- species[i]
  climcf$cwd.cf[i] <- summary(fit)$coeff[2,1] * 120
  climcf$aet.cf[i] <- summary(fit)$coeff[3,1] * 75
  climcf$tminmin.cf[i] <- summary(fit)$coeff[4,1] * 4
  plot(xx$suit[order(xx$daet,xx$dtmintmin,xx$dcwd)])
  xx2 <- xx[xx$scen %in% c('climABA','climBBB','climCBC','climDBD'),]
  fit <- lm(suit~dtmintmin,data=xx2)
  climcf$cwdtmin.cf[i] <- summary(fit)$coeff[2,1]
}
climcf
d <- merge(res[,c('sp','dev.expl','max.fit.value','suit')],climcf)
head(d)
barplot(t(as.matrix(climcf[,-1])),beside=T)

topo <- read.csv('data/pwd_niche_means.csv')
head(topo)
d <- merge(d,topo[,-c(1:2)],by.x='sp',by.y='sci.names')
d

reg <- read.csv('data/regional_niche_stats.csv',as.is=T)
head(reg[reg$var=='cwd',])
d <- merge(d,reg[reg$var=='cwd',],by.x='sp',by.y='species')
head(d)
names(d)
names(d)[15:17] <- c('reg.cwd.mean','reg.cwd.median','reg.cwd.max')
d <- merge(d,reg[reg$var=='aet',],by.x='sp',by.y='species')
names(d)[19:21] <- c('reg.aet.mean','reg.aet.median','reg.aet.max')

d <- merge(d,reg[reg$var=='tminmin',],by.x='sp',by.y='species')
names(d)[23:25] <- c('reg.tminmin.mean','reg.tminmin.median','reg.tminmin.max')
d <- d[,-c(14,18,22)]
d

pairs(d[,c('cwd.mean','reg.cwd.mean','cwd.cf')])
cor(d[,c('cwd.mean','reg.cwd.mean','cwd.cf')])

pairs(d[,c('reg.aet.mean','aet.cf')])
cor(d[,c('reg.aet.mean','aet.cf')])
summary(lm(aet.cf~reg.aet.mean,data=d))

pairs(d[,c('topoid.mean','reg.tminmin.mean','tminmin.cf')])
cor(d[,c('topoid.mean','reg.tminmin.mean','tminmin.cf')])
summary(lm(tminmin.cf~reg.tminmin.mean,data=d))

plot(d$cwd.mean,d$cwd.cf)
summary(lm(d$cwd.cf~d$cwd.mean))
abline(lm(d$cwd.cf~d$cwd.mean))
abline(h=0,lty=2)
sort(d$cwd.cf)

#regional cwd mean
plot(d$reg.cwd.mean,d$cwd.cf)
summary(lm(d$cwd.cf~d$reg.cwd.mean))
abline(lm(d$cwd.cf~d$reg.cwd.mean))
abline(h=0,lty=2)

plot(d$south.mean,d$cwd.cf)
summary(lm(d$cwd.cf~d$south.mean))
abline(lm(d$cwd.cf~d$south.mean))
abline(h=0,lty=2)

# check combined responses to cwd and temp
# not related to landscape position
plot(d$reg.cwd.mean,d$cwdtmin.cf)
summary(lm(d$cwdtmin.cf~d$reg.cwd.mean))
abline(lm(d$cwd.cf~d$reg.cwd.mean))
abline(h=0,lty=2)

plot(d$reg.cwd.mean,d$cwd.mean)
summary(lm(d$cwd.mean~d$reg.cwd.mean))
abline(lm(d$cwd.mean~d$mean))


plot(d$reg.aet.mean,d$aet.cf)
summary(lm(d$aet.cf~d$reg.aet.mean))
abline(lm(d$aet.cf~d$cwd.mean))
abline(h=0,lty=2)
sort(d$aet.cf)

plot(d$mean,d$aet.cf)
summary(lm(d$aet.cf~d$mean))
abline(lm(d$aet.cf~d$mean))
abline(h=0,lty=2)

plot(d$south.mean,d$aet.cf)
summary(lm(d$aet.cf~d$south.mean))
abline(lm(d$aet.cf~d$south.mean))
abline(h=0,lty=2)

plot(d$cwd.mean,d$tminmin.cf)
abline(h=0)

plot(d$south.mean,d$tminmin.cf)
abline(h=0)

plot(d$topoid.mean,d$tminmin.cf)
abline(h=0)

pairs(d[,c('cwd.mean','cwd.cf','aet.cf','tminmin.cf')])
cor(d[,c('cwd.mean','cwd.cf','aet.cf','tminmin.cf')])

### FIGURE 5 for paper
op=par(mfrow=c(1,2))
plot(cwd.cf~mean,data=d,pch=19,xlab='Regional CWD niche mean (mm)',ylab='Response to +120 CWD')
abline(lm(cwd.cf~mean,data=d))
summary(lm(cwd.cf~mean,data=d))
abline(h=0,lty=2)

plot(cwd.cf~cwd.mean,data=d,pch=19,xlab='Topographic CWD niche mean (mm)',ylab='Response to +120 CWD')
abline(lm(cwd.cf~cwd.mean,data=d))
summary(lm(cwd.cf~cwd.mean,data=d))
abline(h=0,lty=2)
par(op)

##### After GAM models have been fit, they can be reloaded here without going through the entire loop below - just to
if (FALSE) {
  for(sp in unique(c(fia$gs, as.character(cch$gs)))){
    outfile <- paste0("data/pwd_distributions/models/",
                      sp, "_", vars, ".rds")
    fit <- readRDS(outfile) 
    message(paste(sp,round(summary(fit)$dev.expl,4)))
  }
}
