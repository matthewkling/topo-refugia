# Pepperwood tree distribution GAM analysis - projection of model results to future scenarios
# futures based on approximate ratio of 40 mm cwd per degree warming

rm(list=ls())
library(tidyverse)
library(mgcv)

# load 10 m hyperspectral tree layer
hyp <- readRDS('data/pwd_hyp_topo.Rdata')
names(hyp)
dim(hyp)

# sum proportional occurrence for the woody plant columns
hyp$woodyt <- apply(hyp[,15:27],1,sum) 

# most pixels that have woody veg are all woody veg
hist(hyp$woodyt)

# subset to pixels that are at least 50% woody plants and data available for southness
env.complete <- intersect(which(!is.na(hyp$southness)) , which(!is.na(hyp$topoid)))
hypv <- hyp[intersect(which(hyp$woodyt>=0.5) , env.complete),]
dim(hypv)
hypv <- hypv[-which(hypv$rock.group=='ultra'),]
dim(hypv)
names(hypv)

rm('hyp')
hypv <- hypv[,1:29]

rsamp <- sample(nrow(hypv),10000)

# explanation for some of the abiotic variables
names(hypv)
#cwd8110 - 30 year cwd from pepperwood 10 m dem and Flints downscaled model
#ptype - geology
#southness: -cos(asp)*sin(slp)
# tpi100, 500, 1k: Topographic position index with 100, 500, 1000 m radii
# topoid: hydrologic topography index
# model3: Stu Weiss's winter min temp model
# janmin: downscaled janmin from Flint's work


species <- c("Shrubland","Maple", "Buckeye", "Madrone","Tanoak", "Doug.fir","Coast.live.oak", "Blue.oak", "Oregon.oak", "Black.oak", "Valley.oak","Redwood","Bay")
all(species==names(hypv)[15:27])
sci.names <- c('Adenostoma fasciculatum','Acer macrophyllum','Aesculus californica','Arbutus menziesii','Notholithocarpus densiflorus','Pseudotsuga menziesii','Quercus agrifolia','Quercus douglasii','Quercus garryana','Quercus kelloggii','Quercus lobata','Sequoia sempervirens','Umbellularia californica')
cbind(species,sci.names)

# gam model variables
vars <- c("cwd8110", "model3")

# read in gam model outputs
gfits <- readRDS('big_data/pwd_gam2_fits.Rdata')

# cspace is an orthogonal matrix spanning the range of cwd and tmin vals, to visualize model fit
summary(hypv$cwd8110)
summary(hypv$model3)
cvals <- seq(min(hypv$cwd8110,na.rm=T),max(hypv$cwd8110,na.rm=T),length.out = 100)
tvals <- seq(min(hypv$model3,na.rm=T),max(hypv$model3,na.rm=T),length.out = 100)
cspace <- data.frame(matrix(NA,nrow=length(cvals)*length(tvals),ncol=(2+length(species))))
dim(cspace)
names(cspace) <- c('cwd8110','model3',paste0(species,'_pred'))
cspace$cwd8110 <- rep(cvals,100)
cspace$model3 <- rep(tvals,each=100)
head(cspace)

# setup futures
hypF <- hypv
hypF$cwd8110 <- hypv$cwd8110 + 120
hypF$model3 <- hypv$model3 + 3

sfits <- list()
i=1
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  fit <- gfits[[i]]
  hypv[,paste0(sp, "_pred")] <- predict(fit, hypv, type="response")
  hypF[,paste0(sp, "_pred")] <- predict(fit, hypF, type="response")
}

# extract summaries and model fits from gam models
# define model predictors
gam_histfut <- data.frame(species,sci.names,gam.Habund=NA,gam.Fabund=NA,cwd.gam.Hmean=NA,cwd.gam.Fmean=NA,tmin.gam.Hmean=NA,tmin.gam.Fmean=NA,south.gam.Hmean=NA,south.gam.Fmean=NA,topoid.gam.Hmean=NA,topoid.gam.Fmean=NA)

i=1
for (i in 1:length(species)) {
  sp <- species[i]
  message(sp)
  gam_histfut$gam.Habund[i] <- sum(hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$gam.Fabund[i] <- sum(hypF[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$cwd.gam.Hmean[i] <- weighted.mean(hypv$cwd8110,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$tmin.gam.Hmean[i] <- weighted.mean(hypv$model3,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$cwd.gam.Fmean[i] <- weighted.mean(hypF$cwd8110,hypF[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$tmin.gam.Fmean[i] <- weighted.mean(hypF$model3,hypF[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$south.gam.Hmean[i] <- weighted.mean(hypv$southness,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$topoid.gam.Hmean[i] <- weighted.mean(hypv$topoid,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$south.gam.Fmean[i] <- weighted.mean(hypF$southness,hypF[,paste0(sp,"_pred")],na.rm=T)
  gam_histfut$topoid.gam.Fmean[i] <- weighted.mean(hypF$topoid,hypF[,paste0(sp,"_pred")],na.rm=T)
}
gam_histfut$gam.del.abund <- gam_histfut$gam.Fabund - gam_histfut$gam.Habund
gam_histfut$gam.del.south <- gam_histfut$south.gam.Fmean - gam_histfut$south.gam.Hmean
gam_histfut$gam.del.topoid <- gam_histfut$topoid.gam.Fmean - gam_histfut$topoid.gam.Hmean

gam_histfut

plot(gam.Fabund~gam.Habund,data=gam_histfut);abline(0,1)
plot(gam.del.abund~cwd.gam.Hmean,data=gam_histfut)
plot(south.gam.Fmean~south.gam.Hmean,data=gam_histfut);abline(0,1)
plot(topoid.gam.Fmean~topoid.gam.Hmean,data=gam_histfut);abline(0,1)
plot(gam.del.topoid~gam.del.south,data=gam_histfut);abline(h=0);abline(v=0)
write.csv(gam_histfut,'data/gamF3.csv')
