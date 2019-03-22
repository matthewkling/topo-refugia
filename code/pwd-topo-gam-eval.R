# Pepperwood tree distribution GAM analysis - evaluation of model results

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

rsamp <- sample(nrow(hypv),5000)

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
species==names(hypv)[15:27]
sci.names <- c('Adenostoma fasciculatum','Acer macrophyllum','Aesculus californica','Arbutus menziesii','Notholithocarpus densiflorus','Pseudotsuga menziesii','Quercus agrifolia','Quercus douglasii','Quercus garryana','Quercus kelloggii','Quercus lobata','Sequoia sempervirens','Umbellularia californica')
cbind(species,sci.names)

# gam model variables
vars <- c("cwd8110", "model3")

# read in gam model outputs
gfits <- readRDS('big_data/gam_fits.Rdata')

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

sfits <- list()
i=1
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  fit <- gfits[[i]]
  hypv[,paste0(sp, "_pred")] <- predict(fit, hypv, type="response")
  cspace[,paste0(sp, "_pred")] <- predict(fit,cspace,type="response")
  sfits[[i]] <- summary(fit)
}

# plot gams
names(cspace)
wm <- which(cspace$model3==unique(cspace$model3)[51])
(maxp <- max(cspace[wm,3:15]))
pall <- rainbow(13)
plot(cspace$cwd8110[wm],cspace[wm,3]/max(cspace[wm,3]),type='l',ylim=c(0,1),col=pall[1])
for (i in 4:15) lines(cspace$cwd8110[wm],cspace[wm,i]/max(cspace[wm,i]),col=pall[i-3])

# extract summaries and model fits from gam models
# define model predictors
gam_niche <- data.frame(species,sci.names,gam.tabund=NA,cwd.cspace.opt=NA,tmin.cspace.opt=NA,cwd.hypv.opt=NA,tmin.hypv.opt=NA,cwd.gam.mean=NA,tmin.gam.mean=NA,hypv.pmax=NA,cspace.pmax=NA,dev.expl=NA,cwd.chisq=NA,tmin.chisq=NA)

i=1
for (i in 1:length(species)) {
  sp <- species[i]
  message(sp)
  gam_niche$gam.tabund[i] <- sum(hypv[,paste0(sp,"_pred")],na.rm=T)
  cwm <- which.max(cspace[,paste0(sp,"_pred")])
  pwm <- which.max(hypv[,paste0(sp,"_pred")])
  gam_niche$cwd.cspace.opt[i] <- cspace$cwd8110[cwm]
  gam_niche$tmin.cspace.opt[i] <- cspace$model3[cwm]
  gam_niche$cwd.hypv.opt[i] <- hypv$cwd8110[pwm]
  gam_niche$tmin.hypv.opt[i] <- hypv$model3[pwm]
  gam_niche$cwd.gam.mean[i] <- weighted.mean(hypv$cwd8110,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_niche$tmin.gam.mean[i] <- weighted.mean(hypv$model3,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_niche$hypv.pmax[i] <- hypv[pwm,paste0(sp,"_pred")]
  gam_niche$cspace.pmax[i] <- cspace[cwm,paste0(sp,"_pred")]
  gam_niche$dev.expl[i] <- sfits[[i]]$dev.expl
  gam_niche[i,c('cwd.chisq','tmin.chisq')] <- sfits[[i]]$chi.sq
}
gam_niche
range(gam_niche$dev.expl)

pairs(gam_niche[,c('cwd.cspace.opt','cwd.gam.mean','tmin.cspace.opt','tmin.gam.mean')])
cor(gam_niche[,c('cwd.cspace.opt','cwd.gam.mean','tmin.cspace.opt','tmin.gam.mean')])
write.csv(gam_niche,'data/pwd_gam1.csv')

species
sp=species[12]
message(sp)
plot(cspace[,c('cwd8110',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('cwd8110',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('cwd8110',paste0(sp))])

plot(cspace[,c('model3',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('model3',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('model3',paste0(sp))])

plot(hypv[rsamp,c('southness',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('topoid',paste0(sp,"_pred"))])
plot(hypv[rsamp,c('southness','topoid')])

# Do predictions sum to close to 100%
head(cspace)
cspace$pwoody <- apply(cspace[,grep('pred',names(cspace))],1,sum)
hist(cspace$pwoody)

# plot GAM models on top of each other
midTmin <- median(cspace)



#### LATER
# models with geology
# construct formula, fit gam, add model predictions to data frame
head(hypv)
table(hypv$rock.group)
table(hypv$rock.group.num)
hypv$rock.group3 <- hypv$rock.group
hypv$rock.group3[hypv$rock.group=='ultra'] <- NA

vars3 <- c("cwd8110", "model3", "rock.group3")

# save models
sfit2 <- list()

# construct formula, fit gam, add model predictions to data frame
i=11
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula2 <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit2 <- gam(formula2, data=hypv, family=binomial(logit))
  formula3 <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + "),' + ','rock.group3'))
  fit3 <- gam(formula3, data=hypv, family=binomial(logit))
  hypv[,paste0(sp, "_pred3")] <- predict(fit3, hypv, type="response")
  
  # fit cwd model on different rock types
  # formula1 <- as.formula(paste0(sp, " ~ ", "s(cwd8110)"))
  # fit.mun <- gam(formula1, data=hypv[hypv$rock.group %in% c('melange','uncon'),], ,family=binomial(logit))
  # hypv[,paste0(sp, "_pred.mun")] <- predict(fit.mun, hypv, type="response")
  # 
  # fit.bl <- gam(formula1, data=hypv[hypv$rock.group %in% c('block'),],family=binomial(logit))
  # hypv[,paste0(sp, "_pred.bl")] <- predict(fit.bl, hypv, type="response")
  
}
summary(fit3)
AIC(fit3)
vis.gam(fit3, c("cwd8110","model3"), type="response")
boxplot(hypv[,paste0(sp, "_pred3")]~hypv$rock.group,)
plot(hypv[rsamp,c('cwd8110',paste0(sp, "_pred.mun"))])
points(hypv[rsamp,c('cwd8110',paste0(sp, "_pred.bl"))])

