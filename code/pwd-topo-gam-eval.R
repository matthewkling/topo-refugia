# Pepperwood tree distribution GAM analysis - evaluation of model results

rm(list=ls())
library(tidyverse)
library(mgcv)

# load 10 m hyperspectral tree layer
hyp <- readRDS('big_data/pwd_hyp_topo.Rdata')
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

# choose random rows before selecting spatial vs non-spatial
rsamp <- sample(nrow(hypv),10000)

# gam model variables
#vars <- c("cwd8110", "model3")
SPATIAL <- FALSE

# cspace is an orthogonal matrix spanning the range of cwd and tmin vals, to visualize model fit
summary(hypv$cwd8110)
summary(hypv$model3)
cvals <- seq(min(hypv$cwd8110,na.rm=T),max(hypv$cwd8110,na.rm=T),length.out = 100)
tvals <- seq(min(hypv$model3,na.rm=T),max(hypv$model3,na.rm=T),length.out = 100)

(hpXmean <- mean(hypv$hpX))
(hpYmean <- mean(hypv$hpY))

if (SPATIAL) cspace <- data.frame(cwd8110=rep(cvals,100),model3=rep(tvals,each=100),hpX=rep(hpXmean,10000),hpY=rep(hpYmean,10000)) else cspace <- data.frame(cwd8110=rep(cvals,100),model3=rep(tvals,each=100))

dim(cspace)
head(cspace)

# extract summaries and model fits from gam models
# define model predictors
# '2' is for cwd,model3 bivariate gam; '1' is for cwd univariate gam
gam_niche <- data.frame(species,sci.names,gam2.tabund=NA,cwd2.cspace.opt=NA,tmin2.cspace.opt=NA,cwd2.hypv.opt=NA,tmin2.hypv.opt=NA,cwd2.gam.mean=NA,tmin2.gam.mean=NA,hypv2.pmax=NA,cspace2.pmax=NA,dev2.expl=NA,cwd2.chisq=NA,tmin2.chisq=NA,cwd1.cspace.opt=NA,cwd1.hypv.opt=NA,cwd1.gam.mean=NA,hypv1.pmax=NA,cspace1.pmax=NA,cwd1.dev.expl=NA)

# read in gam bivariate model outputs
if (SPATIAL) gfits <- readRDS('big_data/pwd_gam2spatial_fits.Rdata') else gfits <- readRDS('big_data/pwd_gam2_fits.Rdata')
#str(gfits)

sfits <- list()
i=1
for (i in 1:length(species)){
  sp <- species[i]
  message(sp)
  fit <- gfits[[i]]
  hypv[rsamp,paste0(sp, "_pred")] <- predict(fit, hypv[rsamp,], type="response")
  cspace[,paste0(sp, "_pred")] <- predict(fit,cspace,type="response")
  sfits[[i]] <- summary(fit)
}

head(cspace)

# plot gams
names(cspace)
wm <- which(cspace$model3==unique(cspace$model3)[51])
if (SPATIAL) minPcol <- 5 else minPcol <- 3
maxPcol <- minPcol+12
maxp <- max(cspace[wm,minPcol:maxPcol])
pall <- rainbow(13)
plot(cspace$cwd8110[wm],cspace[wm,minPcol]/max(cspace[wm,minPcol]),type='l',ylim=c(0,1),col=pall[1])
for (i in (minPcol+1):maxPcol) lines(cspace$cwd8110[wm],cspace[wm,i]/max(cspace[wm,i]),col=pall[i-3])

i=1
for (i in 1:length(species)) {
  sp <- species[i]
  message(sp)
  gam_niche$gam2.tabund[i] <- sum(hypv[,paste0(sp,"_pred")],na.rm=T)
  cwm <- which.max(cspace[,paste0(sp,"_pred")])
  pwm <- which.max(hypv[,paste0(sp,"_pred")])
  gam_niche$cwd2.cspace.opt[i] <- cspace$cwd8110[cwm]
  gam_niche$tmin2.cspace.opt[i] <- cspace$model3[cwm]
  gam_niche$cwd2.hypv.opt[i] <- hypv$cwd8110[pwm]
  gam_niche$tmin2.hypv.opt[i] <- hypv$model3[pwm]
  gam_niche$cwd2.gam.mean[i] <- weighted.mean(hypv$cwd8110,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_niche$tmin2.gam.mean[i] <- weighted.mean(hypv$model3,hypv[,paste0(sp,"_pred")],na.rm=T)
  gam_niche$hypv2.pmax[i] <- hypv[pwm,paste0(sp,"_pred")]
  gam_niche$cspace2.pmax[i] <- cspace[cwm,paste0(sp,"_pred")]
  gam_niche$dev2.expl[i] <- sfits[[i]]$dev.expl
  gam_niche[i,c('cwd2.chisq','tmin2.chisq')] <- sfits[[i]]$chi.sq
}
gam_niche
gam_niche[,c('species','cwd2.cspace.opt')]

if (SPATIAL) write.csv(gam_niche,'data/pwd_gam_Aug19_spatial.csv') else write.csv(gam_niche,'data/pwd_gam_Aug19_nonspatial.csv')

#### COMPARE SPATIAL and NON-SPATIAL outputs
gam_spatial <- read.csv('data/pwd_gam_Aug19_spatial.csv')
gam_nonsp <- read.csv('data/pwd_gam_Aug19_nonspatial.csv')

plot(gam_nonsp$cwd2.cspace.opt,gam_spatial$cwd2.cspace.opt)
cbind(species,gam_nonsp$cwd2.cspace.opt,gam_spatial$cwd2.cspace.opt)

###### END HERE FOR NOW
# read in gam univariate model outputs
gfits <- readRDS('big_data/pwd_gam1cwd_fits.Rdata')

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
unique(cspace$model3)

wm <- which(cspace$model3==unique(cspace$model3)[51])
(maxp <- max(cspace[,3:15]))
pall <- rainbow(13)
plot(cspace$cwd8110[wm],cspace[wm,3]/max(cspace[wm,3]),type='l',ylim=c(0,1),col=pall[1])
for (i in 4:15) lines(cspace$cwd8110[wm],cspace[wm,i]/max(cspace[wm,i]),col=pall[i-3])


i=1
for (i in 1:length(species)) {
  sp <- species[i]
  message(sp)
  
  cwm <- which.max(cspace[,paste0(sp,"_pred")])
  pwm <- which.max(hypv[,paste0(sp,"_pred")])
  gam_niche$cwd1.cspace.opt[i] <- cspace$cwd8110[cwm]
  gam_niche$cwd1.hypv.opt[i] <- hypv$cwd8110[pwm]
  
  gam_niche$cwd1.gam.mean[i] <- weighted.mean(hypv$cwd8110,hypv[,paste0(sp,"_pred")],na.rm=T)
  
  gam_niche$hypv1.pmax[i] <- hypv[pwm,paste0(sp,"_pred")]
  gam_niche$cspace1.pmax[i] <- cspace[cwm,paste0(sp,"_pred")]
  gam_niche$cwd1.dev.expl[i] <- sfits[[i]]$dev.expl
}
gam_niche

plot(gam_niche$cwd1.dev.expl,gam_niche$dev2.expl);abline(0,1)
cbind(gam_niche$cwd1.dev.expl,gam_niche$dev2.expl)


pairs(gam_niche[,c('cwd2.cspace.opt','cwd2.gam.mean','tmin2.cspace.opt','tmin2.gam.mean')])
cor(gam_niche[,c('cwd2.cspace.opt','cwd2.gam.mean','tmin2.cspace.opt','tmin2.gam.mean')])

pairs(gam_niche[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd2.gam.mean','cwd1.gam.mean')])
cor(gam_niche[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd2.gam.mean','cwd1.gam.mean')])

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

#### mantel tests for spatial structure of topoclimate
head(hypv)
rsamp <- sample(nrow(hypv),100)
#rsamp <- c(rsamp,rsamp+1)
edist <- function(x1,y1,x2,y2) sqrt((x1-x1)^2+(y2-y1)^2)

dmatrix <- c()
i=2;j=1
yvar <- hypv$cwd8110
plot(hypv$hpX[rsamp],hypv$hpY[rsamp])
for (i in 2:length(rsamp)) {
  ri <- rsamp[i]
  for (j in 1:(i-1)) {
    rj <- rsamp[j]
    dist <- edist(hypv$hpX[ri],hypv$hpY[ri],
                  hypv$hpX[rj],hypv$hpY[rj])
    dy <- abs(yvar[ri]-yvar[rj])
    dmatrix <- rbind(dmatrix,c(dist,dy))
  }
}
head(dmatrix)
plot(dmatrix,log='x')
cor(dmatrix)
