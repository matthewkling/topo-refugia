# Pepperwood tree distribution GAM analysis - fit models and save for further analysis

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
head(hypv)

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

# look at pairwise correlations of some predictors
#pairs(hypv[rsamp,c('cwd8110','southness','TPI100','TPI500','TPI1k','topoid','model3','janmin')])
cor(hypv[rsamp,c('cwd8110','southness','TPI100','TPI500','TPI1k','topoid','model3','janmin')],use = 'pair')

# pairs that should not be combined in a model:
# cwd and southness
# tpi500 and tpi1k
# topoid and model3

# concepts to capture
# exposure - pick either cwd8110 or southness; cwd has advantage that we can change value for futures
# hilltop-valley bottom - used model3 as units are temp and we can explore futures (if the causal factor influencing hilltop-valley bottom distributions is temp; if it's water accumulation, then use topoid or TIP1k and we would keep it fixed for futures)

species <- c("Shrubland","Maple", "Buckeye", "Madrone","Tanoak", "Doug.fir","Coast.live.oak", "Blue.oak", "Oregon.oak", "Black.oak", "Valley.oak","Redwood","Bay")
sci.names <- c('Adenostoma fasciculatum','Acer macrophyllum','Aesculus californica','Arbutus menziesii','Notholithocarpus densiflorus','Pseudotsuga menziesii','Quercus agrifolia','Quercus douglasii','Quercus garryana','Quercus kelloggii','Quercus lobata','Sequoia sempervirens','Umbellularia californica')
cbind(species,sci.names)

# set up orthogonal environmental space spanning range of variables
cr <- seq(min(hypv$cwd8110,na.rm=T),max(hypv$cwd8110,na.rm=T),length.out = 100)
tr <- seq(min(hypv$model3,na.rm=T),max(hypv$model3,na.rm=T),length.out = 100)
cspace <- data.frame(cwd8110=rep(cr,100),model3=rep(tr,each=100),hpX=rep(NA,10000),hpY=rep(NA,10000))
dim(cspace)

#subsample to reduce spatial autocorrelation
rsamp <- sample(nrow(hypv),10000)

# bivariate gams
vars <- c("cwd8110", "model3")

# save models for predictions and projections on future scenarios
gfits <- list()

# construct formula, fit gam, add model predictions to data frame
i=1
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit <- gam(formula, data=hypv[rsamp,], family=binomial(logit))
  #hypv[,paste0(sp, "_pred_nonsp")] <- predict(fit, hypv[rsamp,], type="response")
  #cspace[,paste0(sp, "_pred_nonsp")] <- predict(fit,cspace,type="response")
  gfits[[i]] <- fit
}
saveRDS(gfits,'big_data/pwd_gam2_fits.Rdata')

# construct formula with xy spatial grid, fit gam, add model predictions to data frame
vars <- c("cwd8110", "model3",'hpX','hpY')
i=1
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit <- gam(formula, data=hypv[rsamp,], family=binomial(logit))
  #hypv[,paste0(sp, "_pred_spat")] <- predict(fit, hypv, type="response")
  #cspace[,paste0(sp, "_pred_spat")] <- predict(fit,cspace,type="response")
  gfits[[i]] <- fit
}
saveRDS(gfits,'big_data/pwd_gam2spatial_fits.Rdata')

# Now fit univariate cwd and southness gams and save models
# univariate gams
vars <- c("cwd8110")

# save models for predictions and projections on future scenarios
gfits <- list()

# construct formula, fit gam, add model predictions to data frame
sp=species[1]
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit <- gam(formula, data=hypv, family=binomial(logit))
  hypv[,paste0(sp, "_pred")] <- predict(fit, hypv, type="response")
  cspace[,paste0(sp, "_pred")] <- predict(fit,cspace,type="response")
  gfits[[i]] <- fit
}
saveRDS(gfits,'big_data/pwd_gam1cwd_fits.Rdata')

# read back in to evaluate
gfits <- readRDS('big_data/pwd_gam2_fits.Rdata')
gfits[1]
summary(gfits[[1]])$dev.expl
for (i in 1:length(gfits)) print(summary(gfits[[i]])$dev.expl)
