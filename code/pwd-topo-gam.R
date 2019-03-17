# Pepperwood tree distribution GAM analysis

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
hypv <- hyp[intersect(which(hyp$woodyt>=0.5) , which(!is.na(hyp$southness))),]
dim(hypv)

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
pairs(hypv[rsamp,c('cwd8110','southness','TPI100','TPI500','TPI1k','topoid','model3','janmin')])
cor(hypv[rsamp,c('cwd8110','southness','TPI100','TPI500','TPI1k','topoid','model3','janmin')],use = 'pair')

# pairs that should not be combined in a model:
# cwd and southness
# tpi500 and tpi1k
# topoid and model3

# concepts to capture
# exposure - pick either cwd8110 or southness; cwd has advantage that we can change value for futures
# hilltop-valley bottom - used model3 as units are temp and we can explore futures (if the causal factor influencing hilltop-valley bottom distributions is temp; if it's water accumulation, then use topoid or TIP1k and we would keep it fixed for futures)

species <- c("Shrubland","Maple", "Buckeye", "Madrone","Bay", 
             "Tanoak", "Doug.fir","Redwood",
             "Coast.live.oak", "Blue.oak", "Oregon.oak", "Black.oak", "Valley.oak")
vars <- c("cwd8110", "model3")

# save model summaries
sfits <- list()

# cspace is an orthogonal matrix spanning the range of cwd and tmin vals, to visualize model fit
cvals <- seq(min(hypv$cwd8110),max(hypv$cwd8110),length.out = 100)
tvals <- seq(min(hypv$model3,na.rm=T),max(hypv$model3,na.rm=T),length.out = 100)
cspace <- data.frame(matrix(NA,nrow=length(cvals)*length(tvals),ncol=(2+length(species))))
dim(cspace)
names(cspace) <- c('cwd8110','model3',paste0(species,'_pred'))
cspace$cwd8110 <- rep(cvals,100)
cspace$model3 <- rep(tvals,each=100)
head(cspace)

# construct formula, fit gam, add model predictions to data frame
sp=species[1]
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit <- gam(formula, data=hypv, family=binomial(logit))
  hypv[,paste0(sp, "_pred")] <- predict(fit, hypv, type="response")
  cspace[,paste0(sp, "_pred")] <- predict(fit,cspace,type="response")
  sfits[[i]] <- summary(fit)
}

# extract summaries and model fits from gam models
# define model predictors
fit_max <- data.frame(species,cwd_opt=NA,tmin_opt=NA,dev.expl=NA)

i=1
for (i in 1:length(species)) {
  sp <- species[i]
  message(sp)
  wm <- which.max(cspace[,paste0(sp,"_pred")])
  fit_max[fit_max$species==sp,'cwd_opt'] <- cspace$cwd8110[wm]
  fit_max[fit_max$species==sp,'tmin_opt'] <- cspace$model3[wm]
  fit_max[fit_max$species==sp,'dev.expl'] <- summary(fits[[i]])$dev.expl
}

fit_max[order(fit_max$cwd_opt),]
fit_max[order(fit_max$tmin_opt),]
write.csv(fit_max,'data/pwd_gam1.csv')
sp=species[11]

plot(cspace[,c('cwd8110',paste0(sp,"_pred"))])
plot(cspace[,c('model3',paste0(sp,"_pred"))])


#### LATER
# models with geology
# construct formula, fit gam, add model predictions to data frame
head(hypv)
table(hypv$rock.group.num)
vars <- c("cwd8110", "model3", "rock.group")

# save models
sfit2 <- list()

# construct formula, fit gam, add model predictions to data frame
i=11
for(i in 1:length(species)){
  sp <- species[i]
  message(sp)
  formula <- as.formula(paste0(sp, " ~ ", paste0("s(", vars, ")", collapse=" + ")))
  fit <- gam(formula, data=hypv, family=binomial(logit))
}
summary(fit)