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





species <- c("Maple", "Buckeye", "Madrone", "Tanoak", "Doug.fir", "Redwood",
             "Coast.live.oak", "Blue.oak", "Oregon.oak", "Black.oak", "Valley.oak")

# define model predictors
vars <- c("cwd8110", "model3")


for(sp in species){
      message(sp)
      
      # construct formula, fit gam, add model predictions to data frame
      spname <- paste0(sp, "_X")
      formula <- as.formula(paste0(sp, "_X ~ ", paste0("s(", vars, ")", collapse=" + ")))
      fit <- gam(formula, data=hypv, family=binomial(logit))
      hypv[,paste0(sp, "_pred")] <- predict(fit, hypv)
      
}
