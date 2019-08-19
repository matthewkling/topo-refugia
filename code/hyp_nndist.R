rm(list=ls())
library(spatstat)
library(spdep)
# load 10 m hyperspectral tree layer
hyp <- readRDS('big_data/pwd_hyp_topo.Rdata')
names(hyp)
dim(hyp)

range(hyp$hpX)
range(hyp$hpY)

rsamp <- sample(nrow(hyp),500,replace = F)
xyppp <- as.ppp(hyp[rsamp,c('hpX','hpY')],c(range(hyp$hpX),range(hyp$hpY)))
plot(xyppp,axes=T)

xydist <- nndist(xyppp)
hist(xydist)
mean(xydist)
