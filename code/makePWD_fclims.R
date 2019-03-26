# prepare a raster surrounding pepperwood to explore discrete climate change scenarios
library(raster)
library(rgdal)
library(rgeos)
library(maptools)

climate <- stack("big_data/climate/historic.gri")
climate$ppt <- log10(climate$ppt)
plot(climate$ppt)

expand_extent <- function(ex,d) {
  ex@xmin <- ex@xmin - d
  ex@ymin <- ex@ymin - d
  ex@xmax <- ex@xmax + d
  ex@ymax <- ex@ymax + d
  return(ex)
}

pwd <- readShapeSpatial('data/PPshapefile-geo/Pepperwood')
plot(pwd,add=T)
pwdx <- extent(pwd)
pwdx@xmin <- -122.82
pwdx@xmax <- -122.6
pwdx@ymin <- 38.5
pwdx@ymax <- 38.67

pclim <- crop(climate,pwdx)
plot(pclim[[1]])
plot(pwd,add=T)

saveRDS(pwdx,'data/PWDexpanded_extent.Rdata')
saveRDS(pclim,'data/climate/historic/historic.Rdata')

# summary of scenarios in terms of range of changes in climate factors
# cwd: historic = 716, max fut = 831
#     for future exploration: 0, +40, +80, +120
# aet: historic 415, min/max fut = 391/450
#     for future exploration:  -25, 0, 25, 50
# tminmin: historic 4.4, max fut 7.6
#     for future exploration: 0, +1, +2, +3

names(pclim)
# name scenarios A-D for each factor, in three character string
# AAA = cwd 0, aet -25, tminmin 0
# ABA = cwd 0, aet 0, tminmin 0
# etc.

dcwds <- c(0,40,80,120)
daets <- c(-25,0,25,50)
dtminmin <- c(0,1,2,3)
chstring <- c("A","B","C","D")

i=1;j=1;k=1
# rerun this loop to remake folders
if (FALSE) {
  for (i in 1:4) {
    ch1 <- chstring[i]
    for (j in 1:4) {
      ch2 <- chstring[j]
      for (k in 1:4) {
        ch3 <- chstring[k]
        chst <- paste0('clim',ch1,ch2,ch3)
        system(paste0('mkdir data/climate/future/',chst))
      }
    }
  }
}

for (i in 1:4) {
  ch1 <- chstring[i]
  for (j in 1:4) {
    ch2 <- chstring[j]
    for (k in 1:4) {
      ch3 <- chstring[k]
      chst <- paste0('clim',ch1,ch2,ch3)
      fclim <- pclim
      fclim$cwd <- pclim$cwd + dcwds[i]
      fclim$aet <- pclim$aet + daets[j]
      fclim$tminmin <- pclim$tminmin + dtminmin[k]
      saveRDS(fclim,paste0('data/climate/future/',chst,'.Rdata'))
    }
  }
}
