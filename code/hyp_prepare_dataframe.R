## assemble data.frame of CWD, geology, and tree ID for 10 m aea layers

rm(list=ls())
library(raster)
library(sp)
library(maptools)
library(nnet)
library(sf)
library(rgdal)
source('/Users/david/Documents/Projects/Toolbox/spatial_tools/expand_extent.R')

#geology
mc <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/geology/McLaughlinMap_ta/PepperwoodGeologyQuadRaster.grd')
plot(mc)
mc_lut <- read.csv('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/geology/McLaughlinMap_ta/PepperwoodGeologyQuadRaster_LUT.csv',row.names=1)
mc_lut

#topography
dem <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/pwd_base/pw_10m_t2.asc')
projection(dem) <- CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000')
slp.rad <- terrain(dem,'slope')
asp.rad <- terrain(dem,'aspect')
southness <- -cos(asp.rad)*sin(slp.rad)
plot(dem)
plot(southness)
slp.deg <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/slpdegree.asc')
marrad <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/marrad.asc')
tpi100 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_100.asc')
tpi500 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_500.asc')
tpi1k <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_1000.asc')
plp100 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/PLP_100.asc')
plp500 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/PLP_500.asc')
tpi1k <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tpi_1000.asc')
topoid <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/topoidx.asc')
model3 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/model3.asc')
janmin <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/PW10m_topoclimate/tmnavejan.asc')

######
#compare extent problem - need to find problem layer
pwd_topo <- stack(dem,slp.deg,southness,marrad,tpi100,tpi500,tpi1k,topoid,model3,janmin)
names(pwd_topo) <- c('elevation','slope.deg','southness','marchRad','TPI100','TPI500','TPI1k','topoid','model3','janmin')
#pwd_topo <- stack(dem,slp,marrad,tpi100,tpi500,tpi1k,plp100,plp500,tpi1k,topoid,model3,janmin)
#######

#cwd
cwd <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/bcm_outputs/aea/cwd1981_2010_ave.asc')
cwd[cwd==-9999] <- NA
plot(cwd)
xx <- getValues(cwd)
length(xx)
nax <- rep(NA,length(xx))
cwdna <- setValues(cwd,nax)
cwdna

#tree map
hp <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_spectral_metrics_class_teale-albers_10m.tif')
hp
hp.ind <- stack(raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_10m.tif',band=1))
for (i in 2:22) hp.ind <- addLayer(hp.ind,raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_10m.tif',band=i))
(rCodes <- read.csv('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/utm10-original/svm_raster_codes.csv',as.is=T))
names(hp.ind) <- rCodes$ShortName
plot(hp.ind)

#pwd shapefile
pwd <- readShapePoly('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/PPshapefiles/PPshapefile-teale-albers/Pepperwood')
plot(pwd,axes=T)

#ogrDrivers()
#dsn <- system.file("/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/PPshapefiles/PPshapefile-teale-albers/Pepperwood.shp",package='rgdal')[1]
#dsn
#ogrListLayers(dsn)
#pwd <- rgdal::readOGR('/Users/david/Documents/Projects/Pepperwood/PWD_GIS/PPshapefiles/PPshapefile-teale-albers/')

pwdx <- extent(pwd)
pwdx <- expand_extent(pwdx,3500)
plot(crop(cwd,pwdx))
plot(hp,add=T)
plot(pwd,add=T)
plot(mc,add=T)
plot(pwd,add=T)

extentIntersect <- function(x1,x2) {
  x3 <- x1
  x3@xmin <- max(x1@xmin,x2@xmin)
  x3@xmax <- min(x1@xmax,x2@xmax)
  x3@ymin <- max(x1@ymin,x2@ymin)
  x3@ymax <- min(x1@ymax,x2@ymax)
  return(x3)
}

#expand geology to merge all together
mc.ex <- merge(mc,cwdna)
plot(mc.ex)
plot(hp,add=T)
plot(cwd,add=T)

hpx <- extent(hp)
#mcx <- extent(pwdx)
#overx <- extentIntersect(hpx,mcx)

cwdc <- crop(cwd,hpx)
hpc <- crop(hp,hpx)
mcc <- crop(mc.ex,hpx)
hpic <- crop(hp.ind,hpx)
pwd_topoc <- crop(pwd_topo,hpx)

plot(cwdc)
plot(hpc,add=T)
plot(mcc,add=T)
plot(hpic$Coast.live.oak,add=T)
plot(pwd,add=T)
#plot(hpic)
#plot(pwd_topoc)

plot(pwd_topoc[[1]])
plot(hpic$Coast.live.oak,add=T)

if (FALSE) {
  hyp.boundary <- drawPoly()
  projection(hyp.boundary) <- projection(hpic)
  
  getwd()
  saveRDS(hyp.boundary,'/Users/david/Google Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/data/HYPshapefile-teale-albers/hyp.boundary.Rdata')
} else {
  hyp.boundary <- readRDS('/Users/david/Google Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/data/HYPshapefile-teale-albers/hyp.boundary.Rdata')
}
plot(hyp.boundary,add=T)

hpXY <- xyFromCell(hpc,1:length(hpc))
dim(hpXY)
head(hpXY)

extent(cwdc)
extent(hpc)
extent(mcc)
extent(hpic)
extent(pwd_topoc)

pwds <- stack(cwdc,mcc,pwd_topoc,hpc)
pwds
names(pwd_topoc)
names(pwds) <- c('cwd8110','ptype',names(pwd_topoc),'hpspp')
#plot(pwds)
names(hpic)
pwds <- addLayer(pwds,hpic[[c(3,5,11:22)]])
names(pwds)
plot(pwds)

pwdv <- data.frame(getValues(pwds))
(names(pwdv) <- names(pwds))
dim(pwdv)
head(pwdv)

# NOW recode rocks as hard/soft (?) and check for differences

# unconsolidated: Qal,Qls,QTg,Qls?,QTg?
# uncon: 14,22,23,28,29
# melange: fcm, fcm? = 4,5
# ultramafic: Jos, Jos? = 10,11
# block: 'fcs2','fcv','Tsb','TSr','Tst','fcs2?','Tsb?','Tst?' = 6,7,8,37,38,42,44,45

mc_lut
pwdv$rock.group <- NA
pwdv$rock.group[pwdv$ptype%in%c(14,22,23,28,29)] <- 'uncon'
pwdv$rock.group[pwdv$ptype%in%c(4,5)] <- 'melange'
pwdv$rock.group[pwdv$ptype%in%c(10,11)] <- 'ultra'
pwdv$rock.group[pwdv$ptype%in%c(6,7,8,37,38,42,44,45)] <- 'block'
table(pwdv$rock.group)

pwdv$rock.group.num <- NA
pwdv$rock.group.num <- c(1:4)[match(pwdv$rock.group,sort(unique(pwdv$rock.group)))]
table(pwdv$rock.group.num)
rock.type <- setValues(hpc,pwdv$rock.group.num)

plot(rock.type,col=c('black','blue','green','brown'))
plot(pwd,add=T,border='white')

# look at distribution of cwd on different rock.types
range(pwdv$cwd8110,na.rm=T)
op=par(mfrow=c(3,1))
hist(pwdv$cwd8110[pwdv$rock.group=='block'],breaks=seq(200,1400,by=50))
hist(pwdv$cwd8110[pwdv$rock.group=='melange'],breaks=seq(200,1400,by=50))
hist(pwdv$cwd8110[pwdv$rock.group=='uncon'],breaks=seq(200,1400,by=50))
par(op)

mc_lut

# make dummy vars for selected veg types
# 2017 Map Key
# 1 Builtup
# 2 Forest (does not exist; replaced by species identification 11-22)
# 3 Grassland
# 4 Orchard
# 5 Shrubland
# 6 Vineyard
# 7 Water
# 8 WetHerbaceous
# 11 Big-leaf maple (Acer macrophylum)
# 12 California buckeye (Aesculus californica)
# 13 Pacific madrone (Arbutus menzeisii)
# 14 Tanoak (Notholithocarpus densiflorus)
# 15 Douglas-fir (Pseudotsuga menzeisii)
# 16 Coast live oak (Quercus agrifolia)
# 17 Blue oak (Quercus douglasii)
# 18 Oregon white oak (Quercus garryana)
# 19 Black oak (Quercus kelloggii)
# 20 Valley oak (Quercus lobata)
# 21 Coast redwood (Sequoia sempvirens)
# 22 California bay laurel (Umbellularia californica)

snames <- c('Builtup','Forest','Grassland','Orchard','Shrubland','Vineyard','Water','WetHerbaceous',NA,NA,'Maple','Buckeye','Madrone','Tanoak','Douglas-fir','Coast live oak','Blue oak','Oregon oak','Black oak','Valley oak','Redwood','California bay')

selspp <- c(3,5,11:22)
rCodes$ShortName[selspp]
dim(pwdv)
for (i in selspp) pwdv <- cbind(pwdv,NA)
dim(pwdv)
nCols <- dim(pwdv)[2]
newCols <- (nCols-length(selspp)+1):nCols
head(pwdv)
names(pwdv)[newCols] <- paste(names(pwdv)[13:26],'_X',sep='')
names(pwdv)

i=1
for (i in 1:length(selspp)) {
  pwdv[which(pwdv$hpspp==selspp[i]),newCols[i]] <- 1
  pwdv[which(pwdv$hpspp!=selspp[i]),newCols[i]] <- 0
  pwdv[which(is.na(pwdv$hpspp)),newCols[i]] <- NA
}

#hptar is vector of species ID for target species - grassland, shrubland, and the 12 tree species
pwdv$hptar <- pwdv$hpspp
pwdv$hptar[pwdv$hpspp %in% c(1,2,4,6,7,8)] <- NA
head(pwdv)
table(pwdv$hpspp)
table(pwdv$hptar)

#add quadratic term for cwd
pwdv$cwd8110sq <- pwdv$cwd8110^2

# dummy variables for geology, only for three major groups
pwdv$block_X <- NA
pwdv$block_X[which(pwdv$rock.group=='block')] <- 1
pwdv$block_X[which(pwdv$rock.group=='uncon')] <- 0
pwdv$block_X[which(pwdv$rock.group=='melange')] <- 0

pwdv$uncon_X <- NA
pwdv$uncon_X[which(pwdv$rock.group=='block')] <- 0
pwdv$uncon_X[which(pwdv$rock.group=='uncon')] <- 1
pwdv$uncon_X[which(pwdv$rock.group=='melange')] <- 0

pwdv$hpX <- hpXY[,1]
pwdv$hpY <- hpXY[,2]

names(pwdv)

write.csv(pwdv,'/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/data/pwd_hyp_topo.csv',row.names=F)
saveRDS(pwdv,'/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/data/pwd_hyp_topo.Rdata')

saveRDS(pwdv,'/Users/david/Google\ Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/big_data/pwd_hyp_topo.Rdata')

