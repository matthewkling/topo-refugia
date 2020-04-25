#hyperspectral to aea and 10m
rm(list=ls())
library(raster)
library(sp)
library(maptools)

PWD <- readShapeSpatial('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/PPshapefiles/PPshapefile-utm10/Pepperwood')
plot(PWD,axes=T)

hp <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/utm10-original/svm_avirisng_merged_spectral_metrics_class.tif')
#plot(hp)
hp
projection(hp) <- CRS('+proj=utm +zone=10N +datum=WGS84')
hpx <- extent(hp)
hpx

rCodes <- read.csv('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/utm10-original/svm_raster_codes.csv',as.is=T)
head(rCodes)
rCodes[rCodes$pall!='grey',c('ShortName','pall')]

source('/Users/david/Documents/Projects/Toolbox/spatial_tools/extent2poly.R')
hpxSP <- extent2poly(hpx)
projection(hpxSP) <- CRS(projection(hp))
plot(hp)
plot(PWD,add=T)
plot(hpxSP,add=T)

#bring in cwd, pwd, and geology in teale-albers as projection to reproject to aea
cwd <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/bcm_outputs/aea/cwd1981_2010_ave.asc')
PWDa <- readShapeSpatial('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/PPshapefiles/PPshapefile-teale-albers/Pepperwood')
cwd[cwd==-9999] <- NA
plot(cwd)
plot(PWDa,add=T)
projection(PWDa) <- CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000')
projection(cwd) <- CRS(projection(PWDa))
mcg <- readRDS('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/PWD_GIS/geology/McLaughlinMap_ta/PepperwoodGeologyQuadShapefile.Rdata')

# reproject hpx extent box
hpxSP.aea <- spTransform(hpxSP,CRS(projection(cwd)))
#plot(cwd)
plot(hpxSP.aea,add=T)


if (FALSE) { # can skip this and read in reprojected hp raster below
  # build empty raster with 2 m resolution and same origin and extent as cwd.hpx
  aea.2m <- raster(ext=extent(cwd.hpx),crs=CRS(projection(cwd)),resolution=2)
  origin(aea.2m) <- origin(cwd.hpx)
  aea.2m
  
  aea.2m <- setValues(aea.2m,1)
  plot(cwd.hpx)
  plot(aea.2m,add=T)
  
  # now have 2m raster with same origin and projection as cwd
  # now reproject hyperspectral map to this new raster with nearest neighbor to preserve canopy ID
  hp.aea <- projectRaster(hp,aea.2m,method='ngb')
  hp.aea
  origin(hp.aea)
  origin(cwd.hpx)
  
  plot(cwd.hpx)
  plot(hp.aea,add=T)
  writeRaster(hp.aea,'/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_spectral_metrics_class_teale-albers_2m.tif',prj=T)
}

hp.aea <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_spectral_metrics_class_teale-albers_2m.tif')
hp.aea

# now crop cwd down to that extent for efficiency
cwd.hpx <- crop(cwd,hpxSP.aea)
plot(cwd.hpx,col=terrain.colors(100))
plot(hpxSP.aea,add=T)
plot(mcg,add=T)
plot(hp.aea,add=T)
plot(PWDa,add=T,border='red',lwd=2)

if (FALSE) {# can skip this and read in aggregate 10 m hp layer below
  #now resample hp.aea to 10m, with most common type assigned to each pixel
  hp.aea10 <- aggregate(hp.aea,5,modal)
  hp.aea10
  cwd.hpx
  origin(hp.aea10)
  origin(cwd.hpx)
  
  plot(cwd.hpx,col=terrain.colors(10))
  plot(hp.aea10,add=T)
  
  writeRaster(hp.aea10,'/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_spectral_metrics_class_teale-albers_10m.tif',prj=T)
}

hp.aea10 <- raster('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_spectral_metrics_class_teale-albers_10m.tif')

#examine close up at boundary
closeup <- extent(cwd.hpx)
closeup@xmin <- -233687.6
closeup@xmax <- -233376.9
closeup@ymin <- 69416.27
closeup@ymax <- 69694.28
plot(crop(cwd.hpx,closeup),col=terrain.colors(10))
plot(crop(hp.aea10,closeup),add=T)

## Perfectly aligned!!!

hyp_poly <- readRDS('/Users/david/Google\ Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/data/hyperspectral-poly.Rdata')

hp.aea10
rCodes <- rCodes[order(rCodes$RasterValue),]

pall <- c('gray','gray','gray','gray','brown','gray','gray','gray','gray','gray','green','blue','red','green','blue','red','green','blue','red','green','blue','red')
length(pall)

#### Figure 2B #######
getwd()
png('/Users/david/Google\ Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/figures/F2b_hyperspectral-map-white-border.png',1500,1570)
op=par(mar=c(1,1,1,5))
plot(hp.aea10,axes=F,col=rCodes$pall)
plot(hyp_poly,add=T,border='red',lwd=3)
#plot(hpxSP.aea,add=T)
plot(PWDa,add=T,lwd=10,border='white')
par(op)
dev.off()
### END FIGURE ####


png('/Users/david/Google\ Drive/Drive-Projects/Topoclimate_Velocity_Project/topo-refugia.git/figures/F2d_PWD-cwd-map.png',1500,1570)
op=par(mar=c(1,1,1,5))
pall <- c("#10002D", "#00008B", "#005954", "#228B22", "#9ACD32", "#CDCD00", "#EEAD0E", "#D2691E", "#8B0000")
cwd.hpx[1] <- 1
cwd.hpx[2] <- 1393
plot(cut(cwd.hpx,breaks=seq(0,1394,length.out=10)),axes=F,col=pall,lwd=3)
plot(hyp_poly,add=T,border='red',lwd=3)
#plot(hp.aea10,add=T)
#plot(hpxSP.aea,add=T)
plot(PWDa,add=T,border='white',lwd=10)
par(op)
dev.off()

## Now create separate 10 m raster for each species, with frequency of cells from 2m as value
# start with 2 m reprojected
hp.aea

calcProportion <- function(x,i) length(which(x==i))/length(x)

table(getValues(hp.aea))
table(getValues(hp.aea10))
if (FALSE) { # skip this and read in result below
  for (i in 1:22) {
    calcProportion <- function(x,v=i,na.rm=T) {
      xx <- x[!is.na(x)]
      length(which(xx==v))/length(xx)
    }
    hp.aea1 <- aggregate(hp.aea,5,calcProportion)
    print(i)
    plot(hp.aea1)
    if (i==1) hp.ind <- stack(hp.aea1) else hp.ind <- addLayer(hp.ind,hp.aea1)
  }
  names(hp.ind) <- rCodes$ShortName
  plot(hp.ind)
  
  writeRaster(hp.ind,'/Users/david/Documents/Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_10m.tif')
}
hp.ind <- stack(raster('/Users/david/Documents/Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_10m.tif',band=1))
for (i in 2:22) hp.ind <- addLayer(hp.ind,raster('/Users/david/Documents/Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_10m.tif',band=i))

# Success - they have identical origin and extent!

#check that new raster has same origin as cwd
origin(hp.aea10)
origin(cwd)
extent(hp.aea10)
extent(cwd)

cwdh <- crop(cwd,hp.aea10)
plot(cwdh)
plot(hp.aea10,add=T)

plot(cwd,col=terrain.colors(20))
plot(hp.aea10,add=T)
hp.aea10

# aggregate to 50m
cwd50 <- aggregate(cwdh,5,mean)
plot(cwd50)

i=1
for (i in 1:22) {
  hp.aea50 <- aggregate(hp.ind[[i]],5,mean)
  if (i==1) hp.ind50 <- stack(hp.aea50) else hp.ind50 <- addLayer(hp.ind50,hp.aea50)
}
names(hp.ind50) <- rCodes$ShortName
plot(hp.ind50)

writeRaster(hp.ind50,'/Users/david/Documents/Projects/Pepperwood/HyperspectralTreeMap/aea/svm_avirisng_merged_ind_class_proportions_teale-albers_50m.tif')
