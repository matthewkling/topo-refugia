# load cwd and southness and compare
cwd10 <- raster('/Users/david/Google Drive/Drive-Projects/Pepperwood/PWD_GIS/BCM_10m/bcm_outputs/aea/cwd1981_2010_ave.asc')
projection(cwd10)
plot(cwd10)
cwd10[cwd10<0] <- NA
hist(getValues(cwd10))

dem <- raster('/Users/david/Google Drive/Drive-Projects/Pepperwood/PWD_GIS/dem/BCM10m/pwd_10m_t1.asc')
plot(dem)

cwd10 <- crop(cwd10,dem)
plot(cwd10)

projection(dem) <- CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000')
slp <- terrain(dem,'slope')
asp <- terrain(dem,'aspect')
south <- -1*cos(asp)*sin(slp)
plot(south)
plot(getValues(cwd10)~getValues(south))
cor(getValues(cwd10),getValues(south),use='pair')
length(getValues(cwd10))