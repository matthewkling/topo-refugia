#########################################################################
# Get species CWD ranges on Pepperwood and plot distribution curves
#########################################################################

library(raster)

# Feb 2017 AVIRIS-NG model GeoTIFFs
#hs.file <- "Z:/Box/Desktop/Projects/Veg Plots/gis/hyperspec/rf_avirisng_merged_spectral_metrics_class.tif"
hs.file <- "Z:/Box/Desktop/Projects/Veg Plots/gis/hyperspec/svm_avirisng_merged_spectral_metrics_class.tif"

# Load the hyperspectral species analysis
hs_species <- raster(hs.file)

# 10m bcm CWD 1981-2010
cwd <- raster("Z:/Box/Desktop/Projects/Climate Models/PWD_10m/cwd1981_2010_ave.asc");cwd[cwd<0] <- NA
crs(cwd) <- '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'
cwd <- projectRaster(cwd,crs=crs(hs_species))

# Extract CWD values for all the hyperspectral map pixels
hs_points <- as.data.frame(rasterToPoints(hs_species))
names(hs_points)[3] <- "class"
coordinates(hs_points) <- 1:2
hs_points$cwd <- extract(cwd,hs_points)

# Throw the CWD values into bins
b<-50
dstr <- data.frame(cwd=seq(from=380+(b/2), to=1380-(b/2), by=b),pixels=as.numeric(table(cut(hs_points$cwd,seq(from=380, to=1380, by=b)))))
dstr$grass <- as.numeric(table(cut(hs_points@data[hs_points$class==3,'cwd'],seq(from=380, to=1380, by=b))))
dstr$shrub <- as.numeric(table(cut(hs_points@data[hs_points$class==5,'cwd'],seq(from=380, to=1380, by=b))))

dstr$Acm <- as.numeric(table(cut(hs_points@data[hs_points$class==11,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Aec <- as.numeric(table(cut(hs_points@data[hs_points$class==12,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Arm <- as.numeric(table(cut(hs_points@data[hs_points$class==13,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Nd <- as.numeric(table(cut(hs_points@data[hs_points$class==14,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Pm <- as.numeric(table(cut(hs_points@data[hs_points$class==15,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Qa <- as.numeric(table(cut(hs_points@data[hs_points$class==16,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Qd <- as.numeric(table(cut(hs_points@data[hs_points$class==17,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Qg <- as.numeric(table(cut(hs_points@data[hs_points$class==18,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Qk <- as.numeric(table(cut(hs_points@data[hs_points$class==19,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Ql <- as.numeric(table(cut(hs_points@data[hs_points$class==20,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Ses <- as.numeric(table(cut(hs_points@data[hs_points$class==21,'cwd'],seq(from=380, to=1380, by=b))))
dstr$Umc <- as.numeric(table(cut(hs_points@data[hs_points$class==22,'cwd'],seq(from=380, to=1380, by=b))))

# Distribtions as proportion of pixels in CWD bins
dstr_prop <- as.data.frame(cbind(cwd=dstr[,1],apply(dstr[,3:16], 2, function(x) x/dstr[,2])))

# Calculate weighted niche means based on proportional occupancy
apply(dstr_prop[,2:15], 2, FUN=function(x) weighted.mean(dstr_prop$cwd,x))

# Color palette for plotting with evergreen-deciduous
classes <- c("Umbcal","Acemac","Quekel","Quedou","Queagr","Psemen","Arbmen","Quegar","Arcman","Notden","Quelob","Aescal","Fralat")
pal <- c("darkgreen","yellow","coral","orange","seagreen3","turquoise4","lawngreen","indianred3","steelblue2","darkolivegreen","violetred","red","yellow2")

# plot(dstr$pixels~dstr$cwd, lwd=2, type="l", main="Distribution over CWD (Pepperwood)", xlab="CWD", ylab="Number of pixels")
# legend("topright", c("all pixels",classes[c(1,3:8,10)]), col=c("black",pal[c(1,3:9)]), lty=1, lwd=2)
# lines(dstr$Pm~dstr$cwd, col=pal[6], lwd=2)

plot(dstr$Pm~dstr$cwd, col=pal[6], lwd=2, type="l", main="Species Distribution over CWD (Pepperwood)", xlab="CWD", ylab="Number of pixels")
lines(dstr$Qd~dstr$cwd, col=pal[4], lwd=2)
lines(dstr$Qk~dstr$cwd, col=pal[3], lwd=2)
lines(dstr$Qa~dstr$cwd, col=pal[5], lwd=2)
lines(dstr$Qg~dstr$cwd, col=pal[8], lwd=2)
lines(dstr$Arm~dstr$cwd, col=pal[7], lwd=2)
lines(dstr$Umc~dstr$cwd, col=pal[1], lwd=2)
lines(dstr$Nd~dstr$cwd, col=pal[9], lwd=2)
# lines(dstr$Ses~dstr$cwd, col=pal[8], lwd=2)
# lines(dstr$Acm~dstr$cwd, col=pal[2], lwd=2)
# lines(dstr$Aec~dstr$cwd, col=pal[12], lwd=2)
# lines(dstr$Ql~dstr$cwd, col=pal[11], lwd=2)
legend("topright", classes[c(1,3:8,10)], col=pal[c(1,3:9)], lty=1, lwd=2)


plot(lowess(dstr_prop$Pm~dstr_prop$cwd,f=0.1), col=pal[6], lwd=2, type="l", main="Species as Percent of CWD Bands (Pepperwood)", xlab="CWD", ylab="Proportion of pixels (lowess fit)", ylim=c(0,0.6))
lines(lowess(dstr_prop$Qd~dstr_prop$cwd,f=0.1), col=pal[4], lwd=2)
lines(lowess(dstr_prop$Qk~dstr_prop$cwd,f=0.1), col=pal[3], lwd=2)
lines(lowess(dstr_prop$Qa~dstr_prop$cwd,f=0.1), col=pal[5], lwd=2)
lines(lowess(dstr_prop$Qg~dstr_prop$cwd,f=0.1), col=pal[8], lwd=2)
lines(lowess(dstr_prop$Arm~dstr_prop$cwd,f=0.1), col=pal[7], lwd=2)
lines(lowess(dstr_prop$Umc~dstr_prop$cwd,f=0.1), col=pal[1], lwd=2)
lines(lowess(dstr_prop$Nd~dstr_prop$cwd,f=0.1), col=pal[9], lwd=2)
# lines(lowess(dstr_prop$Ses~dstr_prop$cwd,f=0.1), col=pal[8], lwd=2)
# lines(lowess(dstr_prop$Acm~dstr_prop$cwd,f=0.1), col=pal[2], lwd=2)
# lines(lowess(dstr_prop$Aec~dstr_prop$cwd,f=0.1), col=pal[12], lwd=2)
# lines(lowess(dstr_prop$Ql~dstr_prop$cwd,f=0.1), col=pal[11], lwd=2)
legend("topright", classes[c(1,3:8,10)], col=pal[c(1,3:9)], lty=1, lwd=2)

