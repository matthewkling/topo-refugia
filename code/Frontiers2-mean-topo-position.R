## Frontiers #2 - mean topo position of species at pepperwood

rm(list=ls())
library(raster)

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

# load 10 m hyperspectral tree layer
hyp <- readRDS('/Users/david/Google\ Drive/Drive-Projects/Pepperwood/HyperspectralTreeMap/data/pwd_hyp_topo.Rdata')
names(hyp)
dim(hyp)

hist(hyp$cwd8110)
hist(hyp$Doug.fir_X)

weighted.mean(hyp$cwd8110,hyp$Doug.fir_X,na.rm=T)
hist(hyp$cwd8110)
hist(hyp$cwd8110[hyp$Doug.fir_X>0.9])
hist(hyp$cwd8110[hyp$Shrubland_X >0.9])
hist(hyp$Doug.fir_X)
table(hyp$Maple_X)
plot(hyp$Maple,hyp$Maple_X)

head(hyp)
hyp$vegt <- apply(hyp[,13:26],1,sum)
hyp$vegXt <- apply(hyp[,29:42],1,sum)
hist(hyp$vegt)

hypv <- hyp[!is.na(hyp$vegXt)&!is.na(hyp$cwd8110),]

hypv$cwd50 <- cut(hypv$cwd8110,breaks=seq(325,1425,by=50))
summary(hypv$cwd50)

summary(hypv$cwd8110)
cwdb <- hist(hypv$cwd8110,breaks=c(seq(325,1425,by=50)))
cwdb$breaks
cwdb$counts

barplot(cwdb$counts)
cwdb

vegCX <- data.frame(matrix(NA,nrow = length(cwdb$counts),ncol=14))
names(vegC) <- names(hypv[13:26])

vegC <- data.frame(matrix(NA,nrow = length(cwdb$counts),ncol=14))
names(vegC) <- names(hypv[13:26])

i=14
for (i in 1:14) {
  #print(i)
  vegCX[,i] <- hist(hypv$cwd8110[hypv[,28+i]>0.9],breaks=cwdb$breaks)$counts
  vegC[,i] <- tapply(hypv[,12+i],hypv$cwd50,sum)
}

## calculate absolute and fractionally weighted cwd means by veg
niche <- data.frame(matrix(NA,14,5))
names(niche) <- c('cwd','cwdX','cwd.count','cwdX.frac','cwd.frac')
row.names(niche) <- names(hypv[,13:26])

i=1
for (i in 1:14) {
  niche$cwd[i] <- weighted.mean(hypv$cwd8110,hypv[,12+i])
  niche$cwdX[i] <- weighted.mean(hypv$cwd8110,hypv[,28+i])
  niche$cwd.count[i] <- weighted.mean(cwdb$mids,vegC[,i])
  niche$cwdX.frac[i] <- weighted.mean(cwdb$mids,vegCX[,i]/cwdb$counts)
  niche$cwd.frac[i] <- weighted.mean(cwdb$mids,vegC[,i]/cwdb$counts)
}
niche
pairs(niche)

pdf('figsFrontiers2019/vegfrachist.pdf',8,20)
op=par(mfrow=c(14,1),mar=c(3,3,0,0))
for (i in 1:14) barplot(vegC[,i]/cwdb$counts)
par(op)
dev.off()

pdf('figsFrontiers2019/veghist.pdf',8,20)
op=par(mfrow=c(14,1),mar=c(3,3,0,0))
for (i in 1:14) barplot(vegC[,i])
par(op)
dev.off()

pdf('figsFrontiers2019/nichepairs.pdf',8,8)
#op=par(mfrow=c(14,1),mar=c(3,3,0,0))
pairs(niche[,c(1,2,4)])
#par(op)
dev.off()

round(niche[,c(1,2,5,4)],0)

## does niche mean shift across rock types?
table(hypv$rock.group)

niche$cwd.block <- NA
niche$cwd.melange <- NA
niche$cwd.uncon <- NA
niche$cwd.ultra <- NA

block.rows <- which(hypv$rock.group=='block')
melange.rows <- which(hypv$rock.group=='melange')
uncon.rows <- which(hypv$rock.group=='uncon')
ultra.rows <- which(hypv$rock.group=='ultra')

i=1
for (i in 1:14) {
  niche$cwd.block[i] <- weighted.mean(hypv$cwd8110[block.rows],hypv[block.rows,12+i])
  niche$cwd.melange[i] <- weighted.mean(hypv$cwd8110[melange.rows],hypv[melange.rows,12+i])
  niche$cwd.uncon[i] <- weighted.mean(hypv$cwd8110[uncon.rows],hypv[uncon.rows,12+i])
  niche$cwd.ultra[i] <- weighted.mean(hypv$cwd8110[ultra.rows],hypv[ultra.rows,12+i])
}
niche
op=par(mfrow=c(2,2))
plot(niche$cwd.block~niche$cwd);abline(0,1)
plot(niche$cwd.melange~niche$cwd);abline(0,1)
plot(niche$cwd.uncon~niche$cwd);abline(0,1)
plot(niche$cwd.ultra~niche$cwd);abline(0,1)
par(op)

# on the face of it, it appears that all veg types shift to hotter locations on uncon and ultra rock, and cooler locations on block and melange. But is it just due to non-random distribution overll? Check cross-tab of cwd by rock types:
cwd.rock <- tapply(hypv$cwd8110,hypv$rock.group,mean)
cwd.rock - mean(hypv$cwd8110)

# replot using these mean displacements to see if this explains the shifts
op=par(mfrow=c(2,2))
plot(niche$cwd.block~niche$cwd);abline(-32.9,1)
plot(niche$cwd.melange~niche$cwd);abline(-7.9,1)
plot(niche$cwd.uncon~niche$cwd);abline(21.3,1)
plot(niche$cwd.ultra~niche$cwd);abline(68.1,1)
par(op)

# LOOKS lime the overall nobn-random pretty much accounts for shifts in each veg type

# what about the overall distribution of vegtypes across rock types
rock.veg.obs <- data.frame(matrix(NA,14,4))
names(rock.veg.obs) <- names(table(hypv$rock.group))
row.names(rock.veg.obs) <- names(hypv[13:26])
for (i in 1:14) rock.veg.obs[i,] <- tapply(hypv[,12+i],hypv$rock.group,sum)
round(rock.veg.obs,0)

rock.veg.exp <- rock.veg.obs
rock.freq <- apply(rock.veg.obs,2,sum)/sum(rock.veg.obs)
rock.freq
veg.freq <- apply(rock.veg.obs,1,sum)/sum(rock.veg.obs)
veg.freq
for (i in 1:14) for (j in 1:4) rock.veg.exp[i,j] <- sum(rock.veg.obs) * rock.freq[j] * veg.freq[i] 
round(rock.veg.exp)
round(rock.veg.obs,0)
round(rock.veg.obs,0) - round(rock.veg.exp)
sum(round(rock.veg.obs,0) - round(rock.veg.exp))
