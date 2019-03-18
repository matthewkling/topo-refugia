## Frontiers paper step #2 - mean topo position of species at pepperwood

rm(list=ls())
library(raster)

# load 10 m hyperspectral tree layer
hyp <- readRDS('data/pwd_hyp_topo.Rdata')
names(hyp)
dim(hyp)

# look at some of the data
# hist(hyp$cwd8110)
# hist(hyp$Doug.fir_X)

# look at some methods to come up with mean topo position
# weighted.mean(hyp$cwd8110,hyp$Doug.fir_X,na.rm=T)
# hist(hyp$cwd8110)
# hist(hyp$cwd8110[hyp$Doug.fir_X>0.9])
# hist(hyp$cwd8110[hyp$Shrubland_X >0.9])
# hist(hyp$Doug.fir_X)
# table(hyp$Maple_X)
# plot(hyp$Maple,hyp$Maple_X)
#tapply(hyp$Maple,hyp$Maple_X,max,na.rm=T)
#tapply(hyp$Grassland,hyp$Grassland_X,max,na.rm=T)

# Look at distribution of values - mostly fall on the 1/25 sequence, but some in between. These would be due to edges and non-native types excluded from denominator in the aggregation step of the hyperspectral
plot(sort(hyp$Grassland[hyp$Grassland!=0]))

## THERE"S SOMETHING WRONG WITH THE '_X' variables - ignoring them for now
# They are supposed to be the binary species assignment, based on most common 

# Checking the sum of relative frequencies of species of interest in each cell
# Because some types are excluded - water, ag, etc. - not all sum to one. This has to be dealt with below in some of the weighting schemes
names(hyp)
hyp$vegt <- apply(hyp[,14:27],1,sum)
hyp$woodyt <- apply(hyp[,15:27],1,sum)
#hyp$vegXt <- apply(hyp[,31:43],1,sum)
hist(hyp$vegt)
hist(hyp$woodyt)
length(which(hyp$woodyt==1))
length(which(hyp$woodyt>=0.5))
#hist(hyp$vegXt)

# subsequent analyses will be done one cells with at least 50% woody native veg types (including Shrubland) and have values for southness and cwd
xx <- complete.cases(hyp[,c('cwd8110','southness')])
length(xx)
length(which(!is.na(hyp$cwd8110)))
length(which(!is.na(hyp$southness)))
length(intersect(which(!is.na(hyp$southness)), which(!is.na(hyp$cwd8110))))
# something weird here - complete.cases is longer than either column. But this shows that southness has fewer cases and is a subset of cwd, so we'll use that to subset data.frame

# subset data frame
hypv <- hyp[intersect(which(hyp$woodyt>=0.5) , which(!is.na(hyp$southness))),]
dim(hypv)
names(hypv)

# set up breaks for histograms, either equally spaced or by even percentiles of the abundances; look at relationship between cwd and southness

summary(hypv$cwd8110)
summary(hypv$southness)
hist(hypv$cwd8110)
hist(hypv$southness,breaks=seq(-0.75,0.7,by=0.05))
rsamp <- sample(nrow(hypv),10000)
plot(hypv[rsamp,c('southness','cwd8110')])
plot(hypv[rsamp,c('topoid','model3')])

# set up bins for niche means
# this is evenly spaced bins from environment factor histogram
summary(hypv$cwd8110)
cwd.even <- hist(hypv$cwd8110,breaks=c(seq(325,1425,by=50)))
cwd.even$breaks
cwd.even$counts
cwd.even
hypv$cwd.even.fac <- cut(hypv$cwd8110,cwd.even$breaks)

summary(hypv$southness)
south.even <- hist(hypv$southness,breaks=c(seq(-0.75,0.7,by=0.05)))
south.even$breaks
south.even$counts
hypv$south.even.fac <- cut(hypv$southness,south.even$breaks)
length(south.even.fac)

summary(hypv$topoid)
topoid.even <- hist(hypv$topoid,breaks=c(seq(0,25,by=0.5)))
topoid.even$breaks
topoid.even$counts


# sum veg occupancy by bins
vegC <- data.frame(matrix(NA,nrow = length(cwd.even$counts),ncol=13))
names(vegC) <- names(hypv[15:27])

vegS <- data.frame(matrix(NA,nrow = length(south.even$counts),ncol=13))
names(vegS) <- names(hypv[15:27])

names(hypv)
i=13
for (i in 1:13) {
  message(i)
  vegC[,i] <- tapply(hypv[,14+i],hypv$cwd.even.fac,sum)
  vegS[,i] <- tapply(hypv[,14+i],hypv$south.even.fac,sum)
}
vegC.tot <- apply(vegC,1,sum)
vegS.tot <- apply(vegS,1,sum)

# now order all pixels by sum of woody plant occupancies, and find breaks that divide them into 25 evenly weighted perentiles

hyp.Csort <- hypv[order(hypv$cwd8110),]
hyp.Csort$wcumsum <- cumsum(hypv$woodyt)
plot(hyp.Csort$cwd8110)

head(hyp.Csort$wcumsum)
tail(hyp.Csort$wcumsum)
wcs.tot <- sum(hyp.Csort$woodyt)

cwd.breaks <- rep(NA,26)
for (i in 1:26) cwd.breaks[i] <- which(hyp.Csort$wcumsum>=(i-1)/25*wcs.tot)[1]
cwd.breaks
hyp.Csort$cwd.breaks.fac <- cut(hyp.Csort$cwd8110,breaks=hyp.Csort$cwd8110[cwd.breaks])
table(hyp.Csort$cwd.breaks.fac)

cwd.occ.perc.mids <- rep(NA,25)
for (i in 1:25) cwd.occ.perc.mids[i] <- mean(hyp.Csort$cwd8110[cwd.breaks[i]:cwd.breaks[i+1]])
cwd.occ.perc.mids

vegC.occ.perc <- data.frame(matrix(NA,nrow = length(cwd.occ.perc.mids),ncol=13))
names(vegC.occ.perc) <- names(hypv[15:27])

i=1
for (i in 1:13) {
  vegC.occ.perc[,i] <- tapply(hyp.Csort[,14+i],hyp.Csort$cwd.breaks.fac,sum)
}
vegC.occ.perc
vegC.occ.perc$tot <- apply(vegC.occ.perc,1,sum)

## calculate absolute and fractionally weighted cwd means by veg
# decision is to use fractional weights as much as possible
# three methods to calculate niche mean. Can do each one for cwd and for southness
#     1 Weighted mean of cwd, weighted by veg fraction
#     2 mean, weighted by percent occupancy within bins, with even bins
#     3 mean, weighted by % occupancy, with percentile bins so each one has same percentage of total distribution - not sure yet how to do this with fractional percentages

niche <- data.frame(matrix(NA,13,8))
names(niche) <- c('cwd.mean','cwd.occ.even','cwd.occ.perc','south.mean','south.occ.even','south.occ.perc','cwd.count','cwd.frac')
row.names(niche) <- names(hypv[,15:27])

i=1
for (i in 1:13) {
  niche$cwd.mean[i] <- weighted.mean(hypv$cwd8110,hypv[,14+i])
  niche$cwd.occ.even[i] <- weighted.mean(cwd.even$mids,vegC[,i]/vegC.tot)
  niche$cwd.occ.perc[i] <- weighted.mean(cwd.occ.perc.mids,vegC.occ.perc[,i]/vegC.occ.perc$tot)
  
  niche$south.mean[i] <- weighted.mean(hypv$southness,hypv[,14+i])
  niche$south.occ.even[i] <- weighted.mean(south.even$mids,vegS[,i]/vegS.tot)
  
  {
    # from previous test calculations
  #niche$cwdX[i] <- weighted.mean(hypv$cwd8110,hypv[,28+i])
    niche$cwd.count[i] <- weighted.mean(cwd.breaks$mids,vegC[,i])
   #niche$cwdX.frac[i] <- weighted.mean(cwdb$mids,vegCX[,i]/cwdb$counts)
    niche$cwd.frac[i] <- weighted.mean(cwd.breaks$mids,vegC[,i]/cwd.breaks$counts)
  }
}
niche
pairs(niche[,c(1,2,4,5)])

pdf('figures/pwd_cwdfrachist.pdf',8,20)
op=par(mfrow=c(14,1),mar=c(3,3,0,0))
for (i in 1:13) barplot(vegC[,i]/cwd.even$counts)
par(op)
dev.off()

pdf('figures/pwd_cwdhist.pdf',8,20)
op=par(mfrow=c(13,1),mar=c(3,3,0,0))
for (i in 1:14) barplot(vegC[,i])
par(op)
dev.off()

pdf('figures/pwd_southfrachist.pdf',8,20)
op=par(mfrow=c(13,1),mar=c(3,3,0,0))
for (i in 1:14) barplot(vegS[,i]/south.even$counts)
par(op)
dev.off()

pdf('figures/pwd_southhist.pdf',8,20)
op=par(mfrow=c(13,1),mar=c(3,3,0,0))
for (i in 1:14) barplot(vegS[,i])
par(op)
dev.off()

pdf('figures/pwd_topo_nichepairs.pdf',8,8)
#op=par(mfrow=c(14,1),mar=c(3,3,0,0))
pairs(niche[,c(1,2,4,5)])
#par(op)
dev.off()

rownames(niche)
sci.names <- c('Adenostoma fasciculatum','Acer macrophyllum','Aesculus californica','Arbutus menziesii','Notholithocarpus densiflorus','Pseudotsuga menziesii','Quercus agrifolia','Quercus douglasii','Quercus garryana','Quercus kelloggii','Quercus lobata','Sequoia sempervirens','Umbellularia californica')

write.csv(cbind(sci.names,round(niche[,c(1,2)],0),round(niche[,c(4,5)],3)),'data/pwd_niche_means.csv')

round(niche[,c(1,2,5,4)],0)


## does niche mean shift across rock types?
{
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
}