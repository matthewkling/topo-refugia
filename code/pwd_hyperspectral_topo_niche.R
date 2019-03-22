## Frontiers paper step #2 - mean topo position of species at pepperwood

rm(list=ls())
library(raster)

# load 10 m hyperspectral tree layer
hyp <- readRDS('data/pwd_hyp_topo.Rdata')
names(hyp)
dim(hyp)
nhyp <- length(which(!is.na(hyp$hpspp)))
# area of hyperspectral image
(nhyp * 4)/10000
length(which(is.na(hyp$hpspp)))
summary(hyp$hpspp)
summary(hyp$elevation[which(!is.na(hyp$hpspp))])
summary(hyp$southness[which(!is.na(hyp$hpspp))])
summary(hyp$topoid[which(!is.na(hyp$hpspp))])

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

# slow plot, uncomment to see
#plot(sort(hyp$Grassland[hyp$Grassland!=0]))

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
env.complete <- intersect(which(!is.na(hyp$southness)) , which(!is.na(hyp$topoid)))
hypv <- hyp[intersect(which(hyp$woodyt>=0.5) , env.complete),]
dim(hypv)
hypv <- hypv[-which(hypv$rock.group=='ultra'),]
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


## calculate niche means by direct averaging of environmental values for all pixels, weighted by specieas occurrence in the pixel. Not using 'X' vars

# names from hyperspectral data set, and then set up scientific name vector
(Hyp.name=names(hypv[,15:27]))
sci.names <- c('Adenostoma fasciculatum','Acer macrophyllum','Aesculus californica','Arbutus menziesii','Notholithocarpus densiflorus','Pseudotsuga menziesii','Quercus agrifolia','Quercus douglasii','Quercus garryana','Quercus kelloggii','Quercus lobata','Sequoia sempervirens','Umbellularia californica')

niche <- data.frame(Hyp.name,sci.names,tot.abund=NA,cwd.mean=NA,south.mean=NA,topoid.mean=NA,model3.mean=NA)
niche

i=1
for (i in 1:13) {
  sp <- Hyp.name[i]
  message(sp)
  niche$tot.abund[i] <- sum(hypv[,sp])
  niche$cwd.mean[i] <- weighted.mean(hypv$cwd8110,hypv[,sp])
  niche$south.mean[i] <- weighted.mean(hypv$southness,hypv[,sp])
  niche$topoid.mean[i] <- weighted.mean(hypv$topoid,hypv[,sp])
  niche$model3.mean[i] <- weighted.mean(hypv$model3,hypv[,sp])
}
niche
pairs(niche[,4:7])
cor(niche[,4:7])


pdf('figures/pwd_topo_nichepairs.pdf',8,8)
pairs(niche[,4:7])
dev.off()
system('open figures/pwd_topo_nichepairs.pdf')
rownames(niche)

write.csv(niche,'data/pwd_niche_means.csv')


## does niche mean shift across rock types?
if (FALSE) {
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