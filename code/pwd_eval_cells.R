# evaluate climate change magnitude under different scenarios

clim <- read.csv('data/pwd_climate_1km.csv',as.is=T)
head(clim)
clims <- data.frame(scenario=unique(clim$scenario),cwd=NA,aet=NA,tminmin=NA,ppt=NA,djf=NA,jja=NA)
for (c in c('cwd','aet','tminmin','ppt','djf','jja')){
  clims[,c] <- tapply(clim[,c],clim$scenario,mean)
}
clims
pairs(clims[,c('cwd','aet','djf')])
cor(clims[,c('cwd','aet','djf')])

# summary of scenarios in terms of range of changes in climate factors
# cwd: historic = 716, max fut = 831
#     for future exploration: 0, +40, +80, +120
# aet: historic 415, min/max fut = 391/450
#     for future exploration:  -25, 0, 25, 50
# tminmin: historic 4.4, max fut 7.6
#     for future exploration: 0, +1, +2, +3
#     same works for djf


head(clim)
hclim <- clim[clim$scenario=='historic',]
dclim <- clim[clim$scenario!='historic',]
c='cwd'
for (c in c('cwd','aet','tminmin','ppt','djf')){
  xx <- hclim[match(dclim$cell,hclim$cell),c]
  dclim[,c] <- dclim[,c]-xx
}
dclim$scencell <- paste0(dclim$scenario,dclim$cell)

deltas <- read.csv("data/pwd_niche_deltas.csv",as.is=T)
deltas$scencell <- paste0(deltas$scenario,deltas$cell)
head(deltas)
head(dclim)

d2d <- match(deltas$scencell,dclim$scencell)
deltas$dcwd <- dclim$cwd[d2d]
deltas$dtminmin <- dclim$tminmin[d2d]
deltas$dppt <- dclim$ppt[d2d]
head(deltas)

unique(deltas$species)
sp <- unique(deltas$species)[12]
plot(delta~dcwd,data=deltas[deltas$species==sp 
                            & deltas$algorithm=='maxent'
                            & deltas$var_set=='cwd_ppt_djf_jja',])
plot(delta~dtminmin,data=deltas[deltas$species==sp 
                                & deltas$algorithm=='maxent'
                                & deltas$var_set=='cwd_ppt_djf_jja',])
plot(delta~dppt,data=deltas[deltas$species==sp 
                                & deltas$algorithm=='maxent'
                                & deltas$var_set=='cwd_ppt_djf_jja',])
