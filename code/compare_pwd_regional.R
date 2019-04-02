## compare PWD and regional means
rm(list=ls())
library(tidyverse)
select <- dplyr::select

names <- read.csv('data/names.csv',as.is=T) %>% 
      select(Sci.name, Common.name, Plot.abb)
pwdgam <- read.csv('data/pwd_gam1.csv',as.is=T) %>% 
      rename(Sci.name = sci.names) %>% 
      select(-X)
pwdtopo <- read.csv('data/pwd_niche_means.csv',as.is=T) %>% 
      rename(Sci.name = sci.names) %>% 
      select(-X, -Hyp.name)
reg <- read.csv('data/regional_niche_stats.csv',as.is=T) %>% 
      rename(Sci.name = species) %>% 
      gather(summary, value, mean:max) %>%
      unite(stat, var, summary) %>%
      mutate(stat = paste0("reg.", stat)) %>%
      spread(stat, value)

d <- names %>%
      left_join(pwdgam) %>%
      left_join(pwdtopo) %>%
      left_join(reg)

names(d)
plot(d$tot.abund,d$gam.tabund)
plot(d$hypv2.pmax,d$cspace2.pmax)
abline(0,1)

## Variable name code
# gam.tabund = sum of predicted values across PWD environmental space
# cwd.cspace.opt = PWD CWD value at peak of GAM model, from orthogonal espace
# tmin.cspace.opt = PWD Tmin (model3) value at peak of GAM model, from orthogonal climate space
# cwd.hypv.opt = PWD CWD value at peak of GAM model, from observed espace
# tmin.hypv.opt = PWD Tmin (model3) value at peak of GAM model, from observed climate space
# cwd.gam.mean = weighted mean CWD value based on GAM predicted values, across PWD espace
# tmin.gam.mean = weighted mean Tmin value based on GAM predicted values, across PWD espace
# hypv.pmax = maximum predicted value across PWD espace
# cspace.pmax = maximum predicted value across orthogonal espace
# dev.expl = deviance explained by gam
# cwd.chisq and tmin.chisq = GAM chi squared values for the two variables
# tot.abund = sum of observed abundances across PWD espace
# cwd.mean = mean of CWD, weighted by observed abundances
# south.mean = mean of southness, weighted by observed abundances
# topoid.mean = mean of topographic wetness index, weighted by observed abundances
# model3.mean = mean of Tmin (model3), weighted by observed abundances
# reg.cwd.mean = reginal CWD niche model, mean
# reg.cwd.median = reginal CWD niche model, median
# reg.cwd.max = reginal CWD niche model, max

# Have a look at scatterplots and correlations
pairs(d[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd1.gam.mean','south.mean','reg.cwd_mean','reg.cwd_median','reg.cwd_max')])
cor(d[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd1.gam.mean','south.mean','reg.cwd_mean','reg.cwd_median','reg.cwd_max')])

# compare regional niche metrics
pairs(d[,c('reg.cwd_mean','reg.cwd_median','reg.cwd_max')])
cor(d[,c('reg.cwd_mean','reg.cwd_median','reg.cwd_max')])

# compare regional niche metrics
pairs(d[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd1.gam.mean','south.mean')])
cor(d[,c('cwd1.hypv.opt','cwd2.hypv.opt','cwd1.gam.mean','south.mean')])
summary(d$cwd1.gam.mean)
sort(d$cwd1.gam.mean)
summary(d$cwd1.hypv.opt)
sort(d$cwd1.hypv.opt)

# is lowest of these correlations significant? YES
op=par(mfrow=c(2,1))
plot(cwd1.hypv.opt~reg.cwd_mean,data=d,type='n',
     xlab='Range-wide CWD mean (mm)',
     ylab='Pepperwood CWD optimum (mm)')
text(d$reg.cwd_mean,d$cwd1.hypv.opt,labels=d$Plot.abb)
fit <- lm(cwd1.hypv.opt~reg.cwd_mean,data=d)
abline(fit)
#abline(0,1,lty=2)
summary(fit)
cor(d$reg.cwd_mean,d$cwd1.hypv.opt)


plot(south.mean~reg.cwd_mean,data=d,type='n',
     xlab='Range-wide CWD mean (mm)',
     ylab='Pepperwood southness')
text(d$reg.cwd_mean,d$south.mean,labels=d$Plot.abb)
fit <- lm(south.mean~reg.cwd_mean,data=d)
abline(fit)
#abline(0,1,lty=2)
summary(fit)
cor(d$reg.cwd_mean,d$south.mean)
par(op)

d$Common.name

#write.csv(cbind(d[,c('Sci.name','Common.name','Plot.abb')],100*round(d$tot.abund/sum(d$tot.abund),3)),'data/table1.csv',quote=F)

names(d)
pairs(d[,c('tmin2.hypv.opt','tmin2.gam.mean','topoid.mean','reg.cwd_mean')])
cor(d[,c('tmin2.hypv.opt','tmin2.gam.mean','topoid.mean','reg.cwd_mean')])


plot(south.mean~topoid.mean,data=d,type='n',
     xlab='Pepperood topoid mean (mm)',
     ylab='Pepperwood southness')
text(d$topoid.mean, d$south.mean,labels=d$Plot.abb)


plot(reg.ppt_mean~topoid.mean,data=d,type='n',
     xlab='Pepperood topoid mean (mm)',
     ylab='Regional ppt mean (log10 mm)')
text(d$topoid.mean, d$reg.ppt_mean,labels=d$Plot.abb)






# niche model deltas vs. pwd topoclimate

deltas <- read.csv("data/pwd_niche_deltas.csv",as.is=T) %>%
      group_by(species, scenario, var_set, algorithm) %>%
      summarize(delta=mean(delta)) %>%
      rename(Sci.name = species) %>%
      left_join(names) %>%
      left_join(pwdgam) %>%
      left_join(pwdtopo) %>%
      separate(scenario, c("year", "rcp", "model"), sep="_")
dim(deltas)
names(deltas)
unique(deltas$var_set)

ggplot(deltas, aes(cwd1.gam.mean, delta)) +
      geom_hline(yintercept=0, linetype=2) +
      geom_smooth(method=lm) +
      geom_smooth(method=lm, level=0.5, alpha=1) +
      facet_grid(algorithm + var_set ~ model + rcp, scales="free") +
      geom_text(aes(label=Plot.abb), size=2) +
      theme_bw()

algs <- unique(deltas$algorithm)
rcps <- unique(deltas$rcp)
models <- unique(deltas$model)
vars <- unique(deltas$var_set)

a=algs[1]
r=rcps[1]
m=models[1]
v=vars[1]
aa <- c()
rr <- c()
mm <- c()
vv <- c()
slps <- c()
pvals <- c()
for (a in algs)
  for (r in rcps)
  for (m in models)
    for (v in vars){
      aa <- c(aa,a);rr <- c(rr,r); mm <- c(mm,m);vv <- c(vv,v)      
      xx <- deltas[deltas$algorithm==a & deltas$rcp == r & deltas$model==m & deltas$var_set == v,]
      fit <- lm(delta~cwd1.gam.mean,data=xx)
      slps <- c(slps,summary(fit)$coeff[2,1])
      pvals <- c(pvals,summary(fit)$coeff[2,4])
      #lmres <- rbind(lmres,c(a,m,v,slp,pval))
    }
lmres <- data.frame(algorithm=aa,rcp=rr,model=mm,var_set=vv,slp=slps,pval=pvals)
lmres$scenario <- paste0('2061-2080_',lmres$rcp,'_',lmres$model)
head(lmres)
plot(lmres$slp,lmres$pval)

# evaluate the changes occurring at pepperwood under different scenarios
clim <- read.csv('data/pwd_climate_1km.csv')
head(clim)
clims <- data.frame(scenario=unique(clim$scenario),cwd=NA,aet=NA,tminmin=NA,ppt=NA,djf=NA,jja=NA)
for (c in c('cwd','aet','tminmin','ppt','djf','jja')){
  clims[,c] <- tapply(clim[,c],clim$scenario,mean)
}
clims

plot(clims$cwd,clims$aet)
abline(v=rev(clims$cwd)[1])
abline(h=rev(clims$aet)[1])

plot(clims$cwd,clims$ppt)
abline(v=rev(clims$cwd)[1])
abline(h=rev(clims$ppt)[1])

# are deltas a function of degree of climate change
lmresc <- merge(clims,lmres)
xx <- which(lmresc$rcp=='rcp85' & lmresc$algorithm=='maxent')
plot(lmresc$cwd[xx],lmresc$slp[xx])
yy <- intersect(xx,lmresc$pval<=0.05)

lmres[order(lmres$pval),]
#### SCRIPT BELOW FROM BEFORE RENAMING VARS
# plot(pwd.reg$cwd.occ.even~pwd.reg$clim)
# fit <- lm(pwd.reg$cwd.occ.even~pwd.reg$clim)
# summary(fit)
# abline(fit)
# 
# pwdg.reg <- merge(pwdgam,pwd.reg,by.x='sci.name',by.y='sci.names')
# pwdg.reg
# plot(pwdg.reg$cwd_opt,pwdg.reg$cwd.occ.even)
# plot(pwdg.reg$cwd_opt,pwdg.reg$cwd.mean.x)
# plot(pwdg.reg$cwd_opt,pwdg.reg$south.mean)
# plot(pwdg.reg$cwd_opt,pwdg.reg$clim)
# 
# names(pwdg.reg)
# pairs(pwdg.reg[,c('cwd_opt','cwd.mean.x','cwd.mean.y','cwd.occ.even','south.mean','south.occ.even','clim')])
# (cwd_corr_matrix <- cor(pwdg.reg[,c('cwd_opt','cwd.mean.x','cwd.mean.y','cwd.occ.even','south.mean','south.occ.even','clim')]))
# min(cwd_corr_matrix)
# 
# # lowest correlation between regional clim mean and any local cwd or southness measure at pepperwood
# plot(pwdg.reg[,c('clim','south.mean')])
# fit=lm(pwdg.reg$south.mean~pwdg.reg$clim)
# abline(fit)
# summary(fit)
# 
# pairs(pwdg.reg[,c('tmin_opt','tmin.mean')])
# (cor(pwdg.reg[,c('tmin_opt','tmin.mean')]))
# 
# fit <- lm(pwdg.reg$clim~pwdg.reg$cwd_opt)
# abline(fit)
# summary(fit)
# 
# plot(pwdg.reg[,c('cwd.mean.x','cwd.mean.y')])
# fit <- lm(pwdg.reg$cwd.mean.y~pwdg.reg$cwd.mean.x)
# abline(fit)
# summary(fit)
# 
# plot(pwdg.reg[,c('cwd_opt','cwd.mean.y')])
# fit <- lm(pwdg.reg$cwd.mean.y~pwdg.reg$cwd_opt)
# abline(fit)
# summary(fit)
