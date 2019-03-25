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

# is lowest of these correlations significant? YES
op=par(mfrow=c(2,1))
plot(cwd.hypv.opt~reg.cwd.mean,data=d,type='n',
     xlab='Range-wide CWD mean (mm)',
     ylab='Pepperwood CWD optimum (mm)')
text(d$reg.cwd.mean,d$cwd.hypv.opt,labels=d$Plot.abb)
fit <- lm(cwd.hypv.opt~reg.cwd.mean,data=d)
abline(fit)
#abline(0,1,lty=2)
summary(fit)
cor(d$reg.cwd.mean,d$cwd.hypv.opt)


plot(south.mean~reg.cwd.mean,data=d,type='n',
     xlab='Range-wide CWD mean (mm)',
     ylab='Pepperwood southness')
text(d$reg.cwd.mean,d$south.mean,labels=d$Plot.abb)
fit <- lm(south.mean~reg.cwd.mean,data=d)
abline(fit)
#abline(0,1,lty=2)
summary(fit)
cor(d$reg.cwd.mean,d$south.mean)
par(op)

d$Common.name

write.csv(cbind(d[,c('Sci.name','Common.name','Plot.abb')],100*round(d$tot.abund/sum(d$tot.abund),3)),'data/table1.csv',quote=F)

pairs(d[,c('tmin.opt','tmin.gam.mean','model3.mean','topoid.mean','reg.cwd.mean')])
cor(d[,c('tmin.opt','tmin.gam.mean','model3.mean','topoid.mean','reg.cwd.mean')])

plot(d$topoid.mean, d$south.mean)




plot(south.mean~topoid.mean,data=d,type='n',
     xlab='Pepperood topoid mean (mm)',
     ylab='Pepperwood southness')
text(d$topoid.mean, d$south.mean,labels=d$Plot.abb)




plot(reg.ppt_mean~topoid.mean,data=d,type='n',
     xlab='Pepperood topoid mean (mm)',
     ylab='Regional ppt mean (log10 mm)')
text(d$topoid.mean, d$reg.ppt_mean,labels=d$Plot.abb)











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
