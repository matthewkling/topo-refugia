## compare PWD and regional means
rm(list=ls())

names <- read.csv('data/names.csv',as.is=T)
names

pwdgam <- read.csv('data/pwd_gam1.csv',as.is=T)
pwdgam

pwdtopo <- read.csv('data/pwd_niche_means.csv',as.is=T)
pwdtopo

pwd <- merge(pwdgam,pwdtopo)
pwd
all(pwd$species == pwd$Hyp.name)

names(pwd)
plot(pwd$tot.abund,pwd$gam.tabund)
plot(pwd$hypv.pmax,pwd$cspace.pmax)
abline(0,1)


reg <- read.csv('data/regional_niche_stats.csv',as.is=T)
reg
regx <- data.frame(sci.names=reg$species[which(reg$stat=='mean')],
                   reg.cwd.mean=reg$clim[which(reg$stat=='mean')],
                   reg.cwd.median=reg$clim[which(reg$stat=='median')],
                   reg.cwd.max=reg$clim[which(reg$stat=='max')])
regx

pwd.reg <- merge(pwd,regx)
names(pwd.reg)

## Variable name code
# gam.tabund = sum of predicted values across PWD environmental space
# cwd.opt = PWD CWD value at peak of GAM model, from orthogonal espace
# tmin.opt = PWD Tmin (model3) value at peak of GAM model, from orthogonal climate space
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
pairs(pwd.reg[,c('cwd.opt','cwd.gam.mean','cwd.mean','south.mean','reg.cwd.mean','reg.cwd.median','reg.cwd.max')])
cor(pwd.reg[,c('cwd.opt','cwd.gam.mean','cwd.mean','south.mean','reg.cwd.mean','reg.cwd.median','reg.cwd.max')])

pairs(pwd.reg[,c('tmin.opt','tmin.gam.mean','model3.mean','topoid.mean','reg.cwd.mean')])
cor(pwd.reg[,c('tmin.opt','tmin.gam.mean','model3.mean','topoid.mean','reg.cwd.mean')])














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
