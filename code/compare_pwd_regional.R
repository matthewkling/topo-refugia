## compare PWD and regional means
rm(list=ls())

names <- read.csv('data/names.csv',as.is=T)
names

pwdgam <- read.csv('data/pwd_gam1.csv',as.is=T)
pwdgam
pwdgam$sci.name <- names$Sci.name[match(pwdgam$species,names$Hyp.name)]

pwd <- read.csv('data/pwd_niche_means.csv',as.is=T)
pwd

pwdm <- merge(pwdgam,pwd,by.x='sci.name',by.y='sci.names')
pwdm
names(pwdm)
pairs(pwdm[,c('cwd_opt','cwd.mean.x','cwd.mean.y','cwd.occ.even','south.mean','south.occ.even')])

reg <- read.csv('data/regional_niche_stats.csv',as.is=T)
reg

pwd.reg <- merge(pwd,reg[which(reg$stat=='mean'),],by.x='sci.names',by.y='species')
pwd.reg
plot(pwd.reg$cwd.occ.even~pwd.reg$clim)
fit <- lm(pwd.reg$cwd.occ.even~pwd.reg$clim)
summary(fit)
abline(fit)

pwdg.reg <- merge(pwdgam,pwd.reg,by.x='sci.name',by.y='sci.names')
pwdg.reg
plot(pwdg.reg$cwd_opt,pwdg.reg$cwd.occ.even)
plot(pwdg.reg$cwd_opt,pwdg.reg$cwd.mean.x)
plot(pwdg.reg$cwd_opt,pwdg.reg$south.mean)
plot(pwdg.reg$cwd_opt,pwdg.reg$clim)

names(pwdg.reg)
pairs(pwdg.reg[,c('cwd_opt','cwd.mean.x','cwd.mean.y','cwd.occ.even','south.mean','south.occ.even','clim')])
(cwd_corr_matrix <- cor(pwdg.reg[,c('cwd_opt','cwd.mean.x','cwd.mean.y','cwd.occ.even','south.mean','south.occ.even','clim')]))
min(cwd_corr_matrix)

# lowest correlation between regional clim mean and any local cwd or southness measure at pepperwood
plot(pwdg.reg[,c('clim','south.mean')])
fit=lm(pwdg.reg$south.mean~pwdg.reg$clim)
abline(fit)
summary(fit)

pairs(pwdg.reg[,c('tmin_opt','tmin.mean')])
(cor(pwdg.reg[,c('tmin_opt','tmin.mean')]))

fit <- lm(pwdg.reg$clim~pwdg.reg$cwd_opt)
abline(fit)
summary(fit)

plot(pwdg.reg[,c('cwd.mean.x','cwd.mean.y')])
fit <- lm(pwdg.reg$cwd.mean.y~pwdg.reg$cwd.mean.x)
abline(fit)
summary(fit)

plot(pwdg.reg[,c('cwd_opt','cwd.mean.y')])
fit <- lm(pwdg.reg$cwd.mean.y~pwdg.reg$cwd_opt)
abline(fit)
summary(fit)
