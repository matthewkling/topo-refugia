## compare PWD and regional means

names <- read.csv('data/names.csv',as.is=T)
names

pwdgam <- read.csv('data/pwd_gam1.csv',as.is=T)
pwdgam
pwdgam$sci.name <- names$Sci.name[match(pwdgam$species,names$Hyp.name)]

pwd <- read.csv('data/pwd_niche_means.csv',as.is=T)
pwd

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
plot(pwdg.reg$cwd_opt,pwdg.reg$cwd.mean)
plot(pwdg.reg$cwd_opt,pwdg.reg$south.mean)
plot(pwdg.reg$cwd_opt,pwdg.reg$clim)
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
