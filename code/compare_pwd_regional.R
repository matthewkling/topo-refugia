## compare PWD and regional means

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
