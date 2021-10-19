library(lme4)
library(tidyverse)
library(reshape2)


dat= read.csv("data/2015_2016_plot_tri_density.csv")
dat= subset(dat, dat$type %in% c('grid','off-grid','plot'))
dat$plot= factor(dat$plot); dat$transect= factor(dat$transect)
d15= subset(dat, year=='2015')[,5:14]
d16= subset(dat, year=='2016')[,6:14]

d= merge(d15, d16, by='plot', all=T, suffixes=c('.15','.16'))


barb <- d %>%
  filter(barb.15 > 0, barb.16 > 0) %>%
  mutate(barb.rgr = log(barb.16) - log(barb.15))

bif <- d %>%
  filter(bif.15 > 0, bif.16 > 0) %>%
  mutate(bif.rgr = log(bif.16) - log(bif.15))

fuc <- d %>%
  filter(fuc.15 > 0, fuc.16 > 0) %>%
  mutate(fuc.rgr = log(fuc.16) - log(fuc.15))

gra <- d %>%
  filter(gra.15 > 0, gra.16 > 0) %>%
  mutate(gra.rgr = log(gra.16) - log(gra.15))

mac <- d %>%
  filter(mac.15 > 0, mac.16 > 0) %>%
  mutate(mac.rgr = log(mac.16) - log(mac.15))

mdon <- d %>%
  filter(mdon.15 > 0, mdon.16 > 0) %>%
  mutate(mdon.rgr = log(mdon.16) - log(mdon.15))


# fit year-to-year growth models ------------------------------------------

# barb
barb.fit <- lmer(barb.rgr ~ barb.15 + bif.15 + fuc.15 + gra.15 + mac.15 + mdon.15 + 
                   (1|transect), data=barb)
summary(barb.fit)
plot(barb.fit)
qqnorm(resid(barb.fit)); qqline(resid(barb.fit))
plot(predict(barb.fit), barb$barb.rgr); abline(a=0, b=1)
barb.ci <- confint(barb.fit, method = 'boot', type = "parametric")

# bif
bif.fit <- lmer(bif.rgr ~ barb.15 + bif.15 + gra.15 + mac.15 + mdon.15 + 
                  (1|transect), data=bif)
summary(bif.fit)
plot(bif.fit)
qqnorm(resid(bif.fit)); qqline(resid(bif.fit))
plot(predict(bif.fit), bif$bif.rgr); abline(a=0, b=1)
bif.ci <- confint(bif.fit, method='boot', type="parametric")

# fuc
fuc.fit <- lmer(fuc.rgr ~ barb.15 + fuc.15 + mac.15 + mdon.15 + 
                  (1|transect), data=fuc)
summary(fuc.fit)
plot(fuc.fit)
qqnorm(resid(fuc.fit)); qqline(resid(fuc.fit))
plot(predict(fuc.fit), fuc$fuc.rgr); abline(a=0, b=1)
fuc.ci <- confint(fuc.fit, method='boot', type="parametric")

# gra
gra.fit <- lmer(gra.rgr ~ barb.15 + bif.15 + gra.15 + mac.15 + mdon.15 + 
                  (1|transect), data=gra)
summary(gra.fit)
plot(gra.fit)
qqnorm(resid(gra.fit)); qqline(resid(gra.fit))
plot(predict(gra.fit), gra$gra.rgr); abline(a=0, b=1)
gra.ci= confint(gra.fit, method='boot', type="parametric")

# mac
mac.fit <- lmer(mac.rgr ~ barb.15 + bif.15 + fuc.15 + gra.15 + mac.15 + mdon.15 + 
                  (1|transect), data=mac)
summary(mac.fit)
plot(mac.fit)
qqnorm(resid(mac.fit)); qqline(resid(mac.fit))
plot(predict(mac.fit), mac$mac.rgr); abline(a=0, b=1)
mac.ci <- confint(mac.fit, method='boot', type="parametric")

# mdon
mdon.fit <- lmer(mdon.rgr ~ barb.15 + bif.15 + fuc.15 + gra.15 + mac.15 + mdon.15 + 
                   (1|transect), data=mdon)
summary(mdon.fit)
plot(mdon.fit)
qqnorm(resid(mdon.fit)); qqline(resid(mdon.fit))
plot(predict(mdon.fit), mdon$mdon.rgr); abline(a=0, b=1)
mdon.ci <- confint(mdon.fit, method='boot', type="parametric")


# combine and plot interaction coefficients -------------------------------

# combine coefficients
fit.list <- list(barb.fit, bif.fit, fuc.fit, gra.fit, mac.fit, mdon.fit)
eff <- names(fixef(barb.fit))[-1]

rgr.coef <- t(sapply(fit.list, function(i) sapply(eff, function(j) fixef(i)[j])))
colnames(rgr.coef) = rownames(rgr.coef) = c('barb','bif','fuc','gra','mac','mdon')

# combine CI's
ci.list <- list(barb.ci, bif.ci, fuc.ci, gra.ci, mac.ci, mdon.ci)

rgr.lb <- t(sapply(ci.list, function(i) sapply(eff, function(j) 
  tryCatch(i[j,1], error= function(e) NA))))
rgr.ub <- t(sapply(ci.list, function(i) sapply(eff, function(j) 
  tryCatch(i[j,2], error= function(e) NA))))

# combine coef estimates, lower bounds, upper bounds
rgr.out <- as.data.frame(melt(rgr.coef))
rgr.out$lb <- melt(rgr.lb)$value
rgr.out$ub <- melt(rgr.ub)$value
names(rgr.out)[1:3]= c("inv","res","coef")

# plot coefficients
dodge <- position_dodge(width=0.5)
ggplot(rgr.out, aes(x=res, y=coef, col=inv)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=dodge) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=0, position=dodge) +
  theme_bw()  

# export interaction coefficients
write_csv(rgr.out, "results/field_interaction_coefficients.csv")

# Self-limitation ---------------------------------------------------------

## calculate self-limitation for each invader-resident pair
#   self lim = res effect on invader - res effect on itself
#   positive values = self limitation

sl <- matrix(NA, nrow=6, ncol=6)
rownames(sl) = colnames(sl) = c('barb','bif','fuc','gra','mac','mdon')
for(i in rownames(sl)) {
  for(j in colnames(sl)) {
    a.jj= rgr.coef[j,j]
    a.ij= rgr.coef[i,j]
    sl[i,j]= a.ij - a.jj
  }
}

## get bootsrapped confidence intervals on interaction coefficients
get.coef <- function(m) fixef(m)[-1]
coef.boot <- lapply(fit.list, function(i) 
  bootMer(i, get.coef, nsim = 10000, type='parametric', parallel='multicore')$t)

names(coef.boot)= c('barb','bif','fuc','gra','mac','mdon')

sl.lb = sl.ub = matrix(NA, nrow=6, ncol=6)
rownames(sl.lb) = colnames(sl.lb) = rownames(sl.ub) = colnames(sl.ub) = names(coef.boot)
for(i in names(coef.boot)) {
  for(j in names(coef.boot)) {
    a.jj <- tryCatch(coef.boot[[j]][,paste(j,'.15',sep='')], error= function(e) rep(NA, 10000))
    a.ij <- tryCatch(coef.boot[[i]][,paste(j,'.15',sep='')], error= function(e) rep(NA, 10000))
    sl.ij <- a.ij - a.jj
    sl.lb[i,j] <- quantile(sl.ij, 0.025, na.rm=T)
    sl.ub[i,j] <- quantile(sl.ij, 0.975, na.rm=T)
  }
}


sl.out <- melt(sl)
sl.out$lb <- melt(sl.lb)$value
sl.out$ub <- melt(sl.ub)$value
names(sl.out)[1:3] <- c('inv','res','self.lim')
sl.out <- subset(sl.out, inv!=res)

# plot self-limitation with confidence intervals
dodge <- position_dodge(width=0.5)
ggplot(sl.out, aes(x=res, y=self.lim, col=inv)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=dodge) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=0, position=dodge) +
  theme_bw()  


# export results ----------------------------------------------------------

write_csv(sl.out, "results/field_self_limitation_bootstrap.csv")
