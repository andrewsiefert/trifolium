library(cplm)
library(tweedie)
library(statmod)
library(MuMIn)
library(reshape2)
library(ggplot2)
library(compiler)
library(tidyverse)

source("scripts/tweedie_functions.r")


### read and prepare data

dat <- read.csv("data/2015_2016_field_expt_data.csv")
dat <- dat[-(which(dat$flag == 1)), ]

# convert shoot weight to 0 for individuals that did not survive
dat$growth <- dat$srv * dat$sht_wt
dat$growth[dat$srv == 0] <-  0

# split data by neighbor removal treatment
dat.split <- split(dat, dat$nr)

# Function to run parametric bootstrap and extract species interaction coefficients
coef.boot <- function(m, n=1000) {
  
  # generate n random draws from fitted model
  y.pred <- predict(m, type='response')  # mu
  y.rand <- replicate(n, rtweedie(length(y.pred), mu = y.pred, power = m$p, phi = m$phi))
  
  # refit model to each bootstrap sample
  eff.boot <- as.data.frame(t(sapply(1:n, function(n) 
    update(m, data = data.frame(growth = y.rand[,n], m@frame[,-1]))$fixef)))
  
  # extract species interaction coefficients from each bootstrap model
  eff.names <- list('s_barb', 's_bif', 's_fuc', 's_gra', 's_mac', 's_mdn')
  sp.eff <- list(NA,NA,NA,NA,NA,NA)
  for(i in 1:length(eff.names)) {
    eff <- eff.boot[,grep(eff.names[[i]], names(eff.boot))]
    sp.eff[[i]] <- data.frame(eff[,1], eff[,1] + eff[,-1])
    names(sp.eff[[i]]) <- c('barb', 'bif', 'fuc', 'gra', 'mac', 'mdn', 'mlum', 'will')
  }
  
  names(sp.eff) <- c('barb', 'bif', 'fuc', 'gra', 'mac', 'mdn')
  return(sp.eff)
}

### Get bootstrapped species interactions coefficients for each NR treament  

## Neighbors present
# model to bootstrap
m.pres <- cpglmm(growth ~ year + bray.p + exch.n + pH + s_forb + s_grass + sand + 
            sp*s_barb + sp*s_bif + sp*s_fuc + sp*s_gra +
            sp*s_mac + sp*s_mdn + (1|site), link = 'log', 
          data = dat.split[[1]], optimizer = 'bobyqa', control = list(max.fun=2e5))

# run bootstrap
coef.boot.pres <- coef.boot(m.pres, n = 10000)

## Neighbors removed
# model to bootstrap
m.rem <- cpglmm(growth ~ year + bray.p + exch.n + pH + s_forb + s_grass + sand + 
            sp*s_barb + sp*s_bif + sp*s_fuc + sp*s_gra +
            sp*s_mac + sp*s_mdn + (1|site), link = 'log', 
          data = dat.split[[2]], optimizer = 'bobyqa', control = list(max.fun=2e5))

# run bootstrap
coef.boot.rem <- coef.boot(m.rem, n = 10000)

############################################################
## Self-limitation - neighbors present
# a.ij - a.jj (positive values mean self-limitation)

# created empty data frames to hold self-limitation extimates and upper and lower bounds of 95% confidence interval
inv.pres.m <- as.data.frame(matrix(NA, nrow = 8, ncol = 6))
names(inv.pres.m) <- levels(dat$sp)[1:6]
rownames(inv.pres.m) <- levels(dat$sp)
inv.pres.ub <- inv.pres.lb <- inv.pres.m

# For each resident-invader combination, calculate self-limitation index using species interaction coefficients from 
# each bootstrap sample and get upper and lower bounds of 95% confidence interval

for(i in levels(dat$sp)) {
  for(j in names(coef.boot.pres)) {
    
    a.jj <- coef.boot.pres[[j]][,j]  # get resident effect on itself
    a.ij = coef.boot.pres[[j]][,i]   # get resident effect on invader
    
    inv = a.ij - a.jj                # calculate self-limitation
    
    inv.pres.m[i,j] = quantile(inv, 0.5)     # get median self-limitation value from bootstrap samples
    inv.pres.ub[i,j] = quantile(inv, 0.975)  # get bootstrapped upper bound
    inv.pres.lb[i,j] = quantile(inv, 0.025)  # get bootstrapped lower bound
  }
}

# combine self-limitation estimates, lower and upper bounds, and reshape data
inv.pres <- melt(inv.pres.m)
inv.pres$res <- rep(levels(dat$sp)[1:6], each=8)
inv.pres$inv <- rep(levels(dat$sp), 6)
inv.pres$self_limitation_pres <- inv.pres$value
inv.pres <- inv.pres[, c(4,3,5)]
inv.pres$lb_pres <- melt(inv.pres.lb)$value
inv.pres$ub_pres <- melt(inv.pres.ub)$value
inv.pres <- subset(inv.pres, inv != res)

# plot self-limitation values and confidence intervals
dodge <- position_dodge(width=0.4)
ggplot(inv.pres, aes(x = res, y = self_limitation_pres, col = inv)) + 
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_point(position = dodge) + 
  geom_errorbar(aes(ymin = lb_pres, ymax = ub_pres), width=0, position=dodge) +
  theme_bw()

write_csv(inv.pres, "results/self_limitation_bootstrap_neighbors_present.csv")

############################################################
## Self-limitation - neighbors removed

# created empty data frames to hold self-limitation extimates and upper and lower bounds of 95% confidence interval
inv.rem.m <- as.data.frame(matrix(NA, nrow = 8, ncol = 6))
names(inv.rem.m) <- levels(dat$sp)[1:6]
rownames(inv.rem.m) <- levels(dat$sp)
inv.rem.ub = inv.rem.lb = inv.rem.m


for(i in levels(dat$sp)) {
  for(j in names(coef.boot.rem)) {
    
    a.jj <- coef.boot.rem[[j]][,j]  # get resident effect on itself
    a.ij <- coef.boot.rem[[j]][,i]  # get resident effect on invader
    
    inv <- a.ij - a.jj              # calculate self-limitation
    
    inv.rem.m[i,j] <- quantile(inv, 0.5)     # get median self-limitation value from bootstrap samples
    inv.rem.ub[i,j] <- quantile(inv, 0.975)  # get bootstrapped upper bound
    inv.rem.lb[i,j] <- quantile(inv, 0.025)  # get bootstrapped lower bound
  }
}

# combine self-limitation estimates, lower and upper bounds, and reshape data

inv.rem <- melt(inv.rem.m)
inv.rem$res <- rep(levels(dat$sp)[1:6], each = 8)
inv.rem$inv <- rep(levels(dat$sp), 6)
inv.rem$self_limitation_rem <- inv.rem$value
inv.rem <- inv.rem[, c(4,3,5)]
inv.rem$lb_rem <- melt(inv.rem.lb)$value
inv.rem$ub_rem <- melt(inv.rem.ub)$value
inv.rem <- subset(inv.rem, inv != res)

dodge <- position_dodge(width=0.4)
ggplot(inv.rem, aes(x=res, y=self_limitation, col=inv)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=dodge) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=0, position=dodge) +
  theme_bw()

write_csv(inv.rem, "results/self_limitation_bootstrap_neighbors_removed.csv")


### Count resident-invader pairs w/ significant feedback
sum(inv.pres$lb_pres >= 0) # significant self-limitation
sum(inv.pres$ub_pres <= 0) # significant invader-limitation

sum(inv.rem$lb_rem >= 0) # significant self-limitation
sum(inv.rem$ub_rem <= 0) # significant invader-limitation


### Permutation test for difference between neighbor removal treatments

perm.test <- function(x, y, n = 1000) {
  
  # combine self-limitation values from neighbor removed and neighbor present treatments
  d <- cbind(x,y)
  
  # randomly permutate the data n times
  d.perm <- lapply(1:n, function(z) t(apply(d, 1, function(i) sample(i,2))))
  
  # generate null distribution of Wilcox test statistic for difference between treatments
  null.dist <- sapply(d.perm, function(i) wilcox.test(i[,1], i[,2], paired=T)$statistic)
  
  # calculate observed test statistic
  obs <- wilcox.test(x,y, paired=T)$statistic
  
  # calculate proportion of values in null distribution that are less than observed value
  p.lt <- sum(null.dist<obs)/n
  
  # calculate proportion of values in null distribution that are greater than observed value
  p.gt <- sum(null.dist>obs)/n
  
  return(data.frame(p.lt, p.gt))
}

# run permutation test 
perm.test(inv.rem$self_limitation_rem, inv.pres$self_limitation_pres, n = 10000)


### Count resident-invader pairs w/ significant feedback
sum(inv.pres$lb>=0) # significant self-limitation
sum(inv.pres$ub<=0) # significant invader-limitation

sum(inv.rem$lb>=0) # significant self-limitation
sum(inv.rem$ub<=0) # significant invader-limitation




