library(tidyverse)
library(cplm)
library(glmmTMB)
library(emmeans)
library(tweedie)
library(glmmTMB)
library(sjPlot)

source("scripts/tweedie_functions.r")


## Read in and prepare data for modeling
dat <- read_csv("data/2015_2016_field_expt_data_6sp.csv") %>%
  mutate(growth = srv*sht_wt %>% ifelse(srv == 0, 0, .)) %>%
  filter(!is.na(growth))


# Model effects of environmental distance on performance ------------------

### Null model (no environmental distance effects)
m0 = cpglmm(growth ~ year + sp*nr + (1|site), 
            link='log', data=dat, optimizer='bobyqa', control=list(max.fun=2e4))

# model diagnostics
qq_tweed(m0)
fit.obs(m0)
summary(m0)

### Environmental distance model

# Global model
m1 <- cpglmm(growth ~ year + sp*nr*env_dist + (1|site), 
             link='log', data=dat, optimizer='bobyqa', control=list(max.fun=2e4))

# model diagnostics
qq_tweed(m1)
fit.obs(m1)

# test significance of distance to home environment effects
anova(m1, m0)


# Test interactions

m1.d1 <- update(m1, .~. - sp:nr:env_dist); anova(m1, m1.d1)


# refit with glmmTMB
m1_tmb <- glmmTMB(growth ~ year + sp*nr*env_dist + (1|site), 
                  family = tweedie(),
                  data = dat)

# get distance to home environment effects within neighbor removal and species
dist_effects <- get_model_data(m1_tmb, 
                               type = "eff", 
                               terms = c("env_dist [0:7]", "nr", "sp")) %>% as_tibble()
write_csv(dist_effects, "results/env_dist_effects.csv")


# Test environmental distance effect within species and neighbor removal
three_way <- emtrends(m1_tmb, ~sp|nr, "env_dist")
plot(three_way) + geom_vline(aes(xintercept = 0), color = "gray")
summary(three_way, infer = T)


# Test for difference in strength of environmental distance effect between neighbor removal treatments within species.
three_way2 <- emtrends(m1_tmb, ~nr|sp, "env_dist")
plot(pairs(three_way2))
pairs(three_way2)


# Test "main" two-way interaction effect between neighbor removal and environmental distance.
nr_env <- emtrends(m1_tmb, ~nr, "env_dist")
plot(nr_env)
pairs(nr_env)

