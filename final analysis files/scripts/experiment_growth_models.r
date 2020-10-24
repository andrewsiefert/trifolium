library(tidyverse)
library(cplm)
library(tweedie)
library(statmod)
library(MuMIn)
library(glmmTMB)
library(sjPlot)
library(emmeans)
library(car)
library(DHARMa)

source("scripts/tweedie_functions.r")


## Read in and prepare data for modeling

dat <- read_csv("data/2015_2016_field_expt_data_6sp.csv", guess_max = 9999) %>%
  mutate(growth = srv*sht_wt %>% ifelse(srv == 0, 0, .), 
         year = year - min(year)) %>%
  filter(!is.na(growth), is.na(flag))


## Trifolium-independent model
m0 <- cpglmm(growth ~ year + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*nr + (1|site), link='log', 
             data=dat, optimizer='bobyqa', control=list(max.fun=2e5))

# model diagnostics
qq_tweed(m0)
fit.obs(m0)


# Test species x neighbor removal interaction
m0_d1 = update(m0, .~. - sp:nr)
anova(m0, m0_d1)


# Overall effect of neighbor removal
m0_d2 <- update(m0, .~. - sp:nr - nr)
anova(m0, m0_d2)


# Refit with glmmTMB

m0_tmb <- glmmTMB(growth ~ year*sp + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*nr + (1|site), 
                  family = glmmTMB::tweedie(), data = dat)


#Plot effect of neighbor removal by species
plot_model(m0_tmb, type = "eff", terms = c("sp", "nr"))


# Test neighbor removal effect within species.
sp_nr <- emmeans(m0_tmb, ~nr|sp)
pairs(sp_nr, type = "response")

# Test overall neighbor removal effect
emmeans(m0_tmb, ~nr, type = "response") %>% pairs()



# Trifolium density model -------------------------------------------------

# Fit full model
m1 <- cpglmm(growth ~ year + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*nr*s_tri + (1|site), 
             link='log', data = dat, optimizer='bobyqa', control=list(max.fun=2e4))

# diagnostics
qq_tweed(m1)
fit.obs(m1)


# Test trifolium density x neighbor removal x species interaction
m1.d1= update(m1, .~. - sp:nr:s_tri)
anova(m1, m1.d1) 


## Test two-way interactions.
#Tri density x species

m1.d2.1 <- update(m1.d1, .~. - sp:s_tri)
anova(m1.d1, m1.d2.1)

#Tri density x neighbor removal
m1.d2.2 <- update(m1.d1, .~. - nr:s_tri); anova(m1.d1, m1.d2.2) 


# Refit with glmmTMB
m1_tmb <- glmmTMB(growth ~ year*sp + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*s_tri + sp*nr + (1|site), 
                  family = glmmTMB::tweedie(), data = dat)

# get Trifolium density marginal effect by species
dens_effects <- get_model_data(m1_tmb, 
                               type = "eff", 
                               terms = c("s_tri [0:6]", "sp")) %>% as_tibble()
write_csv(dens_effects, "results/trifolium_density_marginal_effects.csv")


# Test Trifolium density effect within species
tri_sp <- emtrends(m1_tmb, ~sp, "s_tri")
plot(tri_sp) + geom_vline(xintercept = 0)
summary(tri_sp, infer = T)

# Overall Trifolium density effect
emtrends(m1_tmb, specs = ~1, var = "s_tri") %>% summary(infer = T)


# Frequency dependence model ----------------------------------------------

m2 <- cpglmm(growth ~ year + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*s_tri + sp*nr*s_freq + (1|site), 
             link='log', data = dat, optimizer = 'bobyqa', control = list(max.fun=2e5))

# diagnostics
qq_tweed(m2)
fit.obs(m2)


# Test frequency x neighbor removal x species interaction

m2.d1= update(m2, .~. - sp:nr:s_freq); anova(m2, m2.d1)



# Refit with glmmTMB.

m2_tmb <- glmmTMB(growth ~ year*sp + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*s_tri + sp*nr*s_freq + (1|site), 
                  family = glmmTMB::tweedie(), data = dat)


# get frequency effects within neighbor removal and species
freq_effects <- get_model_data(m2_tmb, 
                               type = "eff", 
                               terms = c("s_freq [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]", "nr", "sp")) %>% as_tibble()
write_csv(freq_effects, "results/frequency_marginal_effects.csv")


# Test frequency effect within species and neighbor removal
twt <- emtrends(m2_tmb, ~nr|sp, "s_freq")
plot(twt)
summary(twt, infer = TRUE)

#Test for difference in strength of frequency effect between neighbor removal treatments within species.
plot(pairs(twt))
pairs(twt)

# Test "main" interaction between frequency and neighbor treatment. 
nr_freq2 <- emtrends(m2_tmb, ~nr, "s_freq")

summary(nr_freq2, infer = T)
plot(nr_freq2)
pairs(nr_freq2)



# Species effects model ---------------------------------------------------

# scale species abundance variables
dat_z <- dat %>% mutate_at(vars(contains("s_")), funs(scale))

m3_tmb <- glmmTMB(growth ~ year*sp + bray.p + exch.n + pH + s_forb + s_grass + sand + sp*nr + 
                    sp*s_barb + sp*s_bif + sp*s_fuc + sp*s_gra + sp*s_mac + sp*s_mdn + (1|site),
                  family = glmmTMB::tweedie(link = "log"), data = dat_z)

#Inspect model
m3_tmb %>% simulateResiduals() %>% plot()

Anova(m3_tmb)


### Test species effects within species

emtrends(m3_tmb, specs = "sp", var = "s_barb") %>% summary(infer = T)
emtrends(m3_tmb, specs = "sp", var = "s_bif") %>% summary(infer = T)
emtrends(m3_tmb, specs = "sp", var = "s_fuc") %>% summary(infer = T)
emtrends(m3_tmb, specs = "sp", var = "s_gra") %>% summary(infer = T)
emtrends(m3_tmb, specs = "sp", var = "s_mac") %>% summary(infer = T)
emtrends(m3_tmb, specs = "sp", var = "s_mdn") %>% summary(infer = T)



## Model comparison

AICc_tab <- function(mods, model_names=NULL) {
  blank_names <-  c('m1','m2','m3','m4','m5','m6','m7','m8')

  if(is.null(model_names)) {
    model_names <- blank_names[1:length(mods)]
  }
  
  out <- tibble(AIC = map_dbl(mods, AIC),
  AICc = map_dbl(mods, AICc),
  k = map_dbl(mods, ~sum(Anova(.)$Df)),
  df = map_dbl(mods, ~summary(.)$AICtab['df.resid'])) %>%
  mutate(dAICc = AICc - min(AICc),
  w = exp(-0.5*dAICc),
  weight = w/sum(w),
  model = model_names) %>%
  select(model, k, df, AIC, AICc, dAICc, weight) %>%
  arrange(desc(weight))
  
  return(out)
}

AICc_tab(list(m0_tmb, m1_tmb, m2_tmb, m3_tmb), model_names = c("Trifolium-independent", "Density dependence", "Frequency dependence", "Species effects"))



