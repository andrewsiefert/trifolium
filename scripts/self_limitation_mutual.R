### Check for "mutual" self-limitation in experiment (neighbors present or removed) and field plots

library(tidyverse)

# read in self-limitation values
pres <- read_csv("results/self_limitation_bootstrap_neighbors_present.csv")
rem <- read_csv("results/self_limitation_bootstrap_neighbors_removed.csv")
field <- read_csv("results/field_self_limitation_bootstrap.csv") %>%
  rename(self_limitation_field = self.lim, 
         lb_field = lb,
         ub_field = ub) %>%
  mutate_at(vars(inv, res), ~str_replace_all(., "mdon", "mdn"))

# combine data and create flags for significant self-limitation
d <- left_join(pres, rem) %>% 
  left_join(field) %>%
  mutate(pres_sig = ifelse(lb_pres > 0, 1, ifelse(ub_pres < 0, -1, 0)),
         rem_sig = ifelse(lb_rem > 0, 1, ifelse(ub_rem < 0, -1, 0)),
         field_sig = ifelse(lb_field > 0, 1, ifelse(ub_field < 0, -1, 0)))

# get list of species pairs
spp <- unique(d$res)
pairs <- combn(spp, 2) %>% t()


### Check for mutual self-limitation in experiment, neighbors present

# rearrange data to get self-limitation value and statistical significance in either direction  for each species pair
pres <- NULL

for(i in 1:nrow(pairs)) {
  sp1 <- pairs[i, 1]
  sp2 <- pairs[i, 2]
  sl1 <- d$self_limitation_pres[d$inv == sp1 & d$res == sp2]
  sl2 <- d$self_limitation_pres[d$inv == sp2 & d$res == sp1]
  sig1 <- d$pres_sig[d$inv == sp1 & d$res == sp2]
  sig2 <- d$pres_sig[d$inv == sp2 & d$res == sp1]
  pres <- rbind(pres, data.frame(sp1, sp2, sl1, sl2, sig1, sig2))
}

# flag species pairs where both species are self-limiting as resident with the other as invader
pres <- pres %>% mutate(both = ifelse(sl1 > 0 & sl2 > 0, 1, 0),
                        both_sig = sig1 * sig2)


### Check for mutual self-limitation in experiment, neighbors removed

# rearrange data to get self-limitation value and statistical significance in either direction  for each species pair
rem <- NULL

for(i in 1:nrow(pairs)) {
  sp1 <- pairs[i, 1]
  sp2 <- pairs[i, 2]
  sl1 <- d$self_limitation_rem[d$inv == sp1 & d$res == sp2]
  sl2 <- d$self_limitation_rem[d$inv == sp2 & d$res == sp1]
  sig1 <- d$rem_sig[d$inv == sp1 & d$res == sp2]
  sig2 <- d$rem_sig[d$inv == sp2 & d$res == sp1]
  rem <- rbind(rem, data.frame(sp1, sp2, sl1, sl2, sig1, sig2))
}

# flag species pairs where both species are self-limiting as resident with the other as invader
rem <- rem %>% mutate(both = ifelse(sl1 >0 & sl2 > 0, 1, 0),
                      both_sig = sig1 * sig2)


### Check for mutual self-limitation in field plots

# rearrange data to get self-limitation value and statistical significance in either direction  for each species pair
field <- NULL

for(i in 1:nrow(pairs)) {
  sp1 <- pairs[i, 1]
  sp2 <- pairs[i, 2]
  sl1 <- d$self_limitation_field[d$inv == sp1 & d$res == sp2]
  sl2 <- d$self_limitation_field[d$inv == sp2 & d$res == sp1]
  sig1 <- d$field_sig[d$inv == sp1 & d$res == sp2]
  sig2 <- d$field_sig[d$inv == sp2 & d$res == sp1]
  field <- rbind(field, data.frame(sp1, sp2, sl1, sl2, sig1, sig2))
}

# flag species pairs where both species are self-limiting as resident with the other as invader
field <- field %>% mutate(both = ifelse(sl1 >0 & sl2 > 0, 1, 0),
                        both_sig = sig1 * sig2)

