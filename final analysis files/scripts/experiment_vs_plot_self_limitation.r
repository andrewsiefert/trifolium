library(tidyverse) 



# read in self-limitation estimates and bootstrap CI's --------------------

# experimental, neighbors present
pres <- read_csv("results/self_limitation_bootstrap_neighbors_present.csv") %>% select(inv, res, self_limitation_pres)

# experimental, neighbors removed
rem <- read_csv("results/self_limitation_bootstrap_neighbors_removed.csv") %>% select(inv, res, self_limitation_rem)

# natural communities
field <- read_csv("results/field_self_limitation_bootstrap.csv") %>% 
  select(inv, res, self.lim) %>%
  rename(self_limitation_field = self.lim) %>%
  mutate_at(vars(inv, res), ~str_replace_all(., "mdon", "mdn"))

# combine neighbors present and natural community data
pres_field <- inner_join(pres, field) %>% na.omit()

# combine neighbors removed and natural community data
rem_field <- inner_join(rem, field) %>% na.omit()


# permutation tests for correlation between self-limitation in exp --------

# function to run permutation test 
perm.test <- function(x, y, n) {
  obs <- cor(x, y, use="complete.obs")
  null <- replicate(n, cor(x, sample(y), use="complete.obs"))
  p <- sum(null>obs)/length(null)
  return(data.frame(obs, p))
}

# test correlation between experiment (neighbors removed) and natural communities
perm.test(x = rem_field$self_limitation_rem, y = rem_field$self_limitation_field, n = 10000)

# test correlation between experiment (neighbors present) and natural communities
perm.test(x = pres_field$self_limitation_pres, y = pres_field$self_limitation_field, n = 10000)



# plot relationships
library(ggplot2)
library(ggthemes)
library(cowplot)

p1 <- ggplot(pres_field, aes(x = self_limitation_pres, y = self_limitation_field)) + 
  geom_hline(yintercept=0, linetype='dashed', col='gray') + 
  geom_vline(xintercept=0, linetype='dashed', col='gray') + 
  stat_smooth(method='lm', se=F) +
  geom_point() + 
  scale_x_continuous(limits=c(-2.25,0.5)) +
  xlab("Stabilization\n(experimental, neighbors present)") +
  ylab("Stabilization\n(natural communities)") +
  theme_few()

p2 <- ggplot(rem_field, aes(x = self_limitation_rem, y = self_limitation_field)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'gray') + 
  geom_vline(xintercept = 0, linetype = 'dashed', col = 'gray') + 
  stat_smooth(method = 'lm', se = F, linetype = 'dashed') +
  geom_point() + 
  xlab("Stabilization\n(experimental, neighbors removed)") +
  ylab("Stabilization\n(natural communities)") +
  theme_few()


cowplot::plot_grid(p2, p1, nrow = 1)

