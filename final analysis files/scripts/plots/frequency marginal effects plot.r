library(tidyverse)
library(ggthemes)

# read in conspecific frequency marginal effects
effects <- read_csv("results/frequency_marginal_effects.csv") %>%
  # remove values outside of range of natural resident frequencies
  filter((facet=='bif' & x > 0.6)==F, (facet=='mdn' & x > 0.7)==F, (facet=='mdn' & x > 0.9)==F)

ggplot(effects) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, group = group), alpha = 0.1) +
  geom_line(aes(x = x, y = predicted, color = group)) + 
  facet_wrap(~facet, scales = "free_y") +
  theme_few() +
  labs(x = "Conspecific frequency", y = "Mean aboveground biomass (g)")

ggsave("results/plots/frequency_marginal_effects_plot.pdf", 
       width = 7.5, height = 4, units = 'in')


