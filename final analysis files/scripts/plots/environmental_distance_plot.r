library(tidyverse)
library(ggthemes)

effects <- read_csv("results/env_dist_effects.csv")

ggplot(effects) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, group = group), alpha = 0.1) +
  geom_line(aes(x = x, y = predicted, color = group)) + 
  facet_wrap(~facet, scales = "free") +
  theme_few() +
  #coord_cartesian(xlim = c(0.1, 6), ylim = c(0, 0.3)) +
  labs(x = "Distance to home environment", y = "Mean aboveground biomass (g)")

ggsave("results/plots/environment_distance_plot.pdf", 
       width = 7.5, height = 4, units = 'in')
