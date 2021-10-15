library(tidyverse)
library(ggthemes)


# read in Trifolium density marginal effects estimates
effects <- read_csv("results/trifolium_density_marginal_effects.csv") %>%
  mutate(x = x / (pi*.075^2)) # convert resident count to density per m2

# plot effects
ggplot(effects) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1, fill = '#2ca25f') +
  geom_line(aes(x = x, y = predicted), color = '#2ca25f') + 
  facet_wrap(~group, scales = "free") +
  theme_few() +
  #coord_cartesian(xlim = c(0.1, 6), ylim = c(0, 0.3)) +
  labs(x = "Trifolium density", y = "Mean aboveground biomass (g)")

# save plot
ggsave("results/plots/trifolium_density_plot.pdf", 
       width = 6.5, height = 4, units = 'in')
