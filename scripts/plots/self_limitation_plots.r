library(tidyverse)
library(ggthemes)
library(lmodel2)
library(cowplot)
library(ggthemes)


# read in self-limitation estimates and bootstrap CI's --------------------

# experimental, neighbors present
pres <- read_csv("results/self_limitation_bootstrap_neighbors_present.csv")

# experimental, neighbors removed
rem <- read_csv("results/self_limitation_bootstrap_neighbors_removed.csv")

# natural communities
field <- read_csv("results/field_self_limitation_bootstrap.csv") %>% 
  rename(self_limitation_field = self.lim, 
         lb_field = lb, 
         ub_field = ub) %>%
  mutate_at(vars(inv, res), ~str_replace_all(., "mdon", "mdn"))

# combine data
dat <- pres %>% left_join(rem) %>% left_join(field) %>%
  rename_all(~str_replace(., "self_limitation", "sl"))
  


# plot experimental self-limitation --------------------------------

dodge <- position_dodge(width=0.6)

dat %>%
  pivot_longer(cols = sl_pres:ub_rem,
               names_to = c(".value", "treatment"),
               names_sep = "\\_") %>%
  ggplot(aes(x=inv, y=sl, group = treatment, color = treatment)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position=dodge, size = 1.25) + 
  geom_errorbar(aes(ymin=lb, ymax=ub), width=0, position=dodge) +
  facet_wrap(~res, scales='free_y', nrow=2) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x='Invader', y='Self-limitation index') +
  theme_few() + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

ggsave("results/plots/self_limitation_exp.svg", width = 8, height = 4.5, units = 'in')


# plot self-limitation, field plots
dat %>% 
  filter(!is.na(sl_field)) %>%
  ggplot(aes(x = inv, y = sl_field)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point(position = dodge, color = '#377eb8') + 
  geom_errorbar(aes(ymin = lb_field, ymax = ub_field), width=0, position=dodge, color = '#377eb8') +
  facet_wrap(~res, scales='free_y', nrow=2) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(x='Invader', y='Self-limitation index') +
  theme_few() + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("results/plots/self_limitation_field.svg", width = 7, height = 4.5, units = 'in')



# self limitation, experimental vs. natural communities

# neighbors removed
fit1 <- lmodel2(sl_field ~ sl_rem, data = dat)
p1 <- ggplot(dat, aes(x = sl_rem, y = sl_field)) + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', alpha = 0.2) +
  geom_abline(aes(intercept = fit1$regression.results[2,2], slope = fit1$regression.results[2,3]), 
              alpha = 0.5, size = 0.75, linetype = 'dashed') +
  geom_point(color = '#1b9e77') + 
  theme_few() +
  labs(x = "Self-limitation index\n(neighbors removed)",
       y = "Self-limitation index\n(natural communities)")

# neighbors present
fit2 <- lmodel2(sl_field ~ sl_pres, data = dat)
p2 <- ggplot(dat, aes(x = sl_pres, y = sl_field)) + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', alpha = 0.2) +
  geom_abline(aes(intercept = fit2$regression.results[2,2], slope = fit2$regression.results[2,3]), alpha = 0.5, size = 0.75) +
  geom_point(color = '#1b9e77') + 
  theme_few() +
  scale_x_continuous(limits = c(-2.2, 0.5)) +
  labs(x = "Self-limitation index\n(neighbors present)",
       y = "Self-limitation index\n(natural communities)")

cowplot::plot_grid(p1, p2, labels = c('a', 'b'), label_size = 16)

ggsave("results/plots/self_limitation_exp_field.svg", 
       width = 7.5, height = 3.5, units = 'in', dpi = 600)
