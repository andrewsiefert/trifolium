library(doBy)
library(vegan)
library(tidyverse)


# Prepare data ------------------------------------------------------------

# read in field experiment data
dat <- read.csv("data/2015_2016_field_expt_data.csv")

# calculate mean growth rate per species/site/year/neighbor removal treatment
p= summaryBy(srv + sht_wt ~ sp + site + nr + year, FUN=mean, na.rm=T, keep.names=T, 
               data=dat)
p$growth= p$srv * p$sht_wt
p$growth[is.na(p$growth)]= 0

# convert data to wide format (one row per site/neighbor removal treatment/year)
p= reshape(p[,c(1:4,7)], v.names = "growth", idvar = c("year", "site", "nr"), timevar = "sp", 
             direction = "wide")
names(p)= gsub('growth.', '', names(p))

year= p$year

# prepare environmental data
env= summaryBy(grass + forb + leg + tec + pH + org + bray.p + ca + mg + k + na + fe + 
                   exch.n + sand ~ site, na.rm=T, keep.names=T, data=dat)
rownames(env)= env$site
env= env[,c(2,6:8,12,14,15)]
names(env)= c('grass','pH','org','P','Na','N','sand')
env= env[p$site,]
env$nr= as.numeric(p$nr=='nbr rm')

p= p[,-(1:3)]


# run CCA  ----------------------------------------------------------------

## Neighbors removed

# get environmental, growth, and year data for neighbor removal treatment
env_nr <- env %>% filter(nr == 1) %>% select(-nr)
p_nr <- p %>% filter(env$nr == 1)
year_nr <- year[env$nr == 1]

# run CCA
ord_nr= cca(p_nr ~ . + Condition(year_nr), data=env_nr)

# plot species and site scores
plot(ord_nr, display= c('sp','cn'), scaling=3)

# test overall significance of ordination
anova(ord_nr, permutations = 9999)

# test significance of environmental factors
anova(ord_nr, by='mar', permutations = 9999)

## Neighbors present

# get environmental, growth, and year data for neighbor present treatment
env_n <- env %>% filter(nr == 0) %>% select(-nr)
p_n <- p %>% filter(env$nr == 0)
year_n <- year[env$nr == 0]

# run CCA
ord_n = cca(p_n ~ . + Condition(year_n), data=env_n)

# plot species and site scores
plot(ord_n, display= c('sp','cn'), scaling=3)

# test overall significance of ordination
anova(ord_n, permutations = 9999)

# test significance of environmental factors
anova(ord_n, by='mar', permutations = 9999)  

## calculate pairwise species distances in CCA space
env.dist = dist(ord$CCA$v)


# plot ordination ---------------------------------------------------------
library(ggpubr)
library(ggthemes)
library(cowplot)

nr_plot <- plot(ord_nr, display= c('sp','cn'), scaling=3)
nr_sp <- nr_plot$species %>% as_tibble(rownames = 'species')
nr_env <- (nr_plot$biplot * 1.95) %>% as_tibble(rownames = 'var')

nr_plot <- ggplot(nr_sp) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  geom_text(data = nr_sp, aes(x = CCA1, y = CCA2, label = species, fontface = 'italic'), color = '#377eb8') + 
  geom_segment(data = nr_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(2, 'mm'))) +
  geom_text(data = nr_env, aes(x = CCA1, y = CCA2, label = var)) + 
  theme_few()

n_plot <- plot(ord_n, display= c('sp','cn'), scaling=3)
n_sp <- n_plot$species %>% as_tibble(rownames = 'species')
n_env <- (n_plot$biplot * 1.95) %>% as_tibble(rownames = 'var')

n_plot <- ggplot(n_sp) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  geom_text(data = n_sp, aes(x = CCA1, y = CCA2, label = species, fontface = 'italic'), color = '#377eb8') + 
  geom_segment(data = n_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(2, 'mm'))) +
  geom_text(data = n_env, aes(x = CCA1, y = CCA2, label = var)) + 
  theme_few()

plot_grid(nr_plot, n_plot, labels = c('a', 'b'), label_size = 16)

ggsave("results/plots/cca_plot.pdf",
       height = 4, width = 8.5, units = "in", dpi = 600)
