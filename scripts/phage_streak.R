# look at phage streaks

library(ggplot2)
library(dplyr)
library(tidyr)
library(MicrobioUoE)

# read in data
d <- read.csv('data/20171206_phage_streaks.csv', stringsAsFactors = FALSE)
meta <- read.csv('data/sample_metadata.csv', stringsAsFactors = FALSE)

d <- merge(d, meta, by = 'sample')

# sample 30 has 2 lacz clones in, only 18 wt clones
d <- gather(d, 'phage', 'num_resist', 2:5) %>%
  mutate(., prop_resist = ifelse(sample == 30, num_resist/18, num_resist/20)) %>%
  group_by(., sample, type, treatment, phage) %>%
  summarise(., sd = sd(prop_resist),
            prop_resist = mean(prop_resist)) %>%
  ungroup()

d_mean <- group_by(d, phage) %>%
  summarise(., prop_resist = mean(prop_resist)) %>%
  ungroup() %>%
  mutate(., group = 1)
   
# plot
p1 <- ggplot(d, aes(phage, prop_resist)) +
  geom_pretty_boxplot(col = 'black', fill = 'black') +
  geom_point(position = position_jitter(width = 0.2, height = 0), shape = 21, fill = 'white') +
  facet_wrap(~ type) +
  ggtitle('Proportional resistance of bacteria to ancestral, contemporary and other community phage') +
  xlab('Phage streak') +
  ylab('Proportion of resistant bacteria')


p2 <- ggplot(d, aes(phage, prop_resist, group = sample)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ sample) +
  ggtitle('Proportional resistance of bacteria to ancestral, contemporary and other community phage') +
  xlab('Phage streak') +
  ylab('Proportion of resistant bacteria') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(d, aes(phage, prop_resist)) +
  geom_line(aes(group = sample), alpha = 0.5) +
  geom_line(aes(group = group), data = d_mean, size = 2) +
  ggtitle('Proportional resistance of bacteria to ancestral, contemporary and other community phage') +
  xlab('Phage streak') +
  ylab('Proportion of resistant bacteria') 

# save plots
pdf('figures/phage_streak.pdf', width = 10, height = 7)
print(p1)
print(p2)
print(p3)
dev.off()
  