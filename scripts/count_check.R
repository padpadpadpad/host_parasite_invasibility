# quick counts check

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(MicrobioUoE)
library(readr)
library(lubridate)

# list count data
files <- list.files('data', pattern = 'sampling_counts', full.names = TRUE)

# load in data
d <- purrr::map_df(files, read_csv, .id = 'num')
files <- data.frame(file = basename(files), num = 1:3, stringsAsFactors = FALSE) %>%
  mutate(date = ymd(parse_number(files)))
d <- merge(d, select(files, -file), by = 'num')

# load in meta data
meta <- read.csv('data/sample_metadata.csv', stringsAsFactors = FALSE)

# add column for CFUs / ml
d <- mutate(d, CFUs_ml = calc_CFU_ml(bact_count, bact_dil, 30)) %>%
  merge(., meta, by = 'sample')

# plot of counts
ggplot(d, aes(interaction(type, treatment), log10(CFUs_ml), fill = type, col = type)) +
  geom_pretty_boxplot() +
  geom_point(position = position_jitter(width = 0.15), shape = 21, fill = 'white') +
  xlab('') +
  ylab('log10 CFUs per mL') +
  scale_color_manual(values = c('black', 'blue', 'black', 'blue')) +
  scale_fill_manual(values = c('black', 'blue', 'black', 'blue')) +
  theme(legend.position = 'none') +
  ylim(c(4, 8)) +
  facet_wrap(~ date)

ggsave('figures/bact_abund_time2.pdf', last_plot(), height = 5, width = 12)

# means
d_mean <- group_by(d, date, type, treatment) %>%
  summarise(sd = sd(log10(CFUs_ml)),
            CFUs_ml = mean(CFUs_ml)) %>%
  ungroup()

ggplot(d, aes(date, log10(CFUs_ml), col = type)) +
  geom_line(aes(group = sample), alpha = 0.1) +
  geom_point(aes(shape = treatment, group = interaction(treatment, type)), size = 3, d_mean, position = position_dodge(width = 1.3)) +
  geom_linerange(aes(x = date, ymin = log10(CFUs_ml) - sd, ymax = log10(CFUs_ml) + sd, group = interaction(treatment, type)), position = position_dodge(width = 1.3),  d_mean) +
  geom_line(aes(group = interaction(treatment, type)), position = position_dodge(width = 1.3), d_mean) +
  scale_color_manual(values = c('black', 'blue')) +
  xlab('Date') +
  ylab('log10 CFUs/ml') +
  ggtitle('Bacteria abundance through time',
          subtitle = 'Means and standard deviations and each microcosms trajectory')

ggsave('figures/bact_abund_time.pdf', last_plot(), height = 5, width = 8)
