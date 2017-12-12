# quick counts check

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(MicrobioUoE)
library(readr)
library(lubridate)
library(purrr)

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

# filter out Time 2
d <- filter(d, date != '2017-11-14')

# load in phage data
d_phage <- read.csv('data/20171120_phage_counts.csv', stringsAsFactors = FALSE) %>%
  mutate(., PFUs_ml = calc_PFU_ml(phage_count, phage_dil, 10),
         date = ymd('20171120'))

# bind together
d <- left_join(d, select(d_phage, sample, date, PFUs_ml), by = c('sample', 'date'))

# add T0 abundance
d_T0 <- meta %>%
  mutate(., date = ymd('20171019'),
         CFUs_ml = ifelse(type == 'WT', calc_CFU_ml(64, 10^-4, 30), calc_CFU_ml(18, 10^-4, 30)),
         PFUs_ml = ifelse(treatment == 'phage', 4e+8, 0))

# bind data with T0
d <- bind_rows(select(d, names(d_T0)), d_T0)

# gather together
d2 <- gather(d, phage_bact, count, c(CFUs_ml, PFUs_ml)) %>%
  mutate(phage_bact = ifelse(phage_bact == 'CFUs_ml', 'bacteria', 'phage'))

# plot of counts
ggplot(filter(d2, phage_bact == 'bacteria'), aes(interaction(type, treatment), log10(count), fill = type, col = type)) +
  geom_pretty_boxplot() +
  geom_point(position = position_jitter(width = 0.15), shape = 21, fill = 'white') +
  xlab('') +
  ylab('log10 CFUs per mL') +
  scale_color_manual(values = c('blue', 'black', 'blue', 'black')) +
  scale_fill_manual(values = c('blue', 'black', 'blue', 'black')) +
  theme(legend.position = 'none') +
  ylim(c(4, 8)) +
  facet_wrap(~ date)

ggsave('figures/bact_abund_time2.pdf', last_plot(), height = 5, width = 12)

# means
d_mean <- group_by(d2, date, type, treatment, phage_bact) %>%
  summarise(sd = sd(log10(count)),
            count = mean(count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(., count = ifelse(is.nan(count), NA, count))

ggplot(filter(d2, !is.na(count) & count != 0), aes(date, log10(count), col = type)) +
  geom_line(aes(group = interaction(sample, phage_bact)), alpha = 0.1) +
  geom_point(aes(shape = phage_bact, group = interaction(treatment, type, phage_bact)), size = 3, filter(d_mean, !is.na(count) & count != 0), position = position_dodge(width = 1.3)) +
  geom_linerange(aes(x = date, ymin = log10(count) - sd, ymax = log10(count) + sd, group = interaction(treatment, type, phage_bact)), position = position_dodge(width = 1.3),  filter(d_mean, !is.na(count) & count != 0)) +
  geom_line(aes(group = interaction(treatment, type, phage_bact)), position = position_dodge(width = 1.3), filter(d_mean, !is.na(count))) +
  scale_color_manual(values = c('blue', 'black')) +
  xlab('Date') +
  ylab('log10 count') +
  ggtitle('Bacteria & Phage abundance through time',
          subtitle = 'Means and standard deviations and each microcosms trajectory') +
  ylim(2.5, 9) +
  facet_wrap(~ treatment)

ggsave('figures/bact_abund_time.pdf', last_plot(), height = 5, width = 11)
