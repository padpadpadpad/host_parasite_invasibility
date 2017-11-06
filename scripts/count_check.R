# quick counts check

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(MicrobioUoE)

# functions


# load in data
d <- read.csv('data/20171101_sampling_counts.csv', stringsAsFactors = FALSE)
meta <- read.csv('data/sample_metadata.csv', stringsAsFactors = FALSE)

# add column for CFUs / ml
d <- mutate(d, CFUs_ml = calc_CFU_ml(bact_count, bact_dil, 30)) %>%
  merge(., meta, by = 'sample')

# plot these 

geom_pretty_boxplot <- function(...){
  list(
    ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.75), ...),
    ggplot2::stat_summary(geom = 'crossbar', position = ggplot2::position_dodge(width = 0.75), fatten = 0, color = 'white', width = 0.4, fun.data = function(x){ return(c(y = stats::median(x), ymin = stats::median(x), ymax = stats::median(x)))})
  )
}

# plot of counts
ggplot(d, aes(interaction(type, treatment), log10(CFUs_ml), fill = type, col = type)) +
  geom_pretty_boxplot() +
  geom_point(position = position_jitter(width = 0.15), shape = 21, fill = 'white') +
  xlab('') +
  ylab('log10 CFUs per mL') +
  scale_color_manual(values = c('black', 'blue', 'black', 'blue')) +
  scale_fill_manual(values = c('black', 'blue', 'black', 'blue')) +
  theme(legend.position = 'none') +
  ylim(c(5, 8))
