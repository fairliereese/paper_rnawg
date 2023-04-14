rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggsci)
library(forcats)
setwd('/gpfs/gibbs/pi/gerstein/bb.shared.projects/RNAWG/S7_C')
## different PSI threshold (0.05/0.95, 0.1/0.9, 0.25/0.75)
summary <- read.csv('lr_human_library_data_summary.tsv', sep = '\t', header = TRUE) %>% select(sample, sample_display) %>% distinct()
AF_combined <- read.table('cerberus_suppa/0.25_0.75/tss_combined.tsv', sep = '\t', header = T) %>% na.omit() %>%
  dplyr::filter(suppa == TRUE | cerberus == TRUE) %>% mutate(type = ifelse((suppa == TRUE & cerberus == TRUE), 'Consensus', ifelse(suppa == TRUE, 'Only suppa', 'Only cerberus')))
AF_combined <- merge(AF_combined, summary, by = 'sample')

AL_combined <- read.table('cerberus_suppa/0.25_0.75/tes_combined.tsv', sep = '\t', header = T) %>% na.omit() %>%
  dplyr::filter(suppa == TRUE | cerberus == TRUE) %>% mutate(type = ifelse((suppa == TRUE & cerberus == TRUE), 'Consensus', ifelse(suppa == TRUE, 'Only suppa', 'Only cerberus')))
AL_combined <- merge(AL_combined, summary, by = 'sample')

p1 <- ggplot(AF_combined, aes(x = fct_infreq(sample_display), fill=factor(type, levels = c('Only suppa', 'Consensus','Only cerberus')))) + geom_bar() + 
  scale_fill_manual(values = c('Consensus' = '#BCBDDC', 'Only suppa' = '#EFEDF5', 'Only cerberus' = '#756BB1')) + 
  theme_classic() + scale_x_discrete(guide = guide_axis(angle = 45)) + labs(x = '', y = 'AS genes (#)', title = 'AF (SUPPA) vs. TSS (Cerberus)') + theme(legend.title=element_blank(), legend.position="bottom")
p2 <- ggplot(AL_combined, aes(x = fct_infreq(sample_display), fill=factor(type, levels = c('Only suppa', 'Consensus','Only cerberus')))) + geom_bar() + 
  scale_fill_manual(values = c('Consensus' = '#BCBDDC', 'Only suppa' = '#EFEDF5', 'Only cerberus' = '#756BB1')) +
  theme_classic() + scale_x_discrete(guide = guide_axis(angle = 45)) + labs(x = '', y = 'AS genes (#)', title = 'AL (SUPPA) vs. TES (Cerberus)') + theme(legend.title=element_blank(), legend.position="bottom")

library(patchwork)
p <- p1 / p2 + plot_layout(guides = 'collect') & theme(legend.position='bottom')
ggsave('cerberus_suppa_psi_0.25_0.75.pdf', width = 10, height = 8)

