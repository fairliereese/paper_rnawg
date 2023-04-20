rm(list = ls())
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ggalluvial)
library(ggsankey)
setwd('/gpfs/gibbs/pi/gerstein/bb.shared.projects/RNAWG/S7_B')
# sankey plot of lapa transcripts and gtex transcripts
## using new transcript table
transcripts <- read.table('cerberus_transcript_novelty.tsv', sep = '\t', header = T, check.names = F, comment.char = '')
transcripts$transcript_id <- gsub(',','_', transcripts$transcript_id)

df <- transcripts %>% 
  filter(source != 'v40' & source != 'v29') %>% 
  select(transcript_id, source, ic_novelty) %>% distinct() %>%
  pivot_wider(names_from = source, values_from = ic_novelty, values_fill = 'Missing') %>%
  group_by(lapa, gtex) %>% summarise(n = n()) %>% as.data.frame()
# relabel
df$lapa <- factor(df$lapa, levels =c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced', 'Missing'))
df$gtex <- factor(df$gtex, levels =c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced', 'Missing'))

cols <- c('Known' = '#009E73', 'ISM' = '#0072B2',
          'NIC' = '#D55E00', 'NNC' = '#E69F00', 
          'Unspliced' = '#F0E442', 'Missing' = 'gray60')

tot <- ggplot(df, aes(y = n, axis1 = lapa, axis2 = gtex,
                      fill = factor(after_stat(stratum), levels = c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced','Missing')))) + 
  stat_alluvium(geom = "flow", lode.guidance = "forward", width = 0.1) + 
  stat_stratum(aes(fill = factor(after_stat(stratum), levels = c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced','Missing'))), size = 0.1, width = 0.1) + 
  scale_x_discrete(limits = c('Observed', 'GTEx'), expand = rep(0.1, 2)) +
  scale_fill_manual(values = cols) + labs(x = 'Source', y = 'Count', title = 'All Transcripts') +
  theme_alluvial(base_size = 14) + theme(legend.title=element_blank()) 
# ggsave('sankey_plot_all_transcripts.pdf', tot, width = 5, height = 5)

## sankey plot of lapa transcripts and gtex transcripts for each event
sankeyPlot_event <- function(transcripts, lapa_ioe_path, gtex_ioe_path, mapping, event_name){
  lapa <- read.table(lapa_ioe_path, header = T)
  lapa_events <- unlist(lapply(lapa$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
  gtex <- read.table(gtex_ioe_path, header = T)
  gtex_events <- unlist(lapply(gtex$alternative_transcripts, function(x) str_split(x, ',')[[1]]))
  gtex_events_mapping <- subset(mapping, mapping$original_transcript_id %in% gtex_events)
  events <- c(lapa_events, gtex_events_mapping$transcript_id)
  
  
  df <- transcripts %>% 
    filter(source != 'v40' & source != 'v29' & transcript_id %in% events) %>% 
    select(transcript_id, source, ic_novelty) %>% distinct() %>%
    pivot_wider(names_from = source, values_from = ic_novelty, values_fill = 'Missing') %>%
    group_by(lapa, gtex) %>% summarise(n = n()) %>% as.data.frame()
  
  df$lapa <- factor(df$lapa, levels =c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced', 'Missing'))
  df$gtex <- factor(df$gtex, levels =c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced', 'Missing'))
  
  cols <- c('Known' = '#009E73', 'ISM' = '#0072B2',
            'NIC' = '#D55E00', 'NNC' = '#E69F00', 
            'Unspliced' = '#F0E442', 'Missing' = 'gray60')
  p <- ggplot(df, aes(y = n, axis1 = lapa, axis2 = gtex,
                      fill = factor(after_stat(stratum), levels = c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced','Missing')))) + 
    stat_alluvium(geom = "flow", lode.guidance = "forward", width = 0.1) + 
    stat_stratum(aes(fill = factor(after_stat(stratum), levels = c('Known', 'ISM', 'NIC', 'NNC', 'Unspliced','Missing'))), size = 0.1, width = 0.1) +
    scale_x_discrete(limits = c('Observed', 'GTEx'), expand = rep(0.1, 2)) +
    scale_fill_manual(values = cols) + labs(x = 'Source', y = 'Count', title = event_name) +
    theme_alluvial(base_size = 14) + theme(legend.title=element_blank()) 
  return(p)
}

gtex_to_cerberus <- read.table('gtex_to_cerberus_id.tsv', header = T)
gtex_to_cerberus$transcript_id <- gsub(',','_', gtex_to_cerberus$transcript_id)

AF <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_AF_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_AF_strict.ioe', 
                       gtex_to_cerberus, 'AF')

AL <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_AL_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_AL_strict.ioe', 
                       gtex_to_cerberus, 'AL')

A3 <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_A3_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_A3_strict.ioe', 
                       gtex_to_cerberus, 'A3')


A5 <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_A5_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_A5_strict.ioe', 
                       gtex_to_cerberus, 'A5')

MX <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_MX_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_MX_strict.ioe', 
                       gtex_to_cerberus, 'MX')

RI <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_RI_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_RI_strict.ioe', 
                       gtex_to_cerberus, 'RI')

SE <- sankeyPlot_event(transcripts, 
                       lapa_ioe_path = 'suppa/localEvents/cerberus.events_SE_strict.ioe', 
                       gtex_ioe_path ='suppa/localEvents_GTEx/flair.events_SE_strict.ioe', 
                       gtex_to_cerberus, 'SE')

library(patchwork)
p <- (tot + AF + AL + A3 + plot_layout(guides = 'collect', ncol = 4)) / (A5 + MX + RI + SE + plot_layout(guides = 'collect', ncol = 4)) + plot_layout(guides = 'collect')
ggsave('sankey_plot.pdf', p, width = 16, height = 8)
