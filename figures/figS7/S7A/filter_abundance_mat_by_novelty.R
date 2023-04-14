setwd('/gpfs/gibbs/pi/gerstein/bb.shared.projects/RNAWG/S7_A')

library(dplyr)
# read the novelty file
df <- read.csv('cerberus_transcript_novelty.tsv', sep = '\t')
tx_id <- df %>% filter(source == 'lapa') %>% pull(transcript_id) %>% unique()
# filter the abundance matrix
abundance <- read.csv('suppa/cerberus_filtered_abundance.tsv', sep = '\t')

df2 <- abundance %>% filter(annot_transcript_id %in% tx_id)
write.table(df2, 'suppa/cerberus_filtered_abundance.novelty.tsv', sep = '\t', quote = F, row.names = F)
