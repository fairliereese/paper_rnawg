# 1. novelties
# 2. filt ab
# 3. output file

# setwd('/gpfs/gibbs/pi/gerstein/bb.shared.projects/RNAWG/S7_A')

library(optparse)

option_list <- list(
  make_option(c('--novelties', type='character', default=NULL, help='')),
  make_option(c('--filtab', type='character', default=NULL, help='')),
  make_option(c('--ofile', type='character', default=NULL, help=''))
)

# parse args
opt_parser <- OptionParser(option_list=option_list)

opt <- parse_args(opt_parser)

# set variables from args
nov <- opt$novelties
filt_ab <- opt$filtab
ofile <- opt$ofile


library(dplyr)
# read the novelty file
df <- read.csv(nov, sep = '\t')
tx_id <- df %>% filter(source == 'lapa') %>% pull(transcript_id) %>% unique()
# filter the abundance matrix
abundance <- read.csv(filt_ab, sep = '\t')

df2 <- abundance %>% filter(annot_transcript_id %in% tx_id)
write.table(df2, ofile, sep = '\t', quote = F, row.names = F)
