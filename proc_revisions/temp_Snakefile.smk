import pandas as pd
import os
import sys
import numpy as np

p = os.getcwd()
sys.path.append(p)

from sm_utils import *
from utils import *
from plotting import *

############# config stuff
configfile: 'config.yml'

datasets_per_talon_run = config['params']['talon']['datasets_per_run']
species=['human', 'mouse']
end_modes = ['tss', 'tes']

# lr stuff
lr_df = process_lr_metadata(config['lr']['meta'],
                            species,
                            datasets_per_talon_run)

# external snakefiles
include: 'download.smk'
include: 'samtools.smk'
include: 'spike_refs.smk'
include: 'refs.smk'
include: 'talon.smk'
include: 'lapa.smk'
include: 'cerberus.smk'
include: 'gtex.smk'
include: 'encode_tss.smk'
include: 'ccre.smk'
include: 'fantom_cage.smk'
include: 'polyasite_atlas.smk'
include: 'pas.smk'
include: 'lrgasp_cage.smk'
include: 'procap.smk'
include: 'swan.smk'
include: 'pol2_chip.smk'
include: 'analysis.smk'
include: 'mane.smk'
include: 'read_lens.smk'

############# snakemake settings stuff
wildcard_constraints:
    dataset='|'.join([re.escape(x) for x in lr_df.dataset.tolist()]),
    species='|'.join([re.escape(x) for x in lr_df.species.unique().tolist()]),
    talon_run='|'.join([re.escape(x) for x in lr_df.talon_run.astype(str).unique().tolist()])

############ rules
rule all:
    input:
        rules.all_read_lens.input
