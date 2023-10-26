from tqdm import tqdm
import functools
import pandas as pd
import pyranges as pr


pr.PyRanges(
    pd.concat([
        pd.read_csv(i, sep='\t', header=None)
        .rename(columns={0: 'Chromosome', 1: 'Start', 2: 'End'})
        for i in tqdm(snakemake.input['beds'])
    ])
).merge(strand=False).drop().to_bed(snakemake.output['bed'])
