from tqdm import tqdm
import functools
import pandas as pd
import pyranges as pr


pr.PyRanges(
    pd.concat([
        pd.read_csv(i, sep='\t', header=None)
        .rename(columns={0: 'Chromosome', 1: 'Start', 2: 'End', 5: 'Strand'})
        for i in tqdm(snakemake.input['beds'])
    ])
).merge(strand=True).drop().to_bed(snakemake.output['bed'])
