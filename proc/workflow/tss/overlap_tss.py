import functools
import numpy as np
import pandas as pd
import pyranges as pr


gr_tss = pr.read_gtf(snakemake.input['gtf']).features.tss()

samples = snakemake.params['samples']
if len(samples):
    express_transcripts = set(
        df_abundance[df_abundance[samples].sum(axis=1) > 0].index)
    gr_tss = gr_tss[gr_tss.transcript_id.isin(express_transcripts)]

pr_cage = pr.PyRanges(
    pd.read_csv(snakemake.input['cage'], sep='\t', header=None).rename(
        columns={0: 'Chromosome', 1: 'Start', 2: 'End', 5: 'Strand'})
)

# pr_cage.End = (pr_cage.Start + pr_cage.End) // 2

# Find overlap between TSS and quantseq polyA cluster
df_overlap = gr_tss.nearest(pr_cage).df
dist = snakemake.params['max_distance']
df_overlap = df_overlap[df_overlap['Distance'] < dist]

# df_overlap = df_overlap[df_overlap['TPM_average'] >
#                         float(snakemake.wildcards['tpm_filter'])]

transcript_overlap = set(df_overlap['transcript_id'])


df_abundance = pd.read_csv(snakemake.input['abundance'], sep='\t') \
                 .rename(columns={'annot_transcript_id': '_transcript_id'}) \
                 .set_index('_transcript_id')

df = gr_tss.df
df['_transcript_id'] = df['transcript_id'].str.split('#').str.get(0)
df = df.set_index('_transcript_id').join(df_abundance, how='inner')

df.loc[df['ISM_subtype'].isin(
    {'None', 'Both'}), 'ISM_subtype'] = 'other'
df['ISM_subtype'] = df['ISM_subtype'].str.lower()

df['transcript_novelty'] = np.where(
    df['transcript_novelty'] == 'ISM',
    df['ISM_subtype'] + ' ' + df['transcript_novelty'],
    df['transcript_novelty']
)


# df_abundance.loc[df_abundance['ISM_subtype'].isin(
#     {'None', 'Both'}), 'ISM_subtype'] = 'other'
# df_abundance['ISM_subtype'] = df_abundance['ISM_subtype'].str.lower()

# df_abundance['transcript_novelty'] = np.where(
#     df_abundance['transcript_novelty'] == 'ISM',
#     df_abundance['ISM_subtype'] + ' ' + df_abundance['transcript_novelty'],
#     df_abundance['transcript_novelty']
# )

df = df.set_index('transcript_id')
transcripts = df.index

# transcripts = df_abundance.index.intersection(gr_tss.transcript_id)

pd.DataFrame({
    'transcript_id': transcripts,
    'support': np.where(transcripts.isin(transcript_overlap), 'yes', 'no'),
    'novelty': df.loc[transcripts, 'transcript_novelty'].tolist()
}).to_csv(snakemake.output['overlap'], index=False)
