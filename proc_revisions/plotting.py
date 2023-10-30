import os
import pandas as pd

def get_biosample_colors(species='human'):
    """
    Get colors for each biosample
    """
    d = os.path.dirname(__file__)
    fname = f'{d}/data/{species}/lr/lr_{species}_library_data_summary.tsv'
    df = pd.read_csv(fname, sep='\t')
    df = df[['sample', 'sample_color_hex_code']].drop_duplicates()

    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry['sample']] = entry.sample_color_hex_code
    order = df['sample'].tolist()

    return c_dict, order
