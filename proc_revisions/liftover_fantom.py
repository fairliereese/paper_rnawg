import logging
from tqdm import tqdm
import pandas as pd
import pyranges as pr
from liftover import get_lifter


logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO)

df = pr.read_bed(snakemake.input['bed'], as_df=True)
converter = get_lifter('hg19', 'hg38')

__import__("pdb").set_trace()


rows = list()

for i, row in tqdm(df.iterrows()):
    start_lift = converter.query(row.Chromosome, row.Start)
    end_lift = converter.query(row.Chromosome, row.End)

    if (len(start_lift) != 1) \
       or (len(end_lift) != 1):
        logging.warning(f'More than one lift location \n {str(row)}')
        continue

    start_lift = start_lift[0]
    end_lift = end_lift[0]

    if (start_lift[0] != end_lift[0]) \
       or (start_lift[2] != end_lift[2]):
        logging.warning(f'Chromosome or strand dont match \n {str(row)}')
        continue

    orig_len = row['End'] - row['Start']
    lift_len = end_lift[1] - start_lift[1]

    if lift_len < 1:
        logging.warning(f'Invalid coordinates \n {str(row)}')
        continue

    if abs(orig_len - lift_len) > 10:
        logging.warning(f'Match length dont match \n {str(row)}')
        continue

    row['Chromosome'] = start_lift[0]
    row['Start'] = start_lift[1]
    row['End'] = end_lift[1]
    row['Strand'] = start_lift[2]
    rows.append(row)

df_lift = pd.DataFrame(rows)
pr.PyRanges(df_lift).to_bed(snakemake.output['bed'])
