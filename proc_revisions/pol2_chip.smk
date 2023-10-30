pol2_df = pd.read_csv(expand(config['pol2']['encode_meta'],
                      species='human')[0],
                      sep='\t').set_index('File accession')
pol2_df['species'] = 'human'
temp = pd.read_csv(expand(config['pol2']['encode_meta'],
                      species='mouse')[0],
                      sep='\t').set_index('File accession')
temp['species'] = 'mouse'
pol2_df = pd.concat([pol2_df, temp], axis=0)

wildcard_constraints:
  encid='|'.join([re.escape(x) for x in pol2_df.index.tolist()])

use rule dl_encid_gz_2 as dl_pol2 with:
  params:
      file_ext = 'bed'
  output:
      out = temporary(config['pol2']['bed_gz'])

use rule gunzip as gz_pol2 with:
  input:
      gz = config['pol2']['bed_gz']
  output:
      out = temporary(config['pol2']['bed'])

rule merge_pol2:
  input:
      beds = expand(config['pol2']['bed'],
                    encid=pol2_df.index.tolist(),
                    species=species)
  resources:
      mem_gb = 8,
      threads = 1
  output:
      bed = config['pol2']['merged']
  script:
      "merge_beds.py"
