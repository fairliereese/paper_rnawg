import pyranges as pr

species = ['human']
output_types = ['bidirectional_peaks', 'unidirectional_peaks']

def get_encid_from_dataset(dataset, meta, file_format):
    m = {'label_bam': 'ENCODE_alignments_id',
         'bam': 'ENCODE_unfiltered_alignments_id',
         'fastq': 'ENCODE_reads_id'}
    if file_format in list(m.keys()):
        id = meta.loc[meta.dataset == dataset, m[file_format]].values[0]
    else:
        id = meta.loc[meta.dataset == dataset, file_format].values[0]
    return id

# procap meta
procap_meta = pd.read_csv(expand(config['procap']['lib_meta'],
                                 species='human')[0],
                                 sep='\t')

wildcard_constraints:
  dataset='|'.join([re.escape(x) for x in procap_meta.dataset.tolist()]),
  output_type='|'.join([re.escape(x) for x in output_types])

rule dl_procap:
  resources:
      mem_gb = 32,
      threads = 1
  params:
      encid = lambda wc:get_encid_from_dataset(wc.dataset,
                                               procap_meta,
                                               wc.output_type)
  output:
      bed = config['procap']['bed_gz']
  shell:
      "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.bed.gz -O {output.bed}"

use rule gunzip as gz_procap with:
  input:
      gz = config['procap']['bed_gz']
  output:
      out = config['procap']['bed']

rule reformat_procap_uni:
  input:
      bed = config['procap']['bed']
  resources:
      mem_gb = 10,
      threads = 1
  output:
      bed = config['procap']['bed_formatted']
  run:
      df = pd.read_csv(input.bed, sep='\t',
                       header=None,
                       names=['Chromosome', 'Start',
                              'End', 'idk1', 'idk2',
                              'Strand'], usecols=[i for i in range(6)])
      df.drop(['idk1', 'idk2'], axis=1, inplace=True)
      pr.PyRanges(df).to_bed(output.bed)

rule merge_procap_cage:
  input:
      beds = expand(config['procap']['bed_formatted'],
             species=species,
             dataset=procap_meta.dataset.tolist(),
             output_type=output_types)
  resources:
      mem_gb = 8,
      threads = 1
  output:
      bed = config['lrgasp_cage']['merged']
  script:
      "merge_beds.py"

rule all_procap:
    input:
        rules.merge_procap_cage.output,
        # expand(rules.merge_procap_cage.output,
        #        species=species,
        #        dataset=procap_meta.dataset.tolist(),
        #        output_type=output_types)
