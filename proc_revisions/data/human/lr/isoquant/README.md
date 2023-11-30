## IsoQuant
/gpfs/gibbs/project/gerstein/yj329/conda_envs/py3/bin/isoquant.py -d pacbio_ccs --bam bam/ENCFF373TKM.sorted.bam bam/ENCFF388HXU.sorted.bam bam/ENCFF985LGZ.sorted.bam -r GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -g gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz --complete_genedb --stranded forward --output wtc11 -t 32 --sqanti_output --check_canonical


## Reference-based
transcript_grouped_tpm.tsv - TSV file with reference transcript expression in TPM;
gene_grouped_tpm.tsv - TSV file with reference gene expression in TPM;
novel_vs_known.SQANTI-like.tsv - discovered novel transcripts vs reference transcripts

## Transcript discovery
transcript_models.gtf - GTF file with discovered expressed transcript (both known and novel transcripts);
transcript_model_grouped_tpm.tsv - expression of discovered transcripts models in TPM (corresponds to transcript_models.gtf)

## GTF
extended_annotation.gtf - GTF file with the entire reference annotation plus all discovered novel transcripts