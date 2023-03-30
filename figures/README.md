Before generating any of the figures, you will need to download data from the ENCODE portal and process it to generate some additional files used.

There is some overlap between the processing code here and in the raw data processing directory.

```bash
snakemake \
  -s snakemake/human_snakefile.smk \
  -j 10 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END, --time=72:00:00" -n

snakemake \
  -s snakemake/Snakefile \
  -j 10 \
  --latency-wait 120 \
  -n
```
