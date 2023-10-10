```bash
snakemake \
  -s Snakefile \
  -j 60 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
  -s Snakefile \
  -j 200 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=5:00:00" -n
```

```bash
snakemake \
  -s Snakefile \
  -j 60 \
  --latency-wait 120 \
  --use-conda -n
```
