```bash
snakemake \
    -s Snakefile \
    -j 200 \
    --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
  -s Snakefile \
  -j 200 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=5:00:00" -n

snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A vswarup_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
-s Snakefile \
-j 200 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A vswarup_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
```

```bash
snakemake \
  -s Snakefile \
  -j 60 \
  --latency-wait 120 \
  --use-conda -n
```

```bash
conda activate snakemake_vis
snakemake --forceall --dag | dot -Tpdf > dag.pdf
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
```


## 10/12/23 Debugging

Bad transcript IDs
* 588882
* 601014
* 601063
SELECT * from OBSERVED WHERE transcript_ID IN (588882,601014,601063)

```bash
```
