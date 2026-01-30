```bash
#
mkdir data
wget link https://zenodo.org/records/18299002/files/human_triplets.h5?download=1 -o data/human_triplets.h5
wget link https://zenodo.org/records/18299002/files/mouse_triplets.h5?download=1 -o data/mouse_triplets.h5
wget link https://zenodo.org/records/18299002/files/human_swan_graph.p?download=1 -o data/human_swan.p
wget link https://zenodo.org/records/18299002/files/mouse_swan_graph.p?download=1 -o data/mouse_swan.p

streamlit run data_viewer/main.py

# run docker
docker run --rm \
  -p 8501:8501 \
  -v $(pwd):/app \
  streamlit-test


# interactive docker
docker run --rm -it -v $(pwd):/app -w /app streamlit-test bash

```

```bash
# human
h5=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/data_viewer/human_triplets.h5
sg=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/swan/swan_graph.p
streamlit run data_viewer/main.py --server.maxUploadSize 2000 -- --h5 $h5 --sg $sg

# mouse
h5=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/data_viewer/mouse_triplets.h5
sg=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/mouse/lr/swan/swan_graph.p
streamlit run data_viewer/main.py --server.maxUploadSize 2000 -- --h5 $h5 --sg $sg



```

scartch

```python
import cerberus
import matplotlib.pyplot as plt
h5 = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/data_viewer/human_triplets.h5'
ca = cerberus.read(h5)
ca.plot_simplex(
    subset={'source': 'sample_det'},
    gene='ELN',
    scatter=True,
    density=False,
    density_scale=100,
    log_density=True,
    density_cmap='Purples',
    # marker size stuff
    size='gene_tpm',
    log_size=True,
    size_scale=.75,
    # marker color
    hue='sample',
    sectors=True,
    legend=True
)
plt.show()
```

```python
import swan_vis as swan
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl

sg = swan.read('/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/swan/swan_graph.p')

gid = 'ENSG00000067225'
subset_sg = sg.subset_on_gene_sg(gid=gid, datasets=['ovary_1_1', 'ovary_2_1', 'ovary_3_1'])
temp = subset_sg.get_transcript_abundance(kind='tpm')
temp.set_index('tid', inplace=True)
min_val, max_val = temp.min().min(), temp.max().max()
temp = temp.loc[temp.sum(axis=1) > 0]
N_T_TPM_META_COLS = 4
tid = 'ENSG00000049540[1,156,3]'
for i,c in enumerate(temp.columns):
  tpm_val = temp.loc[tid, c]
  row_cols[i+N_T_TPM_META_COLS].write(f"{tpm_val:.1f}")


```

# build docker
```bash
docker build -t streamlit-test .

```

```bash
source /opt/conda/etc/profile.d/conda.sh
conda activate app_env

```
```python
import pandas as pd
h5='data/human_triplets.h5'
key='tss'
df = pd.read_hdf(h5, key=key)
```

Making multiarch image
```bash
docker buildx create --use
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -t ghcr.io/fairliereese/encode4-viewer:latest \
    --push .
```

pushing to GHCR
```bash
echo $GHCR_TOKEN
docker login ghcr.io
docker tag streamlit-test ghcr.io/fairliereese/encode4-viewer:latest
docker push ghcr.io/fairliereese/encode4-viewer:latest

# to run using it
docker pull ghcr.io/fairliereese/encode4-viewer:latest
docker run -p 8501:8501 ghcr.io/fairliereese/encode4-viewer:latest

```

<!-- Important) Make it public -->

<!--
By default, GHCR images are private, even if the repo is public.

Go to:

GitHub → Packages → your image

Package settings → Change visibility → Public -->

to see whats in latest
```bash
docker buildx imagetools inspect ghcr.io/fairliereese/encode4-viewer:latest

```
