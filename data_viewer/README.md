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
