```bash
h5=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/cerberus/sample_cerberus_triplets.h5
sg=/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/swan/swan_graph.p
streamlit run data_viewer/main.py --server.maxUploadSize 2000 -- --h5 $h5 --sg $sg

```

scartch

```python
import cerberus
import matplotlib.pyplot as plt
h5 = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/cerberus/sample_cerberus_triplets.h5'
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
sg = swan.read('/Users/fairliereese/Documents/programming/mortazavi_lab/data/paper_rnawg/proc_revisions/data/human/lr/swan/swan_graph.p')
```
