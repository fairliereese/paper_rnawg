import streamlit as st
import cerberus
import matplotlib.pyplot as plt
from utils import *

def render_gene_simplex_tab():
    """
    Original simplex dorito plot
    """

    ca = st.session_state.get('ca', None)

    if ca is not None:

        # layout of page
        col1, col2 = st.columns(2)

        with col2:
            st.markdown("### Simplex View Options")

            # triplet set picker
            triplet_set = st.multiselect(
                label='Triplet set',
                options=ca.triplets.source.unique().tolist(),
                default="sample_det"
            )
            triplet_set = check_multiselect(triplet_set, 'triplet subset')

            # gene picker -- limit to genes in this triplet set
            gnames = (ca.triplets.loc[ca.triplets.source==triplet_set]
                .gname
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            def_gname = None
            if len(gnames) > 1:
                gnames = sorted(gnames)
                if 'ELN' in gnames:
                    def_gname = 'ELN'
                else:
                    def_gname = gnames[0]
                disabled = False
            else:
                disabled = True

            gname = st.multiselect(label='Gene name', options=gnames, default=def_gname, disabled=disabled)
            gname = check_multiselect(gname, 'gene name')

            scatter = st.checkbox("Scatter")
            density = st.checkbox('Density')

            # actually make the plot
            run = st.button('Plot')

        with col1:

            if run:

                # if gname and triplet_set:
                if gname:
                    gid = ca.triplets.loc[ca.triplets.gname==gname].gid.values[0]
                    st.header(f"{triplet_set} triplets for: `{gname}` (Ensembl ID: `{gid}`)")
                else:
                    st.header(f"{triplet_set} triplets")
                # c_dict, order = get_biosample_colors()
                ca.plot_simplex(
                    subset={'source': triplet_set},
                    gene=gname,
                    scatter=scatter,
                    density=density,
                    size_scale=.75,
                    hue='sample',
                    density_cmap='Purples'
                )

                st.pyplot(plt.gcf())
                plt.clf()
                plt.close('all')
