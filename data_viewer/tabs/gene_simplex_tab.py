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
            density = st.checkbox('Density (takes a while to compute)')
            legend = st.checkbox('Legend')

            # marker size
            num_cols = ca.triplets.select_dtypes(include="number").columns.tolist()
            marker_size = st.multiselect(label='Marker size', options=sorted(num_cols))
            marker_size = check_multiselect(marker_size, 'marker size')

            marker_color = st.multiselect(label='Marker color', options=sorted(ca.triplets.columns))
            marker_color = check_multiselect(marker_color, 'marker color')


            # "advanced"
            adv_expander = st.expander("Advanced")
            with adv_expander:

                # sector boundaries
                tss_thresh = st.slider(
                              "TSS-high threshold",
                              0.0, 1.0, 0.5)
                tes_thresh = st.slider(
                              "TES-high threshold",
                              0.0, 1.0, 0.5)
                spl_thresh = st.slider(
                              "Splicing-high threshold",
                              0.0, 1.0, 0.5)

            # log continuous marker colors


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

                legend_scatter = legend_density = False
                if legend and scatter: legend_scatter = True
                if legend and density: legend_density = True

                ca.plot_simplex(
                    subset={'source': triplet_set},
                    gene=gname,
                    scatter=scatter,

                    density=density,
                    density_scale=100,
                    log_density=True,
                    density_cmap='Purples',

                    # marker size stuff
                    size=marker_size,
                    log_size=True,

                    size_scale=.75,

                    # marker color
                    hue=marker_color,

                    sectors=True,
                    legend=legend_scatter,
                    density_cbar=legend_density,

                    # sector boundaries
                    sect_alpha=tss_thresh,
                    sect_beta=tes_thresh,
                    sect_gamma=spl_thresh
                )

                st.pyplot(plt.gcf())
                plt.clf()
                plt.close('all')
