# Gene View Tab Implementation

# borrowed from here
# https://github.com/mortazavilab/swanViewer/blob/main/tabs/gene_view_tab.py

import streamlit as st
import pandas as pd
from swan_vis.swangraph import SwanGraph
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import plotly.express as px
import io
import itertools
from utils import *

N_T_TPM_META_COLS = 4 # number of metadata columns used for the TPM table

def render_swan_tab():
    """Render the Gene View tab content."""

    geneExpander = st.sidebar.expander(
        "Transcript structure / expression view controls",
        expanded=False
    )

    with geneExpander:
        st.markdown("### Transcript structure / expression view options")
        st.markdown("---")
        if not st.session_state.swan_data_loaded:
            st.info("Load data to enable options")
            return

        sg = st.session_state.sg
        if 'tname' in sg.t_df.columns:
            missing_tn_mask = sg.t_df['tname'].isna()
            if missing_tn_mask.any():
                missing_tids = sg.t_df.index[missing_tn_mask]
                sg.t_df.loc[missing_tn_mask, 'tname'] = missing_tids

        all_genes = sg.t_df['gname'].dropna().astype(str).unique().tolist()
        gene_map_lower = {g.lower(): g for g in all_genes}

        def_gname = None
        if len(all_genes) > 1:
            all_genes = sorted(all_genes)
            if 'ELN' in all_genes:
                def_gname = 'ELN'
            else:
                def_gname = all_genes[0]
            disabled = False
        else:
            disabled = True

        st.markdown("### Transcript structure / expression view options")
        gene_name_input = st.selectbox(label='Gene name', options=all_genes)
        # gene_name_input = check_multiselect(gene_name_input, 'gene name')

        sample_names = sg.adata.obs_names.tolist()
        # Build a display-name mapping: prefer `condition_replicate` when available
        # First look for `sg.obs.data` as requested, otherwise fall back to `sg.adata.obs`.
        obs_df = None
        try:
            if hasattr(sg, 'obs') and hasattr(sg.obs, 'data'):
                obs_df = sg.obs.data
            elif hasattr(sg, 'adata') and hasattr(sg.adata, 'obs'):
                obs_df = sg.adata.obs
        except Exception:
            obs_df = None

        sample_display_map = {s: s for s in sample_names}
        if obs_df is not None and 'condition' in obs_df.columns and 'replicate' in obs_df.columns:
            for s in sample_names:
                try:
                    # Try to align by sample name index; fall back to original name on failure
                    row = obs_df.loc[s]
                    cond = '' if pd.isna(row['condition']) else str(row['condition'])
                    rep = '' if pd.isna(row['replicate']) else str(row['replicate'])
                    if cond and rep:
                        sample_display_map[s] = f"{cond} {rep}"
                    elif cond:
                        sample_display_map[s] = cond
                    elif rep:
                        sample_display_map[s] = rep
                except Exception:
                    # leave the original sample name
                    sample_display_map[s] = s
        st.markdown("---")
        ordered_sample_names = st.multiselect(
            "Select and reorder TPM columns:",
            options=sample_names, default=[],
            format_func=lambda x: sample_display_map.get(x, x),
            help="Click to select, deselect, and drag to reorder the TPM columns."
        )
        # add_all_samples = st.button('Add all')
        # rm_all_samples = st.button('Remove all')

        tpm_display_format = st.radio(
            "TPM Display Format:",
            ["Numbers", "Heatmap"],
            help="Choose how to display TPM values in the tables."
        )
        # Allow the user to choose which obs column to group the violin plot by
        grouping_column = None
        try:
            if obs_df is not None:
                obs_columns = [c for c in obs_df.columns.tolist() if obs_df[c].nunique() > 1]
                if 'condition' in obs_columns:
                    default_group = 'condition'
                elif obs_columns:
                    default_group = obs_columns[0]
                else:
                    default_group = None
                if obs_columns:
                    grouping_column = st.selectbox('Group violin plot by (obs column):', options=['(none)'] + obs_columns, index=(obs_columns.index(default_group) + 1) if default_group in obs_columns else 0)
                    if grouping_column == '(none)':
                        grouping_column = None
        except Exception:
            grouping_column = None

        # actually make the plot
        run = st.button('Run')


    if run:
        if gene_name_input:
            gene_name_lower = gene_name_input.lower()
            if gene_name_lower in gene_map_lower:
                gene_name = gene_map_lower[gene_name_lower]
                transcripts_for_gene = sg.t_df[sg.t_df['gname'] == gene_name]
                gid = transcripts_for_gene['gid'].iloc[0] if not transcripts_for_gene.empty else "N/A"

                # Initialize a variable to hold the colorbar figure
                colorbar_fig = None

                st.header(f"{gene_name} ({gid})")

                # gene-level swan graph
                try:
                    sg.plot_graph(gene_name, indicate_novel=True)
                    st.pyplot(plt.gcf())
                    plt.clf()
                    plt.close('all')
                except Exception as e:
                    st.error(f"An error occurred while generating the gene-level SwanGraph for {gene_name}: {e}")
                    st.exception(e)

                def color_df_col(df, col, color_dict):
                    styles = pd.DataFrame('', index=df.index, columns=df.columns)

                    for idx, row in df.iterrows():
                        color_key = row['color']
                        hex_color = color_dict.get(color_key)
                        if hex_color is not None:
                            styles.loc[idx, col] = f'background-color: {hex_color}'

                    return styles

                # section showing individual details
                try:
                    with st.expander("Show Swan Gene Details"):

                        # transcript info
                        keep_cols = ['tname', 'gid', 'gname',
                            'tss_id', 'ic_id', 'tes_id',
                            'loc_path', 'annotation', 'novelty']
                        keep_cols = [c for c in keep_cols if c in sg.pg.t_df.columns]

                        st.write(f"Transcripts for **{gene_name}**:")
                        st.dataframe(sg.pg.t_df[keep_cols])

                        # vertex info
                        keep_cols = ['chrom', 'coord', 'vertex_id', 'annotation',
                                     'TSS', 'internal', 'TES', 'color']
                        keep_cols = [c for c in keep_cols if c in sg.pg.loc_df.columns]
                        sg.pg.loc_df['hex_color'] = sg.pg.loc_df['color'].map(sg.pg.get_color_dict())
                        st.write(f"Vertices for **{gene_name}**:")
                        st.dataframe(sg.pg.loc_df[keep_cols].style.apply(
                            color_df_col,
                            col='vertex_id',
                            color_dict=sg.pg.get_color_dict(),
                            axis=None),
                            hide_index=True)

                        # edge info
                        keep_cols = ['v1', 'v2', 'strand', 'edge_type', 'edge_id',
                                     'annotation', 'color']
                        keep_cols = [c for c in keep_cols if c in sg.pg.edge_df.columns]
                        sg.pg.loc_df['hex_color'] = sg.pg.edge_df['color'].map(sg.pg.get_color_dict())
                        st.write(f"Edges for **{gene_name}**:")
                        st.dataframe(sg.pg.edge_df[keep_cols].style.apply(
                            color_df_col,
                            col='edge_id',
                            color_dict=sg.pg.get_color_dict(),
                            axis=None),
                            hide_index=True)

                except Exception as e:
                    st.error(f"An error occurred while retreiving Swan gene details for {gene_name}: {e}")
                    st.exception(e)

                # dataframe of requested gene expression across samples
                if ordered_sample_names:

                    st.subheader("Gene-level abundance (TPM)")

                    subset_sg = sg.subset_on_gene_sg(gid=gid, datasets=ordered_sample_names)
                    temp = subset_sg.get_transcript_abundance(kind='tpm')
                    gene_temp = (
                        temp.set_index('tid')
                        .sum(axis=0)
                        .to_frame()
                        .transpose()
                    )
                    temp2 = (
                        gene_temp
                            .style
                            .set_properties(**{'text-align': 'center'})
                            .format('{:.1f}')
                    )
                    st.dataframe(temp2, hide_index=True)

                    # violin plot
                    try:
                        if grouping_column:
                            gene_temp = gene_temp.transpose().reset_index()
                            gene_temp.columns = ['dataset', 'tpm']
                            gene_temp = gene_temp.merge(sg.adata.obs[['dataset', grouping_column]].reset_index(drop=True),
                                       how='left',
                                       on='dataset')

                            fig = px.violin(
                                gene_temp,
                                x=grouping_column,
                                y='tpm',
                                box=True,
                                points='all',
                                hover_data=['dataset'],
                                title=f"Gene-Level TPMs for {gene_name}")
                            plotly_config = {
                                "width": 'stretch'
                            }
                            st.plotly_chart(fig, config=plotly_config)

                    except Exception as e:
                        st.warning(f"Could not render violin plot: {e}")

                    st.subheader("Transcript-level abundance (TPM)")

                    # all values are actual tpms
                    temp.set_index('tid', inplace=True)
                    temp['total_tpm'] = temp.sum(axis=1)

                    # filter out unexpressed
                    temp = temp.loc[temp['total_tpm'] > 0]

                    # sort by total expression
                    temp = temp.sort_values(by='total_tpm', ascending=False)

                    # finally remove sum col.
                    temp.drop('total_tpm', axis=1, inplace=True)

                    # columns
                    col_spec = [2, 2, 3, 3] + [1] * len(temp.columns)
                    header_cols = st.columns(col_spec)
                    header_cols[0].markdown("###### Isoform ID")
                    header_cols[1].markdown("###### Genomic Locus")
                    header_cols[2].markdown("###### Swan Graph View")
                    header_cols[3].markdown("###### Browser View")
                    for i, name in enumerate(temp.columns):
                        display_name = sample_display_map.get(name, name)
                        header_cols[N_T_TPM_META_COLS + i].markdown(f"##### {display_name}")
                    st.divider()

                    # colorbar scaling init
                    min_val, max_val = temp.min().min(), temp.max().max()
                    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
                    cmap = plt.get_cmap("Reds")

                    # get info for each transcript
                    for tid in temp.index.tolist():

                      # row setup
                      row_cols = st.columns(col_spec)

                      # col 0: transcript id
                      with row_cols[0]:
                          st.markdown(f"**`{tid}`**")

                      # col 1: human-readable transcript coords
                      with row_cols[1]:
                          coord_min, coord_max = subset_sg.get_transcript_min_max(tid)
                          chrom = subset_sg.get_transcript_chrom(tid)
                          st.markdown(f'`{chrom}:{coord_min}-{coord_max}`')

                      # col 2: swan graph
                      with row_cols[2]:
                          sg.plot_transcript_path(tid, browser=False, indicate_novel=True)
                          buf = io.BytesIO()
                          plt.savefig(buf, format="png", dpi=150, bbox_inches='tight')
                          st.image(buf)
                          plt.clf()
                          plt.close('all')

                      # col 3: browser plot
                      with row_cols[3]:
                          sg.plot_transcript_path(tid, browser=True)
                          buf = io.BytesIO()
                          plt.savefig(buf, format="png", dpi=150, bbox_inches='tight')
                          st.image(buf)
                          plt.clf()
                          plt.close('all')

                      # cols 4+: TPM (either heatmap or numbers)
                      for i,c in enumerate(temp.columns):
                          tpm_val = temp.loc[tid, c]

                          if tpm_display_format == "Numbers":
                              row_cols[i+N_T_TPM_META_COLS].write(f"{tpm_val:.1f}")
                          elif tpm_display_format == 'Heatmap':
                              if max_val > 0:
                                color_hex = mcolors.to_hex(cmap(norm(tpm_val)))
                                text_color = "white" if norm(tpm_val) > 0.6 else "black"
                                html = f"""<div style="background-color:{color_hex}; color:{text_color}; padding: 10px; border-radius: 5px; text-align: center;">{tpm_val:.1f}</div>"""
                                row_cols[i+N_T_TPM_META_COLS].markdown(html, unsafe_allow_html=True)

                # finally add the colorbar
                if tpm_display_format == "Heatmap":
                    col1, col2, col3 = st.columns((1,5,1))
                    if max_val > 0:
                        fig, ax = plt.subplots(figsize=(3, 0.6))
                        fig.subplots_adjust(bottom=0.1)
                        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
                        cb.set_label('TPM')
                        buf = io.BytesIO()
                        plt.savefig(buf, format="png", dpi=150, bbox_inches='tight')
                        with col3:
                            st.image(buf)
                            plt.clf()
