import streamlit as st
import cerberus
import matplotlib.pyplot as plt
import tempfile
import swan_vis as swan
import gc

from utils import *

def reset_data(SESSION_VARS):
    for k in SESSION_VARS:
        st.session_state[k] = None

def load_swan_data(species, SESSION_VARS):
    # fname = f'/data/{species.lower()}_swan.p'
    fname = f'data/{species.lower()}_swan.p'

    # only reload if this is a different file
    if st.session_state.get("sg_path") == fname:
        return

    # reset
    reset_data(['sg'])
    gc.collect()

    # actually load in the data
    try:
        sg = load_swan_from_path(fname)
        st.session_state.sg = sg
        st.session_state.sg_path = fname

    except Exception as e:
        st.error(f"Failed to load data: {e}")
        st.stop()

def load_cerberus_data(species, SESSION_VARS):
    fname = f'data/{species.lower()}_triplets.h5'
    # fname = f'/data/{species.lower()}_triplets.h5'

    # only reload if this is a different file
    if st.session_state.get("ca_path") == fname:
        return

    # reset
    reset_data(['ca'])
    gc.collect()

    # actually load in the data
    try:
        ca = load_cerberus_from_path(fname)
        st.session_state.ca = ca
        st.session_state.ca_path = fname

        # for each triplet set, get the valid gnames,
        # marker colors, and marker sizes
        st.session_state.gnames = {}
        st.session_state.num_cols = {}
        st.session_state.all_cols = {}

        for triplet_set in ca.triplets.source.unique():
            temp = ca.triplets.loc[ca.triplets.source==triplet_set]
            st.session_state.gnames[triplet_set] = sorted(
                temp.gname
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            st.session_state.num_cols[triplet_set] = sorted(
                temp
                .select_dtypes(include="number")
                .dropna(axis=1, how="all")
                .columns
                .tolist()
            )

            st.session_state.all_cols[triplet_set] = sorted(
                temp
                .dropna(axis=1, how='all')
                .columns
                .tolist()
            )

    except Exception as e:
        st.error(f"Failed to load data: {e}")
        st.stop()

@st.cache_resource(show_spinner=True)
def load_cerberus_from_path(path: str):
    return cerberus.read(path)

@st.cache_resource(show_spinner=True)
def load_swan_from_path(path: str):
    return swan.read(path)

def triplets_info():
    # information
    ex = st.expander("Information about gene triplets")

    with ex:
        st.info("""
        ### What are gene triplets?

        Gene triplets summarize **relative isoform usage** for genes with three major transcript categories.
        All triplets shown here are computed **only from polyA genes**.

        Different *triplet sets* correspond to different subsets of transcripts used in the computation.
        Use this guide to understand what each option represents.
        """)

        st.markdown("## Human triplet sets")

        st.markdown("""
        **`v40`**
        All annotated transcripts from **GENCODE v40** genes.

        **`gtex`**
        Triplets computed from the **GTEx long-read RNA-seq dataset**
        (Glinos et al., Nature 2022).

        **`lapa`**
        All transcripts detected **post-LAPA**, unfiltered
        (i.e. before transcript-level filtering in the pipeline).

        **`v29`**
        All annotated transcripts from **GENCODE v29** genes.

        **`obs_det`**
        All **observed transcripts** detected in the dataset
        (aggregated across samples).

        **`obs_major`**
        Only the **major (most highly expressed)** observed transcript per gene
        (aggregated across samples).

        **`sample_det`**
        Transcripts **detected per condition**, computed at the condition level.

        **`sample_major`**
        Major transcripts per gene **within each conditions**.

        **`obs_mm_det`**
        Observed transcripts from **human conditions matched to mouse conditions**
        (for cross-species comparisons).

        **`obs_mm_major`**
        Major observed transcripts from **human conditions matched to mouse conditions**.

        **`all`**
        Aggregate of all transcript sets.
        """)

        st.markdown("## Mouse triplet sets")

        st.markdown("""
        **`vM25`**
        All annotated transcripts from **GENCODE vM25** genes.

        **`vM21`**
        All annotated transcripts from **GENCODE vM21** genes.

        **`lapa`**
        All transcripts detected **post-LAPA**, unfiltered.

        **`obs_det`**
        All **observed transcripts** detected in the dataset
        (aggregated across samples).

        **`obs_major`**
        Only the **major observed transcript** per gene
        (aggregated across samples).

        **`sample_det`**
        Transcripts **detected per conditions**, computed at the conditions level.

        **`sample_major`**
        Major transcripts per gene **within each conditions**.

        **`all`**
        Aggregate of transcript sets.
        """)

def render_landing_tab(SESSION_VARS):
    """
    Choose and load data
    """

    st.markdown('## Welcome to the ENCODE4 long-read RNA-seq data viewer')

    st.markdown('#### Select species to begin')
    species = st.selectbox(
        label='Species',
        options=['Human', 'Mouse']
    )
    st.species = species
    go_species = st.button('Go')

    if go_species:
        load_cerberus_data(species, SESSION_VARS)
        load_swan_data(species, SESSION_VARS)

    triplets_info()
