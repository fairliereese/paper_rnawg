import os
os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
os.environ["STREAMLIT_LOG_LEVEL"] = "debug"


import streamlit as st
import pandas as pd
import tempfile
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cerberus
import swan_vis as swan

SESSION_VARS = ['ca', 'ca_path',
    'sg', 'sg_path',
    'gnames', 'num_cols', 'all_cols']

__version__ = '0.0.1'


d = os.path.dirname(__file__)

st.set_page_config(
    page_title="ENCODE4 Long-read RNA-seq Viewer",
    page_icon=f"{d}/cerberus_logo.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

def init_session_state():
    defaults = {c: None for c in SESSION_VARS}

    for k, v in defaults.items():
        st.session_state.setdefault(k, v)

def main():
    init_session_state()

    st.sidebar.image(f"{d}/cerberus_logo.png", width=300)

    tab_landing, tab_simplex, tab_swan = st.tabs([
        'Home',
        'Simplex view',
        'Swan view'
    ])

    from tabs.landing_tab import render_landing_tab
    from tabs.swan_tab import render_swan_tab
    from tabs.simplex_tab import render_simplex_tab

    with tab_landing:
        render_landing_tab(SESSION_VARS)

    with tab_simplex:
        if st.session_state.ca is None:
            st.info("No data loaded yet.")
        else:
            render_simplex_tab()

    with tab_swan:
        if st.session_state.sg is None:
            st.info("No data loaded yet.")
        else:
            render_swan_tab()

if __name__ == "__main__":
    main()
