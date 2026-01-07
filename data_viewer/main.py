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

__version__ = '0.0.1'

st.set_page_config(
    page_title="ENCODE4 Long-read RNA-seq Viewer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def get_cli_args():
    """Parses command-line arguments for file paths."""
    parser = argparse.ArgumentParser(description='ENCODE4 LR-RNA-seq viewer')
    parser.add_argument('--h5', type=str, help='Path to a saved Cerberus h5 file to load')
    args, _ = parser.parse_known_args()
    return args

def load_data(fname):
    """
    Load CerberusAnnotation into streamlit space
    """

    if fname and st.session_state.get('ca') is None:
        try:
            # check if we've already loaded this session file
            if fname in st.session_state.loaded_session_cache:
                ca = st.session_state.loaded_session_cache[fname]
            else:
                ca = cerberus.read(fname)
                print(f"Loaded CerberusAnnotation from {fname}")
                # Cache the loaded session
                st.session_state.loaded_session_cache[fname] = ca
            st.session_state.ca = ca

            st.session_state.data_loaded = True

            # mark that this session was loaded
            st.session_state.loaded_from_disk = True
            st.session_state.loaded_from_disk_path = fname

            # store a user-friendly filename to display in UI
            try:
                # for CLI loads, display full absolute path to be explicit
                st.session_state.loaded_from_disk_filename = os.path.abspath(fname)
            except Exception:
                st.session_state.loaded_from_disk_filename = fname

        except Exception as e:
            st.error(f"Failed to load file: {e}")
            st.stop()

# main
def main():

    cli_args = get_cli_args()

    ##### INIT SESSION VARS
    st.session_state.data_loaded = False

    # cache for loaded session files to avoid reloading
    if not hasattr(st.session_state, 'loaded_session_cache'):
        st.session_state.loaded_session_cache = {}

    # add cache for UI-loaded session files
    if not hasattr(st.session_state, 'ui_loaded_session_cache'):
        st.session_state.ui_loaded_session_cache = {}

    ##### DATA LOADING
    # load data from cmd line
    load_data(cli_args.h5)

    # if no CerberusAnnotation was loaded, add option to upload file on the sidebar
    if not st.session_state.data_loaded:
        uploaded_file = st.sidebar.file_uploader("Upload CerberusAnnotation", type=["h5"])
        if uploaded_file is not None:
            import tempfile
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5") as tmp:
                tmp.write(uploaded_file.getbuffer())
                file_to_load = tmp.name
                load_data(file_to_load)

    # creating tabs for each visualization
    tab_gene_view, = st.tabs(["Gene View"])

    from tabs.gene_simplex_tab import render_gene_simplex_tab

    with tab_gene_view:
        render_gene_simplex_tab()

    #
    # # testing
    # if st.session_state.data_loaded:
    #     st.write(st.session_state.ca.tss.head(1).source)

    # # if we have a file path, load
    # if file_to_load is not None and st.session_state.ca is None:

    # # testing -- display something from the ca
    # ca = cerberus.read(cli_args.h5)
    # print(f"Loaded CerberusAnnotation from {cli_args.h5}")
    # st.session_state.loaded_session_cache[fname] = ca
    #
    # st.session_state.ca = ca

    st.markdown("---")



if __name__ == "__main__":
    main()
