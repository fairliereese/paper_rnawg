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

__version__ = '0.0.1'


d = os.path.dirname(__file__)

st.set_page_config(
    page_title="ENCODE4 Long-read RNA-seq Viewer",
    page_icon=f"{d}/cerberus_logo.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

def get_cli_args():
    """Parses command-line arguments for file paths."""
    parser = argparse.ArgumentParser(description='ENCODE4 LR-RNA-seq viewer')
    parser.add_argument('--h5', type=str, help='Path to a saved Cerberus h5 file')
    parser.add_argument('--sg', type=str, help='Path to a saved SwanGraph (pickle) file')
    args, _ = parser.parse_known_args()
    return args

def load_swan_data(fname):
    """
    Load SwanGraph into streamlit space
    """

    if fname and st.session_state.get('sg') is None:
        try:
            # check if we've already loaded this session file
            if fname in st.session_state.swan_loaded_session_cache:
                sg = st.session_state.swan_loaded_session_cache[fname]
            else:
                sg = swan.read(fname)
                print(f"Loaded SwanGraph from {fname}")
                # Cache the loaded session
                st.session_state.swan_loaded_session_cache[fname] = sg
            st.session_state.sg = sg

            st.session_state.swan_data_loaded = True

            # mark that this session was loaded
            st.session_state.swan_loaded_from_disk = True
            st.session_state.swan_loaded_from_disk_path = fname

            # store a user-friendly filename to display in UI
            try:
                # for CLI loads, display full absolute path to be explicit
                st.session_state.swan_loaded_from_disk_filename = os.path.abspath(fname)
            except Exception:
                st.session_state.swan_loaded_from_disk_filename = fname

        except Exception as e:
            st.error(f"Failed to load file: {e}")
            st.stop()

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
    st.session_state.swan_data_loaded = False

    # cache for loaded session files to avoid reloading
    if not hasattr(st.session_state, 'loaded_session_cache'):
        st.session_state.loaded_session_cache = {}
    if not hasattr(st.session_state, 'swan_loaded_session_cache'):
        st.session_state.swan_loaded_session_cache = {}


    # add cache for UI-loaded session files
    if not hasattr(st.session_state, 'ui_loaded_session_cache'):
        st.session_state.ui_loaded_session_cache = {}
    if not hasattr(st.session_state, 'swan_ui_loaded_session_cache'):
        st.session_state.swan_ui_loaded_session_cache = {}

    # logo
    st.sidebar.image(f"{d}/cerberus_logo.png", width=300)

    ##### DATA LOADING
    # load data from cmd line
    load_data(cli_args.h5)
    load_swan_data(cli_args.sg)

    # # if no CerberusAnnotation was loaded, add option to upload file on the sidebar
    # upload_expander = st.sidebar.expander("Upload files", expanded=False)
    # with upload_expander:
    #
    #     uploaded_file = st.file_uploader("Upload CerberusAnnotation", type=["h5"])
    #     if not st.session_state.data_loaded:
    #         if uploaded_file is not None:
    #             import tempfile
    #             with tempfile.NamedTemporaryFile(delete=False, suffix=".h5") as tmp:
    #                 tmp.write(uploaded_file.getbuffer())
    #                 file_to_load = tmp.name
    #                 load_data(file_to_load)
    #
    #     st.markdown("---")
    #
    #
    #     # if no SwanGraph was loaded, add option to upload file on the sidebar
    #     uploaded_file = st.file_uploader("Upload SwanGraph", type=[".p"])
    #     if not st.session_state.data_loaded:
    #         if uploaded_file is not None:
    #             import tempfile
    #             with tempfile.NamedTemporaryFile(delete=False, suffix=".p") as tmp:
    #                 tmp.write(uploaded_file.getbuffer())
    #                 swan_file_to_load = tmp.name
    #                 load_swan_data(swan_file_to_load)

    # creating tabs for each visualization
    tab_simplex_view, tab_swan_view, = st.tabs(["Simplex View", "Transcript structure / expression view"])

    from tabs.simplex_tab import render_simplex_tab
    from tabs.swan_tab import render_swan_tab

    with tab_simplex_view:
        render_simplex_tab()

    with tab_swan_view:
        render_swan_tab(st.session_state.sg)

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
