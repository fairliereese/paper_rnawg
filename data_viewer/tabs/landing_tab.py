import streamlit as st
import cerberus
import matplotlib.pyplot as plt
from utils import *

def render_landing_tab():
    """
    Choose between human and mouse
    """

    st.markdown('## Welcome to the ENCODE4 long-read RNA-seq data viewer')

    st.markdown('#### Select species to begin')
    species = st.selectbox(
        label='Species',
        options=['Human', 'Mouse']
    )

    st.species = species
