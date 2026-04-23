import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw
import sys
sys.path.append('../')
from ad_viz_utility import *

max_msg_size = st.get_option("server.maxMessageSize")

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")
# add a comment

st.title("🔬 Proteomics Data Explorer v0.1")

# Load Data
st.sidebar.header("Controls")
cell_type = st.sidebar.selectbox("Select Cell Type", ["astrocytes", "microglia", "pericytes"])
if cell_type == "astrocytes":
    sampcol = 'Samples'
elif cell_type == "microglia":
    sampcol = 'Sample'
elif cell_type == "pericytes":
    sampcol = 'Sample'
with st.spinner("Loading 10,000 proteins..."):
    intensity, samp_meta, entry_meta = load_data(
        f'data/proteomics/{cell_type}_prot_data.csv',
        f'data/proteomics/{cell_type}_prot_samp_meta.csv',
        f'data/proteomics/{cell_type}_prot_meta.csv', 
        sample_col = sampcol)

# --- SIDEBAR FILTERS ---
# Exclude metadata columns so they don't show up in the protein dropdown
meta_cols = ['Parental.line', 'CRISPR.edits', 'Samples']
entry_meta['uniprot_id.mapped_concat'] = entry_meta['id.mapped'] + '_' + entry_meta.index
protein_list = entry_meta['uniprot_id.mapped_concat']
protein_list = sorted(protein_list)

st.sidebar.header("Controls")
id_name = st.sidebar.selectbox("Select or Search Protein", protein_list)
prot_id = entry_meta[entry_meta['uniprot_id.mapped_concat'] == id_name].index[0]

# --- VISUALIZATION ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader(f"Raw Values: {prot_id}")
    protein_data = intensity[[prot_id]].copy()
    protein_data.columns = ["Normalized Abundance"]
    st.dataframe(protein_data)

with col2:
    st.subheader("Expression Boxplot")
    fig_length = 8
    fig_height = 5
    fig, ax = plt.subplots(figsize=(fig_length, fig_height))    
    # boxplots
    ax = sns.boxplot(
        data = intensity,
        x = 'Parental.line',
        y = prot_id,
        hue = 'CRISPR.edits', ax=ax)

    # plot the stripplot
    ax = sns.stripplot(
        data = intensity,
        x = 'Parental.line',
        y = prot_id,
        hue = 'CRISPR.edits',
        ax=ax,
        dodge = True,
        color = 'black',
        legend = False)
    ax.legend(bbox_to_anchor=(1,1), loc= 'upper left', fontsize = 25)
    
    # adds the title
    prot_name = lookup_by_protid(prot_id, entry_meta, 'id.mapped')
    ax.set_title(prot_name)
    prot_desc = lookup_by_protid(prot_id, entry_meta, 'entry_name')
    ax.set_title(tw.fill(f'{prot_name}\n{prot_desc}',
                         width = 32,
                         replace_whitespace = False,
                         drop_whitespace=False), fontsize = 25)
    ax.set_ylabel('normalized levels', fontsize = 25)
    ax.set_xlabel('')
 
    # changes font size of axis ticks
    ax.tick_params(axis = 'y', direction = 'inout', length = 5, width = 1, labelsize = 15)
    ax.tick_params(axis = 'x', direction = 'inout', length = 5, width = 1, labelsize = 20)

    st.pyplot(fig)