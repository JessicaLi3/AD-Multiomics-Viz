import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw

max_msg_size = st.get_option("server.maxMessageSize")
st.write(f"Current Max Message Size: {max_msg_size} MB")

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")
# add a comment
# --- DATA LOADING (CACHED) ---
@st.cache_data
def load_data(data_path, samp_meta_path, entry_meta_path, sample_col):
    # load CSVs
    data = pd.read_csv(data_path, index_col=0)
    samp_meta = pd.read_csv(samp_meta_path, index_col=0)
    entry_meta = pd.read_csv(entry_meta_path, index_col=0)
    
    ## formatting
    # transposes the table so that samples are rows
    data = data.T

    # merge with relevant sample metadata columns
    data = data.merge(samp_meta[[
        'Parental.line', 'CRISPR.edits', sample_col]],
        left_index = True, right_index = True)

    return data, samp_meta, entry_meta
    
def lookup_by_protid(protein_id: str, protein_metadata: pd.DataFrame,
                     name_col: str):
    '''
    Uses a protein ID to look up the rows of the `protein_metadata` table
    where that ID matches the index. Returns the contents of `name_col`
    for that row.
    For proteomics, this is the uniprot ID.
    Args:
    - protein_id (str): The protein identifier to look up in the index
    - protein_metadata (DataFrame): the protein metadata dataframe
    - name_col (str): Column name of the identifier/descriptor to look up
    Returns:
    - str
    '''

    # returns np.nan if the protein_id is not found.
    if protein_id not in protein_metadata.index:
        return np.nan
    # otherwise returns the contents of name_col for that row
    # assumes only one row will be returned (index is unique)
    return protein_metadata.loc[protein_id, name_col]
    
st.title("🔬 Proteomics Data Explorer v0.1")

# Load Data
st.sidebar.header("Controls")
cell_type = st.sidebar.selectbox("Select Cell Type", ["astrocytes", "microglia"])
if cell_type == "astrocytes":
    sampcol = 'Samples'
elif cell_type == "microglia":
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