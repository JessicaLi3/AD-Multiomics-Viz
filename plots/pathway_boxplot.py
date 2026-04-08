import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw
import tarfile
import pdb
import math

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")

# Function that goes from pathway name to a list of uniprot IDs

filtered_reactome_gene_sets = pd.read_csv('data/gene_sets/filtered_reactome_gene_sets.csv', index_col=0)
filtered_GO_gene_sets = pd.read_csv('data/gene_sets/filtered_GO_gene_sets.csv', index_col=0)
intensity = pd.read_csv('data/proteomics/astrocytes_prot_data.csv', index_col=0)
sample_metadata = pd.read_csv('data/proteomics/astrocytes_prot_samp_meta.csv', index_col=0)
entry_meta = pd.read_csv('data/proteomics/astrocytes_prot_meta.csv', index_col=0) 

def read_targz_to_pandas(file, tar_path):
    with tarfile.open(tar_path, "r:gz") as tar:
        file_obj = tar.extractfile(file)
        dataframe = pd.read_csv(file_obj)
    return dataframe

def pathway_to_uniprot(database, db_name, pathway_identifier, pathway_name, protein_identifier):
    if db_name == 'GO':
        filtered = database[database[pathway_identifier] == pathway_name]
        proteins = list(filtered[protein_identifier])
    elif db_name == 'reactome':
        filtered = database[database[pathway_identifier] == pathway_name]
        genes = filtered[protein_identifier]  
        match = entry_meta[entry_meta['id.mapped'].isin(genes)]
        proteins = list(match.index.values)
    return proteins

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

# # --- SIDEBAR FILTERS --- gene ontology or reactome
st.sidebar.header("Controls")
comparison = st.sidebar.selectbox("Select Comparison", ["ABCA7LOF vs APOE3", "APOE4 vs APOE3"])
enrichment = st.sidebar.selectbox("Select Enrichment Type", ["GSEA", "ORA"])
database = st.sidebar.selectbox("Select Database", ["reactome", "GO"])

if comparison == "ABCA7LOF vs APOE3":
    comp_key = "ABCA7LOFvnormal"
elif comparison == "APOE4 vs APOE3":
    comp_key = "APOE4v3"

if enrichment == "GSEA":
    csv_ending = "logFC_pval_table"
elif enrichment == "ORA":
    csv_ending = "enrichment_table"

enrichment_df = pd.read_csv(f'data/enrichment/astrocytes_proteome/{comp_key}_{enrichment}_{database}_{csv_ending}.csv')
pathway_list = enrichment_df['Description']

pathway_name = st.sidebar.selectbox("Select or Search Pathways", pathway_list)

if database == "reactome":
    pathway_identifier = "clean_name"
    gene_set = filtered_reactome_gene_sets
    protein_identifier = "gene"
elif database == "GO":
    pathway_identifier = "TERM"
    gene_set = filtered_GO_gene_sets
    protein_identifier = "UNIPROT"

# Setting up pages
id_list = pathway_to_uniprot(gene_set, database, pathway_identifier, pathway_name, protein_identifier)
plots_per_page = 9

if 'current_page' not in st.session_state:
    st.session_state.current_page = 0

total_pages = math.ceil(len(id_list) / plots_per_page)
start_idx = st.session_state.current_page * plots_per_page
end_idx = start_idx + plots_per_page
current_proteins = id_list[start_idx:end_idx]

st.write(f"### Showing proteins {start_idx + 1} to {min(end_idx, len(id_list))}")

# --- PLOTS ---
fig, axs = plt.subplots(3,3, figsize = (24,18), gridspec_kw = {'hspace':0.7,'wspace':0.4})
axs = axs.flatten()

for i, prot_id in enumerate(current_proteins):
# boxplots
    df_plot = intensity.loc[[prot_id]].T
    df_plot = df_plot.merge(sample_metadata[['Parental.line', 'CRISPR.edits', 'Samples']], left_index = True, right_index = True)
    ax = axs[i]
    ax = sns.boxplot(
        data = df_plot,
        x = 'Parental.line',
        y = prot_id,
        hue = 'CRISPR.edits', ax=ax)
    
# plot the stripplot
    ax = sns.stripplot(
        data = df_plot,
        x = 'Parental.line',
        y = prot_id,
        hue = 'CRISPR.edits',
        ax=ax,
        dodge = True,
        color = 'black',
        legend = False)

    prot_name = lookup_by_protid(prot_id, entry_meta, 'id.mapped')
    prot_desc = lookup_by_protid(prot_id, entry_meta, 'entry_name')
    prot_id = lookup_by_protid(prot_id, entry_meta, 'id.description')
    ax.set_title(tw.fill(f'{prot_name} ({prot_id})\n{prot_desc}',
                         width = 38,
                         replace_whitespace = False,
                         drop_whitespace=False), fontsize = 16)

    ax.set_ylabel("Normalized levels", fontsize=13)
    ax.set_xlabel("", fontsize=13)
    ax.tick_params(axis = 'y', direction = 'inout', length = 8, width = 3, labelsize = 13)
    ax.tick_params(axis = 'x', direction = 'inout', length = 8, width = 3, labelsize = 13)
    ax.legend(bbox_to_anchor=(1,1), loc= 'upper left', fontsize=16)
    if i!=2:
        ax.legend().remove()
    for j in range(len(current_proteins), len(axs)):
        axs[j].axis('off')

st.pyplot(fig)

# Page controls
col1, col2, col3 = st.columns([1,2,1])

with col1:
    if st.button("Previous") and st.session_state.current_page > 0:
        st.session_state.current_page -= 1
        st.rerun()

with col2:
    target_page = st.number_input(
        "Jump to page",
        min_value=1,
        max_value=total_pages,
        value=st.session_state.current_page + 1,
        step=1,
        label_visibility="collapsed"
    )
    if target_page-1 != st.session_state.current_page:
        st.session_state.current_page = target_page-1
        st.rerun()

with col3:
    if st.button("Next") and st.session_state.current_page < total_pages-1:
        st.session_state.current_page += 1
        st.rerun()
