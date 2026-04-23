import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap as tw
import tarfile
import pdb
import math
import sys
sys.path.append('../')
from ad_viz_utility import *

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="wide")

# # --- SIDEBAR FILTERS --- 
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

comparison = st.sidebar.selectbox("Select Comparison", ["ABCA7LOF vs APOE3", "APOE4 vs APOE3"])
enrichment = st.sidebar.selectbox("Select Enrichment Type", ["GSEA", "ORA"])
database = st.sidebar.selectbox("Select Database", ["reactome", "GO"])
 
comp_key, csv_ending, pathway_identifier, gene_set, protein_identifier = file_keys(comparison, enrichment, database, cell_type)

diff = pd.read_csv(f'data/differential_analysis/{cell_type}_proteome/{comp_key}_covDonor_genotype.csv', index_col = 0)
enrichment_df = pd.read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_{database}_{csv_ending}.csv')

#sorting by significance 
diff['log_adjpval'] = diff['adj.P.Val'].apply(lambda x: np.log10(x) * -1)
diff['logFC_logpval'] = np.abs(diff['logFC']*diff['log_adjpval']) # rank proteins by this

# selecting a pathway
pathway_list = enrichment_df['Description']
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select or Search Pathways", pathway_list)

# Setting up pages
id_list = pathway_to_uniprot(gene_set, database, pathway_identifier, pathway_name, protein_identifier, entry_meta)
id_list = [p for p in id_list if p in entry_meta.index]
ranked_id_list = (diff.loc[id_list].sort_values(by='logFC_logpval', ascending=False).index.tolist())

plots_per_page = 9

if 'current_page' not in st.session_state:
    st.session_state.current_page = 0

total_pages = math.ceil(len(ranked_id_list) / plots_per_page)
start_idx = st.session_state.current_page * plots_per_page
end_idx = start_idx + plots_per_page
current_proteins = ranked_id_list[start_idx:end_idx]
valid_proteins = [p for p in current_proteins if p in intensity.columns]

st.write(f"### Showing proteins {start_idx + 1} to {min(end_idx, len(ranked_id_list))} "
        f"out of {len(ranked_id_list)} proteins (Page {st.session_state.current_page + 1} of {total_pages})")

with st.spinner("Loading pathway boxplots..."):
    
    pathway_info = enrichment_df[enrichment_df['Description'] == pathway_name].copy()
    if 'setSize' not in pathway_info.columns:
        pathway_info['setSize'] = np.nan
    if 'NES' not in pathway_info.columns:
        pathway_info['NES'] = np.nan

    summary_table = pathway_info[[
        'Description',
        'p.adjust',
        'setSize',
        'NES']]
    st.dataframe(summary_table)
    
    # --- PLOTS ---
    fig, axs = plt.subplots(3,3, figsize = (24,18), gridspec_kw = {'hspace':0.7,'wspace':0.4})
    axs = axs.flatten()

    for i, prot_id in enumerate(valid_proteins):
    # boxplots
        df_plot = intensity[[prot_id]].copy()
        df_plot = df_plot.merge(samp_meta[['Parental.line', 'CRISPR.edits', sampcol]], left_index = True, right_index = True)
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
