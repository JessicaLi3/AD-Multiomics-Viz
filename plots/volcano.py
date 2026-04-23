import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import re
from adjustText import adjust_text
import textwrap as tw
import sys
sys.path.append('../')
from ad_viz_utility import *

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="centered")

# # --- SIDEBAR FILTERS --- gene ontology or reactome
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
gsea_df = pd.read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_{database}_{csv_ending}.csv')

# selecting a pathway 
pathway_list = gsea_df['Description']
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select Pathway", pathway_list)
id_list = pathway_to_uniprot(gene_set, database, pathway_identifier, pathway_name, protein_identifier, entry_meta)
id_list = [p for p in id_list if p in diff.index]

# adding some adjustments
diff['id.mapped'] = diff.index.map(entry_meta['id.mapped'])
diff['entry_name'] = diff.index.map(entry_meta['entry_name'])
diff['log_adjpval'] = diff['adj.P.Val'].apply(lambda x: np.log10(x) * -1) # plot this in volcano
diff['logFC_logpval'] = np.abs(diff['logFC']*diff['log_adjpval']) # rank proteins by this

ranked_id_list = (diff.loc[id_list].sort_values(by='logFC_logpval', ascending=False).index.tolist())

# --- VISUALIZATION ---
fig, ax = plt.subplots(figsize=(9,7))

# plot all the proteins in the experiment
ax = sns.scatterplot(data = diff, x = 'logFC', y = 'log_adjpval',
                     s = 20, alpha = 0.4, linewidth = 0, color = 'lightgrey')

# plot all the proteins in the pathway
ax = sns.scatterplot(data = diff.loc[ranked_id_list], x = 'logFC', y = 'log_adjpval',
                     s = 20, alpha = 1, linewidth = 0, color = 'mediumaquamarine')
ax.set_title(tw.fill(f'{pathway_name}', width = 50), fontsize = 25)
ax.set_xlabel(f'logFC\n({comp_key})', fontsize = 15)
ax.set_ylabel(f'log10(adj. pval)', fontsize = 15)

# label the top 20 proteins in the pathway (ranked by logFC_logpval)
annotations = []
for prot_id in ranked_id_list[0:20]:
    x = diff.loc[prot_id, 'logFC']
    y = diff.loc[prot_id, 'log_adjpval']
    text = diff.loc[prot_id, 'id.mapped']
    annotations.append(ax.annotate(text, (x,y), fontsize = 12))

# moves the labels around
adjust_text(
    annotations,
    force_text = (0.5,0.5),
    force_static = (0.01,0.01),
    max_move = (30,20),
    force_explode = (0.5,0.5),
    arrowprops = dict(arrowstyle = '-'),
    ax = ax);

st.pyplot(fig)
