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
from matplotlib.lines import Line2D

# Page Config
st.set_page_config(page_title="Proteomics v1.0", layout="centered")

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

enrichment_df = pd.read_csv(f'data/enrichment/{cell_type}_proteome/{comp_key}_{enrichment}_{database}_{csv_ending}.csv')

# selecting a pathway
enrichment_df = enrichment_df.sort_values(by = 'p.adjust')
pathway_list = enrichment_df['Description']
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select or Search Pathways", pathway_list)

# Setting up pages
id_list = pathway_to_uniprot(gene_set, database, pathway_identifier, pathway_name, protein_identifier, entry_meta)
id_list = [p for p in id_list if p in entry_meta.index]

# subset the dataframe to the proteins of interest
data_sub = intensity.T.loc[id_list]

# melt the dataframe
df = data_sub.reset_index(names = ['uniprot_id']).melt(id_vars = 'uniprot_id',
                                                       value_name = 'protein_level',
                                                       var_name = 'Sample_Rep')
df['Sample'] = df['Sample_Rep'].apply(lambda x: samp_meta.loc[x,'Samples'])

# take the mean for every sample type
mean_df = df.groupby(['Sample', 'uniprot_id'])[['protein_level']].agg(np.mean)

# pivot the data
wide = mean_df.reset_index(level=[0,1], names = ['Sample','uniprot_id']).pivot(index='Sample', columns = 'uniprot_id').T
wide.index = wide.index.droplevel(0)
wide.columns = [c.replace('.', '_') for c in wide.columns]

if comparison == "ABCA7LOF vs APOE3":
        gene_var = "ABCA7_G2"
elif comparison == "APOE4 vs APOE3":
        gene_var = "APOE4"

# --- VISUALIZATION ---
fig, ax = plt.subplots(figsize=(11,9))

ax = sns.scatterplot(data = wide, x = 'AG09173_APOE3', y = f'AG09173_{gene_var}',
                    s = 10, alpha = 0.5, linewidth = 0, label = 'AG09')
ax = sns.scatterplot(data = wide, x = 'AG10788_APOE3', y = f'AG10788_{gene_var}', s = 10, alpha = 0.5, linewidth = 0,
                    label = 'AG10')
ax.legend(bbox_to_anchor = (1,1), fontsize = 15)
data_lim = wide.max().max()
lims = (-1*data_lim-0.2,data_lim+0.2)
ax.set_xlim(lims)
ax.set_ylim(lims)
identity_line = Line2D(lims, lims, color="grey", linestyle="--", alpha=0.75)
ax.add_line(identity_line)

ax.set_xlabel('Normalized protein level\nin APOE3 sample', fontsize = 15)
ax.set_ylabel(f'Normalized protein level\nin {gene_var} sample', fontsize = 15)
ax.set_title(pathway_name, fontsize = 25)

st.pyplot(fig)