import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import re
from adjustText import adjust_text
import textwrap as tw
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
sys.path.append('../')
from ad_viz_utility import *
from matplotlib.lines import Line2D

# Page Config
set_config()

# # --- SIDEBAR FILTERS --- 
st.sidebar.header("Controls")

cell_type = st.sidebar.selectbox("Select Cell Type", ["astrocytes", "microglia", "pericytes"])
if cell_type == "astrocytes":
    sampcol = 'Samples'
elif cell_type == "microglia":
    sampcol = 'Sample'
elif cell_type == "pericytes":
    sampcol = 'Sample'

intensity, samp_meta, entry_meta = load_data(
f'data/proteomics/{cell_type}_prot_data.csv',
f'data/proteomics/{cell_type}_prot_samp_meta.csv',
f'data/proteomics/{cell_type}_prot_meta.csv', 
sample_col = sampcol)

comparison, enrichment, database = sidebar_selectors()

comp_key, csv_ending = file_keys(comparison, enrichment)

selection = {
    "cell_type": cell_type,
    "comparison": comparison,
    "enrichment": enrichment,
    "database": database
}

enrichment_df = load_enrichment_df(cell_type, comp_key, enrichment, database, selection)
enrichment_df["Description"] = "(" + enrichment_df["Database"] + ") " + enrichment_df["Description"]

filtered_enrichment_df = pathway_filters(enrichment_df, enrichment)

# selecting a pathway
dropdown_df = filtered_enrichment_df.sort_values(by = 'p.adjust').drop_duplicates(subset='Description', keep='first')
pathway_list = dropdown_df['Description'].unique()
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select or Search Pathways", pathway_list)

selected_db = filtered_enrichment_df[filtered_enrichment_df['Description'] == pathway_name]['Database'].iloc[0]
if selected_db == "GO":
    pathway_identifier, gene_set, protein_identifier = "TERM", load_GO_data(cell_type), "UNIPROT"
elif selected_db == "reactome":
    pathway_identifier, gene_set, protein_identifier = "clean_name", filtered_reactome_gene_sets, "gene"

# Setting up pages
clean_pathway_name = pathway_name.replace("(GO) ", "").replace("(reactome) ", "")
id_list = pathway_to_uniprot(gene_set, selected_db, pathway_identifier, clean_pathway_name, protein_identifier, entry_meta)

# subset the dataframe to the proteins of interest
data_sub = intensity.T.loc[id_list]

# melt the dataframe
df = data_sub.reset_index(names = ['uniprot_id']).melt(id_vars = 'uniprot_id',
                                                       value_name = 'protein_level',
                                                       var_name = 'Sample_Rep')
df['Sample'] = df['Sample_Rep'].apply(lambda x: samp_meta.loc[x, sampcol])

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

title_text = ["<br>".join(
    tw.wrap(pathway_name, width=40)
)]

font_size = get_font_size(title_text)

with st.spinner("Loading pathway boxplots..."):
    
    pathway_info = enrichment_df[enrichment_df['Description'] == pathway_name].copy() 
    possible_summary_cols = [
        'Description',
        'p.adjust',
        'setSize',
        'NES',
        'enrichment_type']
    summary_table = pathway_info[[x for x in possible_summary_cols if x in pathway_info.columns]]
    st.dataframe(summary_table)

# --- VISUALIZATION ---
fig = make_subplots(
    rows=1,
    cols=1,
    subplot_titles=title_text
)

base = px.scatter(
    wide,
    x='AG09173_APOE3',
    y=f'AG09173_{gene_var}'
)

base.data[0].name = "AG09"
base.data[0].showlegend = True

fig.add_traces(base.data)

fig.update_traces(
    marker=dict(size=5, color="orange", opacity=0.5),
    customdata=np.stack([
        wide['AG09173_APOE3'],
        wide[f'AG09173_{gene_var}']
    ], axis=-1),
    hovertemplate=(
        "Protein level in APOE3 (AG09): %{customdata[0]:.3f}<br>"
        f"Protein level in {gene_var} (AG09): %{{customdata[1]:.3f}}<extra></extra>"
    )
)

fig.add_scatter(
    x=wide['AG10788_APOE3'],
    y=wide[f'AG10788_{gene_var}'],
    mode='markers',
    marker=dict(size=5, color="blue", opacity=0.5),
    name="AG10",
    customdata=np.stack([
        wide['AG10788_APOE3'],
        wide[f'AG10788_{gene_var}']
    ], axis=-1),
    hovertemplate=(
        "Protein level in APOE3 (AG10): %{customdata[0]:.3f}<br>"
        f"Protein level in {gene_var} (AG10): %{{customdata[1]:.3f}}<extra></extra>"
    )
)

data_lim = wide.max().max()
lims = (-1*data_lim-0.2,data_lim+0.2)

fig.add_scatter(
    x=lims,
    y=lims,
    mode='lines',
    line=dict(color='grey', dash='dash'),
    name='y = x',
    showlegend=False
)

fig.update_layout(
    xaxis_title=dict(
        text="Normalized protein level<br>in APOE3 sample",
        font=dict(size=27, color="black")
    ),
    yaxis_title=dict(
        text=f"Normalized protein level<br>in {gene_var} sample",
        font=dict(size=27, color="black")
    ),
    xaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
    yaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
    legend=dict(
        font=dict(size=21)
    ),
    height=600,
    width=700,
)

fig.update_annotations(
    font=dict(size=font_size, color="black"),
    align="center"
)

fig.update_xaxes(
    showline=True,
    linecolor="black",
    tickcolor="black",
    tickfont=dict(color="black", size=21),
    zeroline=False)

fig.update_yaxes(
    showline=True,
    linecolor="black",
    tickcolor="black",
    tickfont=dict(color="black", size=21),
    zeroline=False)

st.plotly_chart(fig, use_container_width=True)