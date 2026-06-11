import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import re
import io
from adjustText import adjust_text
import textwrap as tw
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
sys.path.append('../')
from ad_viz_utility import *

# Page Config
set_config()

# # --- SIDEBAR FILTERS --- gene ontology or reactome
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

diff = check_read_csv(f'data/differential_analysis/{cell_type}_proteome/{comp_key}_covDonor_genotype.csv', index_col = 0, selection=selection)

# sliders
st.sidebar.header("Pathway Filters") 

enrichment_df = load_enrichment_df(cell_type, comp_key, enrichment, database, selection)
enrichment_df["Description"] = "(" + enrichment_df["Database"] + ") " + enrichment_df["Description"]
    
filtered_enrichment_df = pathway_filters(enrichment_df, enrichment)
    
# selecting a pathway 
dropdown_df = filtered_enrichment_df.sort_values('p.adjust').drop_duplicates(subset='Description', keep='first')
pathway_list = dropdown_df['Description'].unique()
pathway_list = sorted(pathway_list)
pathway_name = st.sidebar.selectbox("Select Pathway", pathway_list)

selected_db = filtered_enrichment_df[filtered_enrichment_df['Description'] == pathway_name]['Database'].iloc[0]
if selected_db == "GO":
    pathway_identifier, gene_set, protein_identifier = "TERM", load_GO_data(cell_type), "UNIPROT"
elif selected_db == "reactome":
    pathway_identifier, gene_set, protein_identifier = "clean_name", filtered_reactome_gene_sets, "gene"
    
clean_pathway_name = pathway_name.replace("(GO) ", "").replace("(reactome) ", "")
id_list = pathway_to_uniprot(gene_set, selected_db, pathway_identifier, clean_pathway_name, protein_identifier, entry_meta)

# adding some adjustments
diff['id.mapped'] = diff.index.map(entry_meta['id.mapped'])
diff['entry_name'] = diff.index.map(entry_meta['entry_name'])
diff['log_adjpval'] = diff['adj.P.Val'].apply(lambda x: np.log10(x) * -1) 
diff['logFC_logpval'] = np.abs(diff['logFC']*diff['log_adjpval']) 

ranked_id_list = (diff.loc[id_list].sort_values(by='logFC_logpval', ascending=False).index.tolist())
top20 = ranked_id_list[:20]
top20_df = diff.loc[top20]
pathway = diff.loc[ranked_id_list]

title_text = ["<br>".join(
    tw.wrap(pathway_name, width=40)
)]

font_size = get_font_size(title_text)

with st.spinner("Loading pathway boxplots..."):
    
    pathway_info = filtered_enrichment_df[filtered_enrichment_df['Description'] == pathway_name].copy() 
    possible_summary_cols = [
        'Description',
        'p.adjust',
        'setSize',
        'NES',
        'enrichment_type']
    summary_table = pathway_info[[x for x in possible_summary_cols if x in pathway_info.columns]]
    st.dataframe(summary_table)

col1, col2 = st.columns([2, 1])
with col1:
    options = ["Presentation Ready", "Interactive"]
    view = st.segmented_control("Visualization", options, selection_mode="single", default="Presentation Ready")

with col2: 
    show_labels = st.checkbox("Show top 20 protein labels", value=True)

if view == "Presentation Ready": 

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
    if show_labels:
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

elif view == "Interactive":
    fig = make_subplots(
        rows=1,
        cols=1,
        subplot_titles=title_text
    )
    
    base = px.scatter(
        diff,
        x='logFC',
        y='log_adjpval'
    )
    
    base.update_traces(
        selector=dict(mode='markers'),
        marker=dict(
            size=2,
            color="lightgrey",
            opacity=0.3,
            line=dict(width=0)
        ),
        hoverinfo="none",
        hovertemplate=None
    )
        
    fig.add_traces(base.data)
    
    top_scatter = px.scatter(
        diff.loc[ranked_id_list],
        x='logFC',
        y='log_adjpval'
    )
    
    top_scatter.update_traces(
        selector=dict(mode='markers'),
        marker=dict(
            size=4,           
            color="teal",
            opacity = 1,
            line=dict(color="teal", width=1)
        ),
        customdata=np.stack([
            pathway['id.mapped'],
            pathway['entry_name'],
            pathway.index,
        ], axis=-1),
        hovertemplate=(
            "Gene Symbol: %{customdata[0]}<br>"
            "Full name: %{customdata[1]}<br>"
            "UniProt: %{customdata[2]}<br>"
            "logFC: %{x:.2f}<br>"
            "adj p-value: %{y:.2f}<extra></extra>"
        )
    )
    
    fig.add_traces(top_scatter.data)
    fig.update_layout(template="plotly_white")
    
    if show_labels:
        for prot_id in top20:
            fig.add_annotation(
                x=diff.loc[prot_id, 'logFC'],
                y=diff.loc[prot_id, 'log_adjpval'],
                text=diff.loc[prot_id, 'id.mapped'],
                showarrow=False,
                xshift=10,
                yshift=10,
                align="left",
                font=dict(size=13, color="black")
            )
        
    fig.update_layout(
        xaxis_title=dict(
            text="logFC",
            font=dict(size=27, color="black")),
        yaxis_title=dict(
            text="log10(adj. pval)",
            font=dict(size=27, color="black")),
        plot_bgcolor='white',
        xaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
        yaxis=dict(showgrid=False, gridcolor='black', showline=True, zeroline=False),
        height=600,
        width=700,
        showlegend=False
    )
    
    fig.layout.annotations[0].update(
        font=dict(size=font_size, color="black")
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
    

